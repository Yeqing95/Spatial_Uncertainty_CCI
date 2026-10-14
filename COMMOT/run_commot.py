"""Small COMMOT checks; edit CONFIG first. Each stage runs in a fresh process.
Usage: python run_commot.py --stage all
Stages: official_toy, z_test, fixed3d, official_basic, all.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CONFIG = {
    "input_dir": str(ROOT / "example_inputs"),
    "output_dir": str(ROOT / "outputs"),
    "seed": 888,
    "normalization_target": 10000.0,
    "fixed_cutoff_mode": "per_lr",  # per_lr reads declared radius_um; or common
    "common_cutoff_um": 36.0,
    "cot_nitermax": 10000,
    "cot_eps_p": 0.1,
    "cot_rho": 10.0,
}

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def run_stage(stage):
    start = time.perf_counter()
    import resource
    import importlib.metadata as md
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    from scipy.spatial import distance_matrix
    import anndata as ad
    import scanpy as sc
    import commot as ct

    np.random.seed(CONFIG["seed"])
    out = Path(CONFIG["output_dir"]) / stage
    out.mkdir(parents=True, exist_ok=True)
    calls = []
    validations = []

    def peak_rss_mib():
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return peak / (2**20 if sys.platform == "darwin" else 1024)

    def check(condition, description):
        if not bool(condition):
            raise AssertionError(description)
        validations.append(description)

    def infer(a, lr, db, threshold, heteromeric=False, name=None):
        # Caller supplies a fresh AnnData, or explicitly controls cached distances.
        before = peak_rss_mib()
        t0 = time.perf_counter()
        try:
            ct.tl.spatial_communication(
                a, database_name=db, df_ligrec=lr, dis_thr=threshold,
                heteromeric=heteromeric, pathway_sum=True, smooth=False,
                cot_nitermax=CONFIG["cot_nitermax"],
                cot_eps_p=CONFIG["cot_eps_p"], cot_rho=CONFIG["cot_rho"],
            )
        finally:
            elapsed = time.perf_counter() - t0
        calls.append({
            "call": name or db, "inference_seconds": elapsed,
            "process_peak_RSS_before_inference_MiB": before,
            "process_peak_RSS_after_inference_MiB": peak_rss_mib(),
            "memory_scope": "OS process lifetime high-water mark, including imports/loading; not inference-only allocation",
        })

    def export(a, lr, db, ids=None):
        rows, summary = [], []
        a.obs_names.to_series(name="cell_id").to_csv(out / "cell_order.tsv", sep="\t", index=False)
        for k, (ligand, receptor, pathway) in enumerate(lr.itertuples(index=False, name=None)):
            lr_id = ids[k] if ids is not None else f"LR{k+1:03d}"
            key = f"commot-{db}-{ligand}-{receptor}"
            check(key in a.obsp, f"LR output present: {key}")
            m = sp.csr_matrix(a.obsp[key])
            m.eliminate_zeros()
            check(m.shape == (a.n_obs, a.n_obs), f"matrix shape: {lr_id}")
            check(np.isfinite(m.data).all() and (m.data >= 0).all(), f"finite nonnegative scores: {lr_id}")
            np.testing.assert_allclose(np.asarray(m.sum(axis=1)).ravel(),
                a.obsm[f"commot-{db}-sum-sender"][f"s-{ligand}-{receptor}"], atol=1e-8)
            np.testing.assert_allclose(np.asarray(m.sum(axis=0)).ravel(),
                a.obsm[f"commot-{db}-sum-receiver"][f"r-{ligand}-{receptor}"], atol=1e-8)
            validations.append(f"row=sender and column=receiver sums match: {lr_id}")
            sp.save_npz(out / f"{lr_id}.npz", m)
            coo = m.tocoo()
            frame = pd.DataFrame({
                "sender_id": a.obs_names.to_numpy()[coo.row],
                "receiver_id": a.obs_names.to_numpy()[coo.col],
                "lr_id": lr_id, "ligand": ligand, "receptor": receptor,
                "score": coo.data, "is_self": coo.row == coo.col,
            })
            xyz = np.asarray(a.obsm["spatial"])
            frame["distance_coordinate_units"] = np.linalg.norm(xyz[coo.row] - xyz[coo.col], axis=1)
            if "section_id" in a.obs:
                sec = a.obs["section_id"].to_numpy()
                frame["sender_section"] = sec[coo.row]
                frame["receiver_section"] = sec[coo.col]
                frame["cross_section"] = sec[coo.row] != sec[coo.col]
            rows.append(frame)
            summary.append({"lr_id": lr_id, "ligand": ligand, "receptor": receptor,
                "n_nonzero_edges": int(m.nnz), "n_self_edges": int((coo.row == coo.col).sum()),
                "n_cross_section_edges": int(frame["cross_section"].sum()) if "cross_section" in frame else None})
        edges = pd.concat(rows, ignore_index=True)
        edges.to_csv(out / "lr_edges.tsv.gz", sep="\t", index=False)
        pd.DataFrame(summary).to_csv(out / "lr_summary.tsv", sep="\t", index=False)
        return edges

    if stage == "official_toy":
        # The 4-cell example in the official spatial_communication API documentation.
        x = np.array([[1,0,0,0], [2,0,0,0], [0,1,1,4], [0,3,4,1]], dtype=float)
        a = ad.AnnData(sp.csr_matrix(x), obs=pd.DataFrame(index=["c0","c1","c2","c3"]),
            var=pd.DataFrame(index=["LigA","RecA","RecB","RecC"]))
        a.obsm["spatial"] = np.array([[0,0], [1,0], [0,1], [1,2]], dtype=float)
        lr = pd.DataFrame([["LigA","RecA_RecB","PathwayA"], ["LigA","RecC","PathwayA"]])
        infer(a, lr, "toyDB", 1.1, heteromeric=True)
        e = export(a, lr, "toyDB")
        check(len(e) == 2 and set(zip(e.sender_id, e.receiver_id)) == {("c0","c2")}, "official toy edge IDs")
        np.testing.assert_allclose([a.obsp["commot-toyDB-LigA-RecA_RecB"][0,2],
                                   a.obsp["commot-toyDB-LigA-RecC"][0,2]], [0.6,0.9], rtol=1e-5)
        validations.append("official toy published scores 0.6 and 0.9 reproduced")

    elif stage == "z_test":
        lr = pd.DataFrame([["LigA", "RecA", "test_pathway"]])
        xy = np.array([[0,0], [10,0], [0,0.5], [10,0.5]], dtype=float)
        xyz0 = np.column_stack([xy, np.zeros(4)])
        xyz = xyz0.copy()
        xyz[2,2] = 5.0  # Move only receiver r1 along z. The s2->r2 pair stays available.
        results, matrices = [], {}
        scenarios = {"xy": (xy, None), "xyz_zero": (xyz0, None),
            "xyz_raised": (xyz, None),
            "xyz_explicit_3d_distance": (xyz, distance_matrix(xyz, xyz)),
            "xyz_stale_2d_distance": (xyz, distance_matrix(xy, xy))}
        for name, (coord, distances) in scenarios.items():
            a = ad.AnnData(sp.csr_matrix([[5.,0.], [5.,0.], [0.,5.], [0.,5.]]),
                obs=pd.DataFrame(index=["s1","s2","r1","r2"]), var=pd.DataFrame(index=["LigA","RecA"]))
            a.obsm["spatial"] = coord.copy()
            if distances is not None:
                a.obsp["spatial_distance"] = distances
            infer(a, lr, "zDB", 1.0, name=name)
            m = a.obsp["commot-zDB-LigA-RecA"].toarray()
            matrices[name] = m
            results.append({"scenario": name, "s1_to_r1": float(m[0,2]), "s2_to_r2": float(m[1,3])})
        np.testing.assert_allclose(matrices["xy"], matrices["xyz_zero"], rtol=1e-7, atol=1e-9)
        np.testing.assert_allclose(matrices["xyz_raised"], matrices["xyz_explicit_3d_distance"], rtol=1e-7, atol=1e-9)
        np.testing.assert_allclose(matrices["xy"], matrices["xyz_stale_2d_distance"], rtol=1e-7, atol=1e-9)
        check(matrices["xy"][0,2] > 1e-6 and matrices["xyz_raised"][0,2] == 0, "changing z removes only the targeted admissible pair")
        check(matrices["xyz_raised"][1,3] > 1e-6, "positive control remains nonzero")
        validations.extend(["zero-z agrees with 2D", "implicit XYZ agrees with explicit 3D distances", "cached 2D distances override XYZ"])
        pd.DataFrame(results).to_csv(out / "z_dimension_test.tsv", sep="\t", index=False)

    elif stage == "fixed3d":
        src = Path(CONFIG["input_dir"])
        names = ["fixed_input/counts_raw.tsv.gz", "fixed_input/coordinates.tsv", "ligand_receptor_table.tsv"]
        dump(out / "input_manifest.json", [{"path": n, "sha256": hashlib.sha256((src/n).read_bytes()).hexdigest(),
            "size_bytes": (src/n).stat().st_size} for n in names])
        counts = pd.read_csv(src / names[0], sep="\t", index_col="cell_id")
        coords = pd.read_csv(src / names[1], sep="\t", index_col="cell_id")
        lr_meta = pd.read_csv(src / names[2], sep="\t")
        check(counts.index.is_unique and coords.index.is_unique and counts.columns.is_unique, "unique cell and gene IDs")
        check(set(counts.index) == set(coords.index), "count and coordinate cell sets match")
        coords = coords.loc[counts.index].copy()
        check(set(coords.units) == {"um"}, "coordinates use declared synthetic micrometer units")
        check(coords.specimen_id.nunique() == 1 and coords.draw_id.nunique() == 1, "one specimen and complete fixed coordinate state")
        x = counts.to_numpy(dtype=float)
        check(np.isfinite(x).all() and (x >= 0).all() and np.equal(x, np.floor(x)).all(), "fixed raw counts are finite nonnegative integers")
        xyz = coords[["x","y","z"]].to_numpy(dtype=float)
        check(np.isfinite(xyz).all() and xyz.shape == (len(counts),3), "complete finite XYZ coordinates")
        required = set(lr_meta.ligand) | set(lr_meta.receptor)
        missing = sorted(required - set(counts.columns))
        dump(out / "missing_lr_genes.json", missing)
        check(not missing, "all synthetic LR genes found")
        check(not lr_meta[["ligand","receptor"]].duplicated().any(), "unique LR pairs")
        a = ad.AnnData(sp.csr_matrix(x), obs=coords.copy(), var=pd.DataFrame(index=counts.columns))
        a.layers["counts"] = a.X.copy()
        a.obsm["spatial"] = xyz
        a.uns["coordinate_units"] = "um"
        sc.pp.normalize_total(a, target_sum=CONFIG["normalization_target"])
        sc.pp.log1p(a)
        # gene_103 is ONE synthetic gene; underscores must not be parsed as complexes.
        lr = pd.DataFrame(lr_meta[["ligand","receptor","lr_id"]].to_numpy())
        if CONFIG["fixed_cutoff_mode"] == "per_lr":
            # COMMOT 0.0.3 indexes cutoffs over the full ligand x receptor grid,
            # even for prohibited LR combinations. Fill these unused entries too;
            # only df_ligrec pairs are enabled in the optimization.
            cutoff = {(l,r): float(lr_meta.radius_um.max()) for l in lr_meta.ligand for r in lr_meta.receptor}
            for row in lr_meta.itertuples():
                cutoff[(row.ligand,row.receptor)] = float(row.radius_um)
        elif CONFIG["fixed_cutoff_mode"] == "common":
            cutoff = CONFIG["common_cutoff_um"]
        else:
            raise ValueError("fixed_cutoff_mode must be per_lr or common")
        check("spatial_distance" not in a.obsp, "no stale cached distances")
        infer(a, lr, "synthetic", cutoff, heteromeric=False)
        edges = export(a, lr, "synthetic", lr_meta.lr_id.to_list())
        limits = lr_meta.set_index("lr_id").radius_um
        if CONFIG["fixed_cutoff_mode"] == "common":
            limits[:] = CONFIG["common_cutoff_um"]
        check((edges.distance_coordinate_units <= edges.lr_id.map(limits) + 1e-8).all(), "every inferred edge respects its 3D distance cutoff")
        check(edges.cross_section.any(), "nonzero cross-section CCC outputs exist")
        # Preserve native self edges; they participate in COMMOT optimization.
        # Nonself truth evaluation must explicitly account for this mismatch.
        lr_meta.to_csv(out / "lr_parameters.tsv", sep="\t", index=False)
        coords.to_csv(out / "coordinates.tsv", sep="\t")
        # AnnData cannot serialize tuple-key dictionaries. Replace only the cutoff
        # metadata representation; do not modify matrices or solver results.
        a.uns["commot-synthetic-info"]["distance_threshold"] = pd.DataFrame({
            "ligand": lr_meta.ligand.to_numpy(), "receptor": lr_meta.receptor.to_numpy(),
            "distance_um": lr_meta.lr_id.map(limits).to_numpy()})
        a.write_h5ad(out / "fixed3d_commot.h5ad", compression="gzip")
        dump(out / "dataset_dimensions.json", {"cells": a.n_obs, "genes": a.n_vars,
            "sections": int(a.obs.section_id.nunique()), "LR_pairs": len(lr)})

    elif stage == "official_basic":
        # Official Basic usage: data loading, preprocessing, inference, and
        # sender/receiver scatter plot (not the later vector-field panels).
        sc.settings.datasetdir = ROOT / "download_cache"
        a = sc.datasets.visium_sge(sample_id="V1_Mouse_Brain_Sagittal_Posterior")
        a.var_names_make_unique()
        sc.pp.normalize_total(a, inplace=True)
        sc.pp.log1p(a)
        lr = pd.DataFrame([["Tgfb1","Tgfbr1_Tgfbr2","Tgfb_pathway"],
            ["Tgfb2","Tgfbr1_Tgfbr2","Tgfb_pathway"], ["Tgfb3","Tgfbr1_Tgfbr2","Tgfb_pathway"],
            ["Fgf1","Fgfr1","FGF_pathway"], ["Fgf1","Fgfr2","FGF_pathway"]])
        infer(a, lr, "user_database", 500, heteromeric=True)
        export(a, lr, "user_database")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1,2,figsize=(10,4))
        pts = a.obsm["spatial"]
        for ax, label, short, cmap in zip(axes, ["sender","receiver"], ["s","r"], ["Blues","Reds"]):
            values = a.obsm[f"commot-user_database-sum-{label}"][f"{short}-Fgf1-Fgfr1"]
            ax.scatter(pts[:,0], pts[:,1], c=values, s=5, cmap=cmap)
            ax.set_title(label.capitalize())
        fig.savefig(out / "official_Fgf1_Fgfr1.png", dpi=180, bbox_inches="tight")
        plt.close(fig)
        validations.append("official Basic usage loading, preprocessing, 5-pair inference and sender/receiver plot executed")
        dump(out / "dataset_dimensions.json", {"cells": a.n_obs, "genes": a.n_vars,
            "distance_cutoff": 500, "units": "native tutorial coordinates; physical scale not independently verified"})
    else:
        raise ValueError(stage)

    pd.DataFrame(calls).to_csv(out / "runtime_memory.tsv", sep="\t", index=False)
    peak_mib = peak_rss_mib()
    dump(out / "run_report.json", {"status": "PASS", "stage": stage,
        "checks": validations, "stage_seconds_including_imports_io": time.perf_counter()-start,
        "process_peak_RSS_MiB_including_imports_io": peak_mib,
        "peak_scope": "fresh stage process including imports, loading, inference and export; not solver-only",
        "python": sys.version, "platform": platform.platform(), "machine": platform.machine(),
        "config": CONFIG, "versions": {k: md.version(k) for k in ["commot","numpy","scipy","pandas","anndata","scanpy","karateclub","psutil"]}})
    print(stage, "PASS", json.dumps(calls), "process peak MiB", round(peak_mib,2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["all","official_toy","z_test","fixed3d","official_basic"], default="all")
    parser.add_argument("--input-dir")
    parser.add_argument("--output-dir")
    args = parser.parse_args()
    if args.input_dir:
        CONFIG["input_dir"] = str(Path(args.input_dir).expanduser().resolve())
    if args.output_dir:
        CONFIG["output_dir"] = str(Path(args.output_dir).expanduser().resolve())
    if args.stage != "all":
        run_stage(args.stage)
    else:
        out = Path(CONFIG["output_dir"])
        out.mkdir(parents=True, exist_ok=True)
        statuses = []
        for stage in ["official_toy", "z_test", "fixed3d", "official_basic"]:
            cmd = [sys.executable, str(Path(__file__).resolve()), "--stage", stage,
                   "--input-dir", CONFIG["input_dir"], "--output-dir", CONFIG["output_dir"]]
            env = dict(os.environ, PYTHONHASHSEED=str(CONFIG["seed"]), MPLBACKEND="Agg",
                OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
            result = subprocess.run(cmd, env=env)
            statuses.append({"stage": stage, "exit_code": result.returncode,
                "status": "PASS" if result.returncode == 0 else "FAILED"})
            dump(out / "stage_status.json", statuses)
            if result.returncode:
                sys.exit(result.returncode)
