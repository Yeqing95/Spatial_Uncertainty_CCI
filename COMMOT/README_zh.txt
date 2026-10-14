COMMOT 3D 最小复现包

本包用于 Aim 3A 的接口准备与初始验证，未运行大规模 benchmark，也未评价准确率或不确定性覆盖率。
测试环境为 Linux x86_64、Python 3.10.22；Mac 上的安装与运行尚需你确认。

1. 安装和运行

解压后，在终端 cd 到包含 run_commot.py 的 commot_3d_pilot 文件夹。
建议创建独立环境，避免改变刚修好的 Spateo / spatial-cci 环境：

    conda env create -f environment.yml
    conda activate commot3d
    python -m pip check
    python run_commot.py --stage all

如果你的 conda 平台暂时没有 Python 3.10.22 构建，可以先用 Python 3.10 的可用补丁版本；
这种情况下记录实际 Python 版本，不能把它称为完全相同的已测试环境：

    conda create -n commot3d python=3.10 pip -y
    conda activate commot3d
    python -m pip install -r requirements-lock.txt
    python -m pip check

可选 Jupyter kernel：
    python -m ipykernel install --user --name commot3d --display-name "Python (commot3d)"

在此 kernel 的 Notebook 中运行全部阶段：
    import subprocess, sys
    subprocess.run([sys.executable, "run_commot.py", "--stage", "all"], check=True)

全部阶段由脚本在独立进程中顺序执行，以便进程峰值内存不受前一阶段影响。
默认重新运行会更新 outputs 下的文件。保留本包实测输出时，请指定新目录：
    python run_commot.py --stage all --output-dir ./my_outputs

2. 参数和输入

所有科学参数集中在 run_commot.py 顶部 CONFIG 中。
example_inputs 内含此前 fixed 3D 数据的三个原始输入文件，未重新模拟数据：
    fixed_input/counts_raw.tsv.gz：cell × gene，含 cell_id 列。
    fixed_input/coordinates.tsv：含 cell_id、specimen_id、section_id、draw_id、x/y/z、units。
    ligand_receptor_table.tsv：含 lr_id、ligand、receptor、mechanism_id、radius_um 等。

使用你已有的完整模拟数据目录：
    python run_commot.py --stage fixed3d --input-dir /你的路径/aim3a_controlled_v1 --output-dir ./my_outputs

当前输入：400 cells、200 genes、5 sections、2 synthetic LR pairs。
这是此前以 scMultiSim 表达背景加自定义三维通讯机制生成的半模拟版本。
它可以验证接口、细胞对应关系与三维距离处理；不应标成 scMultiSim 原生 3D CCI benchmark。

推断输入为 raw counts 按每细胞总量归一化至 10000 后 log1p，原始 counts 保留在 layers['counts']。
所有 200 个基因均用于计算归一化分母，不进行 HVG 筛选，不输入真值边或真值通讯分数。
列名是 synthetic gene IDs，必须使用自定义 LR 表，不能直接用 mouse/human CellChatDB 替换。
gene_103 等名字中的下划线属于单个基因名，fixed3d 阶段必须 heteromeric=False。
official_toy / official_basic 中的 RecA_RecB 等真正复合物则使用 heteromeric=True。

默认 fixed_cutoff_mode='per_lr'，使用已声明的模拟参数：LR001 为 13 um，LR002 为 36 um。
这些半径是本次受控接口验证中已知的生成设定；并非从 LR 数据库推导，也不是估计出的生物距离。
若后续做盲评或未知半径实验，需重新规定参数获取方式。
COMMOT 0.0.3 支持 dict cutoff，但实现会索引完整 ligand × receptor 网格。
脚本补齐禁止 LR 组合的无效字典位置；它们仍由 df_ligrec 禁止，并未增加 LR pairs。
如果设为 common，则两个 LR 使用 common_cutoff_um；这将改变 LR001 的候选边范围。

3. 各阶段含义

official_toy：复现官方 spatial_communication API 文档的 4-cell 小例子，核对 c0→c2 两条
LR 边的分数为 0.6 和 0.9，验证复合物处理及矩阵方向。

z_test：4-cell 控制实验。保持表达、细胞 ID、LR、阈值及所有 x/y 坐标相同。
两组发送-接收细胞平面距离均为 0.5，距离阈值为 1.0。
仅将第一组接收细胞的 z 从 0 改为 5，第二组作为保持不变的阳性对照。
包含五种输入：xy、xyz_zero、xyz_raised、显式 3D distance、误用缓存的 2D distance。
输出显示目标边 5→0，而对照边维持 5；显式 3D 距离与自动 XYZ 计算一致。
缓存二维距离时目标边仍为 5，证明 obsp['spatial_distance'] 会覆盖 obsm['spatial'] 的几何。

fixed3d：对完整 400-cell fixed 输入运行一次 COMMOT，导出逐 LR 稀疏矩阵与非零边列表。
验证 ID 集合和排序、缺失 LR 基因、XYZ 有限值、坐标单位、稀疏矩阵形状、分数非负性、
发送/接收边际和，以及每个 LR 的三维距离限制。未进行真值准确率评价。

official_basic：复现官方 Basic usage 教程的数据加载、归一化、5-LR 推断与 Fgf1–Fgfr1
sender/receiver scatter 图（原教程前 6 个代码块）。没有运行后续向量场绘图与下游 DEG 分析。
使用原教程的 Visium V1_Mouse_Brain_Sagittal_Posterior 数据（3355 spots、32285 genes）。
保留教程的 native coordinates 和 dis_thr=500，没有额外宣称其物理单位已验证。
运行这一阶段会通过 Scanpy 下载约 30 MB 的官方示例输入；数据缓存在 download_cache。
它与我们的 fixed3d 测试是两个独立阶段，不将自定义模拟运行冒充官方教程复现。

4. 输出与解释

每个适用阶段：
    LR001.npz 等：scipy CSR 稀疏矩阵，行=发送细胞，列=接收细胞。
    cell_order.tsv：矩阵两轴使用的同一稳定 cell ID 顺序。
    lr_edges.tsv.gz：非零有向边，含 sender_id、receiver_id、lr_id、ligand、receptor、score。
    lr_summary.tsv：每个 LR 的非零边、自边和跨层边计数。
    runtime_memory.tsv：每次推断调用的 wall time 和调用前后进程高水位内存。
    run_report.json：通过的验证项、平台、包版本、完整阶段耗时及进程峰值 RSS。
fixed3d 还导出 fixed3d_commot.h5ad，包含 raw counts、归一化 X、固定 XYZ 和 COMMOT 输出。
tuple-key 距离阈值字典仅在保存 h5ad 时转换成 ligand/receptor/distance_um 表，不改数值结果。

score 是 COMMOT 的通讯/运输量分数，不是边存在的后验概率，也不是生成器中的真值信号 W。
稀疏文件中没有存储的条目表示数值为零；仅凭输出零不能区分结构不允许、表达缺失或分配结果。
本次保留 COMMOT 原生 self edges，并以 is_self 标记；模拟真值的 universe 不含 self edges。
以后评价时须明确处理此差异。事后移除自边不等价于在优化前禁止自边，因为资源竞争会改变。
分子与几何扰动、性能指标、多种 seed 和大 benchmark 均未开始。

5. 内存与复现范围

内存以 resource.getrusage 的 ru_maxrss 记录；Linux 使用 KiB，Mac 使用 bytes，脚本统一转 MiB。
stage 的峰值包含导入、读取、推断和保存；不是纯 solver 的新增分配量。
runtime_memory.tsv 中每次 inference_seconds 只计 spatial_communication 调用。
数值和耗时来自本包报告记录的 Linux 环境，不能当作你 Mac 的已测结果。
requirements-lock.txt 固定了本次环境中的全部 Python 包版本；未使用 monkey patch 或跳过依赖。
环境依赖检查通过。官方图已检查为非空且包含 sender/receiver 两个面板。

6. 来源

COMMOT PyPI 0.0.3：https://pypi.org/project/commot/0.0.3/
官方安装：https://commot.readthedocs.io/en/latest/installation.html
官方 Basic usage：https://commot.readthedocs.io/en/latest/notebooks/Basic_usage.html
官方 API toy：https://commot.readthedocs.io/en/latest/api/commot.tl.spatial_communication.html
COMMOT 论文：https://doi.org/10.1038/s41592-022-01728-4
源码快照 main commit：d117445bc07eaa19109c7609e97d9e35b26e99ca
本包同时提供安装版本实际执行的 _spatial_communication.py 副本与原 MIT license。
输入校验值见 fixed3d/input_manifest.json；其他来源和文件哈希见 provenance.json、SHA256SUMS.txt。
