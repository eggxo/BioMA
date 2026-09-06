# BioMA 用户与方法手册（中文）

本手册适用于 BioMA 0.8.x，统一说明项目模式和 7 个科学模块的功能、计算
原理、输入、全部可配置字段、运行示例、主要输出及解释边界。配置解析代码和
锁定依赖版本决定 BioMA 实际执行的计算；文末论文用于说明方法背景，不代表
BioMA 的结果已经获得独立生物学验证。

[English manual](USER_GUIDE.en.md) | [安装说明](../INSTALL.md) |
[统一输入清单](../INPUTS.md) | [软件结构](../ARCHITECTURE.md)

## 1. 先理解结果边界

- GF、RONA 和 MaxEnt 都是基于当前空间关系向未来环境外推。高值表示相对更大
  的变化或失配，不是灭绝概率、适合度下降百分比或因果效应。
- loadM/loadD 是 BioMA 脚本定义的 derived-allele 比值代理指标，不等于直接
  测得的个体适合度，也不是经典选择系数加权的遗传负荷。
- MAR 描述面积与多样性的静态幂律关系；WFmoments 描述栖息地损失后中性核苷酸
  多样性随时间的动态。二者回答不同问题。
- 综合脆弱性图把 5 个连续指标在当前项目的群体间相对缩放为 1-3 个方块。
  它是排序和可视化工具，不是具有跨物种绝对意义的新风险指数。
- 不同模块的原始数值尺度不能直接比较。跨 GCM 汇总使用未加权算术均值或配置
  的中位数，不代表模型概率权重。

## 2. 安装、检查和最小示例

```bash
conda-lock install --mamba --name bioma conda-lock.yml
mamba activate bioma
bin/install-bioma-r-deps.sh
python -m pip install --no-deps .
bioma --version
bioma doctor
```

填写项目配置后，再用 `bioma doctor project.ini --strict` 对已启用模块执行严格
检查；不带配置的 `--strict` 也会把未配置的可选组件计为失败。

生成不依赖服务器私有路径的 demo，并检查 7 个模块的输入连接：

```bash
python tests/data/demo/make_demo.py --outdir /tmp/bioma-demo
bioma project /tmp/bioma-demo/project.ini --dry-run
```

`--dry-run` 只做输入、参数、场景和依赖关系验证，不完成科学计算。真正运行：

```bash
bioma project project.ini
```

所有相对路径均相对于其所在 INI 文件解析。布尔值推荐使用 `true` 或 `false`。
源 VCF、表格、栅格和 mask 只读；输出 manifest 记录内容级 SHA-256。

## 3. 项目模式：一次连接全部模块

复制 `workflow.project.example.ini`，在 `[inputs]` 中填写共享输入，在
`[modules]` 中启用模块：

```bash
cp workflow.project.example.ini project.ini
bioma doctor project.ini
bioma project project.ini --dry-run
bioma project project.ini
```

命令行参数：

| 参数 | 含义 |
| --- | --- |
| `config` | 项目 INI 路径。 |
| `--dry-run` | 只完成预检查和有效配置快照。 |
| `--module NAME` | 只运行指定的已启用模块，可重复给出；名称为 `gf`、`rona`、`mar`、`load`、`niche`、`wfmoment`、`vulnerability`。 |
| `--overwrite` | 输入或配置变化时，把旧模块和项目元数据移动到 `_bioma_backups/` 后重跑。 |
| `--no-resume` | 本次忽略 `[project] resume`，不复用已完成模块。 |

### 3.1 `[project]`

| 字段 | 默认值 | 含义 |
| --- | --- | --- |
| `name` | INI 文件名 | 项目名称，写入报告和 manifest。 |
| `output_dir` | 必填 | 项目结果根目录。 |
| `resume` | `true` | 仅在有效配置和输入内容签名完全一致时复用完成结果。 |
| `stop_on_error` | `true` | 一个模块失败后立即停止；`false` 允许无依赖模块继续，并把依赖失败模块标为 blocked。 |
| `report` | `true` | 生成 `report/index.html`。 |

### 3.2 `[modules]`

`gf`、`rona`、`mar`、`load`、`niche`、`wfmoment` 和 `vulnerability` 的值是
各模块 INI 路径。使用 `false`、`no`、`off`、`disabled`、`none`、`0` 或留空
可禁用。模块顺序固定为上述顺序，综合脆弱性依赖 GF、RONA、load 和 niche。

### 3.3 `[inputs]` 统一输入字段

| 字段 | 用途及要求 |
| --- | --- |
| `adaptive_vcf` | GF/RONA 共用的二等位适应性位点 VCF，必须含 `GT`。群体 ALT 频率由 GF 自动生成。 |
| `adaptive_frequency` | 旧版兼容入口；新项目不要填写。仅在不运行 GF 时提供已有频率表。 |
| `whole_genome_vcf` | MAR 使用的全基因组、无缺失或已填补 VCF；不要与适应性 VCF 或 load VCF 混用。 |
| `load_vcf_dir` | derived 推断和 SIFT 注释后的 4 个 load VCF 所在目录。 |
| `population_samples` | 每个个体一行的主表；至少含 `sample_id`、`population_id`，可含 `lon`、`lat`。用于派生 GF sample 表、MAR 坐标、load keep-list 和 MaxEnt occurrence。 |
| `population_environment` | 每群体一行；含群体 ID、经纬度、可选分组及 `bio1`-`bio19`。GF 与 RONA 共用。 |
| `rona_unld_dir` | 19 个 `LD_BIO1.prune.in` 至 `LD_BIO19.prune.in`。这是 RONA 的真实位点选择输入。 |
| `climate_current` | 当代 BIO TIFF 目录，GF 和 niche 共用。 |
| `climate_future` | 未来情景根目录，GF、RONA 和 niche 共用；各模块仍检查自己的文件命名契约。 |
| `species_mask` | 物种范围 shapefile；同名 `.shx/.dbf/.prj/.cpg` 等 sidecar 一并计入哈希。 |
| `population_dir` | 可选的 load 群体 keep-list 目录；留空时从 `population_samples` 生成。 |
| `load_predictors` | 可选当前群体预测变量 CSV；留空时从统一群体表生成兼容文件。 |
| `load_future_dir` | load 未来 BIO CSV 目录。当前版本不能直接从 TIFF 自动生成，仍需提供。 |
| `mar_lonlat` | 可选逐样本 `ID,LONGITUDE,LATITUDE` 表；留空时从统一个体/群体表生成。 |
| `occurrence_csv` | 可选 `species,lon,lat` 点表；留空时从统一坐标生成去重点。 |
| `maxent_jar` | 合法取得的 `maxent.jar`；不随 BioMA 分发。也可由安装的 `dismo` 提供。 |
| `wfmoment_current_raster` | 当前二值栖息地/适生性栅格，推荐等面积 CRS。 |
| `pi_file` | 至少含 `species,pi_obs`。 |
| `structure_file` | 至少含 `species,fst_global_est`；`migration=auto` 时必需。 |
| `species_parameters` | 至少含 `species,z_gdar`，也可保存其他物种参数。 |
| `area_file` | 含 `Scenario,Area_km2`，必须有 `Scenario=current`；与 `future_masks_json` 二选一。 |
| `future_masks_json` | 情景名到未来二值 mask 路径的 JSON 对象；与 `area_file` 二选一。 |

### 3.4 `[shared]`

| 字段 | 默认值 | 含义 |
| --- | --- | --- |
| `species` | 空 | 覆盖 WFmoments 等模块的物种行键。 |
| `models` | 空 | 逗号分隔 GCM，覆盖 GF/RONA；为空时使用模块配置或发现值。 |
| `ssps` | 空 | 逗号分隔 SSP，可写 `245` 或 `ssp245`。 |
| `periods` | 空 | 逗号分隔时期，如 `2061-2080,2081-2100`。 |
| `seed` | 空 | 非空时覆盖支持随机种子的模块。 |
| `rscript` | 空 | 非空时覆盖模块 Rscript；模块需要不同环境时留空。 |
| `compute_python` | 空 | 非空时覆盖 WFmoments 的 Python。 |

### 3.5 `[integration]`

| 字段 | 默认值 | 含义 |
| --- | --- | --- |
| `period` | 条件必填 | 综合脆弱性自动连接的时期；若 `[shared] periods` 只有一个值可自动确定。 |
| `ssp` | 条件必填 | 综合脆弱性自动连接的 SSP。 |
| `ensemble_method` | `mean` | niche 自动结果名使用的 `mean` 或 `median`。 |
| `rona_variables` | 空 | 覆盖综合模块字段；`auto` 发现全部 `*_RONA`。 |
| `rona_summary` | 空 | 覆盖为 `mean`、`max` 或 `weighted`。 |
| `rona_weights` | 空 | `field:weight` 逗号列表，仅用于加权 RONA。 |
| `group_order` | 空 | 显式分组顺序；空值按输入首次出现顺序。 |
| `exclude_groups` | 空 | 要排除的分组；软件默认不排除任何特定物种群体。 |

### 3.6 `[doctor]`

`python`、`rscript`、`gdalinfo`、`gdallocationinfo`、`ogrinfo`、`java` 和
`maxent_jar` 都是可选的运行时覆盖路径。留空时从 PATH、模块 INI 和已安装
`dismo` 中发现。`bioma doctor project.ini --json` 输出机器可读报告；
`--strict` 把可选依赖缺失也视为失败。

## 4. Gradient Forest（GF）

### 4.1 功能、原理和指标

GF 先计算每个群体在适应性位点上的 ALT allele frequency，再以各位点频率为
响应、环境变量为预测因子训练 `gradientForest`。随机森林中沿环境梯度发生的
重要分裂被累积为每个环境变量的非线性变换。BioMA 在该变换空间中计算欧氏
距离：

- `local offset`：同一格点当前与未来环境的 GF 变换距离，表示原地维持所需的
  多位点组成变化。
- `forward offset`：从当前格点出发，在给定迁移半径内寻找未来环境中 GF 距离
  最小的目的地。值高表示即使迁移也难找到匹配的未来环境。
- `reverse offset`：对未来格点寻找当前范围内 GF 距离最小的来源地。值高表示
  未来栖息地缺少相近的当前遗传来源。

这三个值是模型空间距离，不是 allele-frequency 百分比、迁移概率或适合度。
BioMA 对每个坐标先要求所有选定 GCM 有有效值，再计算未加权算术均值。

### 4.2 输入与分析字段

| 字段 | 默认值 | 含义 |
| --- | --- | --- |
| `[inputs] vcf` | 必填 | 二等位适应性位点 VCF/VCF.GZ，必须含 `GT`；ALT 被当作频率统计等位基因，不自动判断 derived。 |
| `[inputs] samples` | 与下一项二选一 | TSV，列为 `sample_id,population_id`。 |
| `[inputs] sample_groups_dir` | 与上一项二选一 | 每群体一个文件、每行一个样本 ID；文件名作为群体名。 |
| `[inputs] coordinates` | 必填 | 空白分隔表；接受 `ID`/`population_id`、`lon`/`longitude`、`lat`/`latitude`，可有 `pop`/`group` 和已知 BIO 值。 |
| `[inputs] present_climate` | 必填 | 当代 BIO1-BIO19 TIFF 目录。 |
| `[inputs] future_climate` | 必填 | 完整 period x SSP x GCM 场景目录。 |
| `[inputs] current_mask` | 必填 | 当前分布/搜索 mask。 |
| `[inputs] future_mask` | `current_mask` | forward/reverse 搜索范围；扩大它会允许未来新栖息地。 |
| `[analysis] output_dir` | 必填 | 模块输出目录。 |
| `[analysis] models` | 自动发现 | 选定 GCM，逗号分隔。 |
| `[analysis] ssps` | 自动发现 | 选定 SSP，逗号分隔。 |
| `[analysis] periods` | 自动发现 | 选定时期，逗号分隔。 |
| `[analysis] expected_sites` | 空 | 独立 QC 的期望位点数；为空时从 VCF 推断，不应写死为测试数据数量。 |

### 4.3 全部 `[parameters]`

| 字段 | 代码默认值 | 含义与约束 |
| --- | --- | --- |
| `predictors` | `bio1,...,bio19` | GF 环境预测列；必须在提取后的群体环境和栅格中存在且在群体间有变异。 |
| `min_population_samples` | `3` | 每群体最少样本数，低于即停止。 |
| `warn_population_samples` | `5` | 低样本量警告阈值，必须不小于最小样本数。 |
| `ntree` | `500` | 每个响应位点随机森林的树数，至少 1；增大通常更稳定但更慢。 |
| `nbin` | `1001` | 累积分裂重要度的环境梯度 bin 数，至少 2；控制变换曲线分辨率。 |
| `corr_threshold` | `0.5` | `gradientForest` 的响应相关性阈值，范围 0-1；这是包内部响应筛选参数，不是 MaxEnt 的 BIO 共线性阈值。 |
| `max_level` | 自动 | 树的最大层级；自动为 `log2(0.368*n_populations/2)`，显式值必须大于 0。 |
| `seed` | `1` | GF 随机种子。 |
| `forward_radii` | `100,250,500,1000,inf` | forward 搜索最大球面距离（km）；`inf`/`unlimited` 表示不限距离。 |
| `initial_knn_k` | `64` | GF 空间最近邻搜索的初始候选数；算法会扩大候选集直到得到精确答案。 |
| `batch_size` | `10000` | offset 查询批次上限；影响内存和速度，不改变预期结果。 |
| `verify_sample` | `10` | 用暴力搜索复核的查询格点数；`0` 关闭复核。 |
| `tie_tolerance` | `1e-12` | GF 距离并列判定容差；并列时先选地理距离更近，再选稳定索引。 |
| `map_radius` | `unlimited` | 综合地图使用哪个 forward 半径，必须对应 `forward_radii` 输出标签。 |
| `minimum_models` | `2`（模板为 `3`） | 完整场景和 ensemble 至少需要的不同 GCM 数，至少 2。 |
| `supplied_bio_tolerance` | `1e-6` | 坐标表已有 BIO 与栅格重新提取值的允许差异；超出写警告，不替换栅格值。 |
| `rscript` | 自动发现 | 包含 `gradientForest`、`FNN`、`data.table`、绘图包的 Rscript。 |
| `gdalinfo` | `gdalinfo` | 栅格几何检查程序。 |
| `gdallocationinfo` | `gdallocationinfo` | 群体坐标栅格提取程序。 |
| `ogrinfo` | `ogrinfo` | mask 验证程序。 |

### 4.4 运行和输出

```bash
bioma run workflow.example.ini --dry-run
bioma run workflow.example.ini
```

主要输出：`01_frequency/population_alt_frequency.tsv` 是 GF/RONA 共用频率；
`03_model/all_gfmod.data` 是模型；`predictor_importance.tsv` 和
`response_performance.tsv` 分别给出变量重要性和位点 OOB R²；
`04_offsets/*/all_offsets.tsv.gz` 保存三类 offset、匹配坐标、距离和方向；
`05_offset_plots/*/ensemble_mean_offsets.tsv.gz` 是综合脆弱性默认读取的 GCM
均值；`06_forward_distance_plots/` 展示迁移半径敏感性。优先检查低正 OOB R²
位点比例、NoData、`no_candidate` 及模型覆盖率。

## 5. RONA

### 5.1 功能、原理和指标

RONA（risk of non-adaptedness）逐 BIO、逐位点拟合当前群体频率与环境的一元
线性回归 `p = alpha + beta*E`。对群体当前频率 `p_obs` 和未来环境
`E_future`，BioMA 计算 `abs(alpha + beta*E_future - p_obs)`，再以各位点回归
R²为权重求平均。高 RONA 表示模型认为需要更大的平均 allele-frequency 改变；
低值表示所需改变较小。它不建模选择强度、基因流、漂变、连锁或真实变化速率。

每个 BIO 使用自己的 `LD_BION.prune.in` 位点列表并与频率列取交集。当前实现
不会从 VCF 自动运行 PLINK，也不能用群体环境表替代 LD 列表。`SE` 是有效位点
绝对差值的普通标准误 `sd/sqrt(n)`，不是 GCM 间不确定性，也不是回归预测区间。

### 5.2 全部配置字段

| 字段 | 默认值 | 含义与约束 |
| --- | --- | --- |
| `[inputs] alt_frequency` | 必填 | 群体 x 适应性位点频率表；第一列可为 `Pop`，项目模式自动使用 GF 的 ALT 频率。 |
| `[inputs] unld_dir` | 必填 | 19 个非空 `LD_BIO1.prune.in` 至 `LD_BIO19.prune.in`。 |
| `[inputs] environment` | 必填 | 含 `ID,pop,lon,lat,bio1,...,bio19`；`ID` 与频率表群体匹配。 |
| `[inputs] future_climate` | 必填 | 场景目录名 `<period>-ssp<code>-<model>`，每目录恰有一个 `bioN.cut.tif`。 |
| `[inputs] mask` | 空 | 可选绘图 mask；不影响群体 RONA 计算。 |
| `[analysis] output_dir` | 必填 | 输出目录。 |
| `[analysis] models` | 空 | GCM 筛选；空值保留全部发现模型。 |
| `[analysis] ssps` | 空 | SSP 筛选；可带或不带 `ssp`。 |
| `[analysis] periods` | 空 | 时期筛选。 |
| `[analysis] expected_populations` | 空 | 可选的预期重叠群体数；正整数。dry-run 和正式运行均要求环境表 `ID` 与频率表首列实际交集数完全等于该值。 |
| `[parameters] rscript` | `Rscript` | 含 `data.table`、`terra`、`plotrix`、`raster`、`mgcv` 和绘图依赖的解释器。 |
| `[parameters] interpolate` | `true` | 是否用 GAM 生成连续空间图；不改变群体表。 |
| `[parameters] grid_step` | `0.1` | 插值网格经纬度步长，必须大于 0；越小越细、耗时和文件越大，但不增加观测信息。 |
| `[parameters] grid_k` | `15` | `mgcv::gam` 二维平滑基维数，至少 3；过大可能过拟合，过小可能过平滑。 |

```bash
bioma rona workflow.rona.example.ini --dry-run
bioma rona workflow.rona.example.ini
```

`weighted/` 保存每个 GCM 的 R²加权 RONA，`SE/` 保存位点标准误，
`weighted_SE/` 保存展示字符串，`ensemble_mean/` 保存所有选定 GCM 都有效时的
逐坐标算术均值，`maps/` 和 `boxplots/` 保存图。解释时应同时报告每个 BIO 的
有效位点数、回归 R²分布、LD 规则及 GCM 数；不要把 RONA 直接解释为存活概率。
运行前会检查环境表完整包含 `bio1`-`bio19`、频率表与环境表至少有 3 个
重叠群体，并确认每个 LD 列表至少有一个位点出现在频率表中。

## 6. MAR

### 6.1 功能、原理和指标

MAR（mutations-area relationship）把空间抽样或模拟栅格灭绝后的遗传多样性
与面积拟合为幂律 `S = c*A^z`。BioMA 用锁定的 `mar 0.2.0` 生成空间基因组
对象和灭绝轨迹，并分别拟合：`M`（观察到 ALT 的位点数）、`E`（只在保留区域
出现的 endemic 位点数）、`thetaW`（Watterson theta）和 `thetaPi`（核苷酸
多样性估计）。归一化剩余多样性为
`D_remaining = A_remaining^z`。在 `0<A_remaining<1` 时，较大的正 `z` 表示同等
面积损失对应更大的预测多样性损失。

MAR 是由当前空间样本拟合的静态关系，主要用于面积减少后的短期/即时代理；
不包含未来世代的漂变动态。采样密度、空间覆盖、VCF ascertainment 和灭绝方案
会改变 `z`。

### 6.2 全部配置字段

| 字段 | 默认值 | 含义与约束 |
| --- | --- | --- |
| `[inputs] vcf` | 必填 | 二等位全基因组 VCF，样本顺序与坐标一致；推荐 Beagle 后无缺失数据。 |
| `[inputs] lonlat` | 必填 | 仅含唯一的 `ID,LON/LONGITUDE,LAT/LATITUDE`，逐样本且顺序与基因型一致。 |
| `[inputs] scenario_file` | 空 | 可选情景表；含 `scenario,A_remaining,geom_job_id` 时向正式曲线添加点。 |
| `[analysis] output_dir` | 必填 | MAR 工作和输出目录。 |
| `[analysis] name` | `bioma_mar` | 输出对象基础名称。 |
| `[analysis] geom_id` | `7` | 从 `scenario_file` 选择 `geom_job_id`，不影响 MAR 拟合。 |
| `[parameters] maxsnps` | `auto` | `auto/all/空` 统计全部 VCF 记录；正整数表示上限，超限时由 `mar` 随机抽取，受 `randseed` 控制。 |
| `[parameters] scheme` | `random` | `random`、`inwards`、`outwards`、`northsouth` 或 `southnorth`。前两类中心方案围绕样本最密集格点；南北方案模拟方向性丧失。 |
| `[parameters] nrep` | `10` | 空间抽样/灭绝重复数，至少 1；方向性过程也可重复。 |
| `[parameters] xfrac` | `0.01` | 每步范围/栅格抽样比例，范围 `(0,1]`；栅格少于 100 格时包至少每步处理 1 格。 |
| `[parameters] quorum` | `true` | MAR 空间抽样时尽量要求抽样框含样本；不改变 extinction 栅格删除逻辑。 |
| `[parameters] randseed` | `123` | 位点下采样和空间重复的随机种子。 |
| `[parameters] marsteps` | `data,gm,sfs,mar,ext,plot` | `mar` 包阶段。BioMA 正式曲线读取 `extdflist`，因此常规完整运行必须保留 `data,gm,ext`；仅在确认已有兼容 RDA 时才跳过前置阶段。 |
| `[parameters] rscript` | `Rscript` | 含 `mar 0.2.0`、`SeqArray`、`sars` 和绘图依赖的解释器。 |

```bash
bioma mar workflow.mar.example.ini --dry-run
bioma mar workflow.mar.example.ini
```

`MAR_official_fit_params.tsv` 给出各指标的 `c` 与 `z`；
`MAR_display_curves.tsv` 是归一化曲线；`MAR_scenario_points.tsv` 是可选情景点；
`MAR_habitat_loss_curves.png/.pdf` 是正式图。RDA/GDS、日志和 manifest 用于审计。
报告结果时必须写清 VCF 筛选、`scheme`、重复数和面积定义。

## 7. loadM/loadD

### 7.1 功能、原理和指标

输入 VCF 的 ALT 必须已被极化为 derived allele。对每个二倍体个体，BioMA
把基因型转为 0/1/2 derived dosage，并分别计算每个位点类别的平均 derived
allele frequency：`Ps`（synonymous）、`Pn`（nonsynonymous）、`Pd`
（deleterious）。当前实现定义：

```text
loadM = Pn / (Pn + Ps)
loadD = Pd / (Pn + Ps)
```

loadM 表示 nonsynonymous derived burden 相对 nonsynonymous + synonymous
背景的比值；loadD 表示 SIFT deleterious derived burden 使用同一分母的比值。
这两个名称是 BioMA/来源脚本的代理指标，不应写成真实适合度负荷。VCF 未正确
极化、不同类别 callable 位点不一致或 SIFT 分类偏差都会影响值。

随后模块使用当前 `bio1`-`bio19`，先仅在训练集拟合标准化/SVD 变换，再用
`ranger` 随机森林预测 relaxed loadM 和 relaxed loadD。每次 80/20 外层划分，
训练集内重复 v-fold CV 以 RMSE 选择 `mtry`、`min_n`、`trees`，最终按
`cv_rsq_mean` 优先、`test_rsq` 次优选择一个 split 做未来预测。未来值是统计
外推，不表示未来突变真正积累的进化过程。

### 7.2 全部配置字段

| 字段 | 代码默认值 | 含义与约束 |
| --- | --- | --- |
| `[inputs] vcf_dir` | 必填 | 必须含 `strict.synonymous.vcf`、`strict.nonsynonymous.vcf`、`strict.deleterious_plain.vcf`、`relaxed.deleterious_with_warning.vcf`；可为 gzip 内容。 |
| `[inputs] population_dir` | 必填 | 每群体一个 keep-list，每行一个样本 ID。 |
| `[inputs] predictors` | 必填 | 当前 CSV，至少含 `pop,bio1,...,bio19`；地图还使用 `lon,lat` 和可选 `cluster`。 |
| `[inputs] future_dir` | 必填 | 未来预测 CSV 目录，每个文件至少含 `bio1`-`bio19`。 |
| `[inputs] mask` | 空 | 可选地图 mask；不影响 load 计算或 RF 拟合。 |
| `[analysis] output_dir` | 必填 | 输出目录。 |
| `[parameters] n_splits` | `20`（模板 `1000`） | 独立 80/20 外层划分数量，至少 1；越多越稳定且计算近似线性增加。 |
| `[parameters] cv_v` | `5` | 每个训练集的 CV 折数，至少 2。 |
| `[parameters] cv_repeats` | `2`（模板 `10`） | v-fold 重复数，至少 1。 |
| `[parameters] grid_size` | `30`（模板 `50`） | 每个 split 的 Latin hypercube 超参数组合数，至少 1。 |
| `[parameters] workers` | `4`（模板 `16`） | 并行 worker 数，至少 1；应服从调度器分配。 |
| `[parameters] scheme` | `svd` | 当前版本只允许 `svd`。 |
| `[parameters] seed` | `1234` | 外层划分、CV、网格和未来重建种子。 |
| `[parameters] rscript` | `Rscript` | 含 `tidymodels`、`ranger`、`future` 及绘图包的解释器。 |
| `[parameters] calc_python` | 当前 BioMA Python | load 计算脚本使用的 Python。 |
| `[parameters] calc_script` | 内置脚本 | 高级覆盖；替换脚本会作为科学输入记录 SHA-256。 |
| `[parameters] rf_script` | 内置脚本 | 高级覆盖 RF 调参实现。只有经过审查、接口完全相同的脚本可用。 |
| `[parameters] predict_script` | 内置脚本 | 高级覆盖未来预测实现。 |
| `[parameters] future_pattern` | `future_env_ssp*_mean.csv` | 在 `future_dir` 中选择未来文件的 glob。 |
| `[parameters] expected_future_files` | `0` | 正整数时强制文件数量；`0` 接受任意非空集合。 |

```bash
bioma load workflow.load.example.ini --dry-run
bioma load workflow.load.example.ini
```

`loads/` 保存个体/群体 strict 和 relaxed load 及样本、位点交集 QC；
`pop_geo_niche_predictors_from_TSS.csv` 是不修改原表的合并副本；两个 `RF_*`
目录保存每个 split 的 CV、测试、重要性和模型；`Figure_S24_loadM_tuning.*` 与
`Figure_S24_loadD_tuning.*` 分开输出；`future_predictions/`、`maps/` 和
`Figure_S25_load_maps.*` 保存未来预测。手稿中的 CV/test R²必须来自正式参数
运行，并报告均值、变异、样本量和空间外推限制，不能引用 smoke test 数值。

## 8. MaxEnt 生态位

### 8.1 功能、原理和指标

MaxEnt 用 presence 点与背景环境估计在约束条件下的最大熵分布。BioMA 先按用户
指定范围计算 BIO 相关性并迭代删除高度相关变量，再组合候选变量子集、feature
class 和 beta multiplier，通过 presence/background 的 k-fold AUC 或 TSS 选择
模型。最终模型用全部 occurrence 与背景点重拟合并投影到每个 GCM。

输出是 cloglog 连续适生性，不是出现概率的直接观测。`maladaptation = current
suitability - future suitability`：正值表示预测适生性下降，负值表示预测改善，
0 表示不变。当前版本明确不二值化。

### 8.2 全部配置字段

| 字段 | 默认值 | 含义与约束 |
| --- | --- | --- |
| `[inputs] work_dir` | 必填存在 | 兼容/工作目录标识；当前内置脚本把正式中间结果写入 `output_dir`。 |
| `[inputs] occurrence_csv` | 必填 | 至少含 `lon,lat`，mask 内至少 5 个具有完整 BIO 的唯一点；`species/cluster` 可用于标注。 |
| `[inputs] current_env_dir` | 必填 | 当代 BIO TIFF，文件名必须唯一识别 BIO 编号。 |
| `[inputs] future_root` | 必填 | `<period>-<scenario>-<gcm>` 目录集合。 |
| `[inputs] mask_shp` | 必填 | 校准背景、裁剪和绘图范围。 |
| `[inputs] maxent_jar` | 空 | 用户提供 jar；为空时使用 `dismo/java/maxent.jar`，两者至少有一个。 |
| `[variables] correlation_threshold` | `0.8` | 删除绝对相关性高于阈值的变量，严格在 0-1 之间。 |
| `[variables] correlation_method` | `pearson` | `pearson`、`spearman` 或 `kendall`。 |
| `[variables] correlation_scope` | `occurrence` | 在 `occurrence`、`background` 或最多 10,000 个 `mask` 样点上计算相关。 |
| `[variables] candidate_subset_sizes` | `4,6,8` | 相关过滤后候选变量组合大小；正整数，可大于剩余变量数（自动截断）。 |
| `[tuning] max_models` | `500` | 变量集 x feature x beta 的最大候选数，至少 1；超出时按种子抽样。 |
| `[tuning] background_n` | `10000` | 随机背景点目标数，至少 10；去除 occurrence 且要求完整 BIO。 |
| `[tuning] cv_folds` | `5` | presence/background CV 折数，至少 2，不超过 occurrence 数。 |
| `[tuning] feature_classes` | `L,LQ,LQH,LQHP` | MaxEnt feature 组合；字母 `L/Q/H/P/T` 分别为 linear/quadratic/hinge/product/threshold。 |
| `[tuning] beta_multipliers` | `0.5,1,2,3,4` | 正则化 multiplier 候选；一般越大模型越平滑。建议只用正值。 |
| `[tuning] selection_metric` | `auc` | `auc` 或 `tss`；均为 presence/background 判别指标，不直接衡量校准度。 |
| `[tuning] seed` | `123` | 背景点、候选子集、网格截断和 CV 分配种子。 |
| `[projection] periods` | `auto` | 时期列表；`auto` 从目录名发现。 |
| `[projection] scenarios` | `auto` | 如 `ssp245,ssp585`；`auto` 发现。 |
| `[projection] gcms` | `auto` | GCM 列表；`auto` 发现。 |
| `[projection] ensemble_method` | `mean` | GCM 栅格逐像元 `mean` 或 `median`；一个 GCM 时直接使用该层。 |
| `[projection] binary_output` | `false` | 当前版本必须为 `false`；不生成阈值化适生地。 |
| `[analysis] output_dir` | 必填 | 输出目录。 |
| `[analysis] rscript` | `Rscript` | MaxEnt 调参和投影环境。 |
| `[analysis] plot_rscript` | `Rscript` | 绘图环境，可与计算环境不同。 |

```bash
bioma niche workflow.niche.example.ini --dry-run
bioma niche workflow.niche.example.ini
```

`tables/correlation_matrix.csv`、`variables_after_correlation.txt` 和
`selected_variables.txt` 记录变量筛选；`tuning_metrics.csv` 与
`selected_model.tsv` 记录模型选择；`models/` 保存 fold 和最终模型；`rasters/`
保存 current、逐 GCM、ensemble 和 maladaptation；`figures/` 保存调参、当代、
未来和不适应性图。应检查 occurrence 空间偏倚、背景定义、外推/clamping、AUC/TSS
方差和 GCM 一致性。适生性变化不等于遗传不适应性。

## 9. WFmoments 2-D deme

### 9.1 功能、原理和指标

模块把当前二值栖息地聚合为 `nx * ny` 个 deme，并用四邻接二维 stepping-stone
Wright-Fisher 模型表示相邻 deme 间迁移。锁定的 `wfmoments` 数值程序利用等位
基因频率联合分布的一、二阶矩组成的常微分方程系统，计算移除 deme 后的
species-wide 中性核苷酸多样性 `pi`。当前零损失状态被归一化为 100%，输出
immediate、配置的 `time3`、`time5` 和可选 equilibrium 曲线。

独立的 GDAR 参考线为 `A_remaining^z_gdar`。它是面积幂律，不替代 WF immediate
曲线。正常连通的 edge contraction 中，较早时间通常高于较晚时间；random
fragmentation 下 species-wide pi 包含 deme 间分化，局部非单调并不自动表示
错误。

### 9.2 全部配置字段

| 字段 | 默认值 | 含义与约束 |
| --- | --- | --- |
| `[inputs] current_raster` | 必填 | 当前二值 mask，正值为适生/占据；推荐等面积 CRS。 |
| `[inputs] pi_file` | 必填 | `species,pi_obs`，`pi_obs>0`。 |
| `[inputs] structure_file` | 空 | `species,fst_global_est`；自动迁移率校准时必需。 |
| `[inputs] param_file` | 空 | `species,z_gdar`；也可改用参数中的 `z_gdar`。 |
| `[inputs] area_file` | 二选一 | `Scenario,Area_km2`，含 `current`；未来面积会限制到当前面积比例 0-1。 |
| `[inputs] future_masks_json` | 二选一 | 情景名到二值 TIFF 的 JSON；网格形状必须与当前栅格一致。 |
| `[analysis] species` | `my_species` | 在物种参数表中精确匹配的键。 |
| `[analysis] output_dir` | 必填 | 输出目录。 |
| `[analysis] plot_title` | 空 | 图标题；空时使用默认。 |
| `[parameters] compute_python` | 自动发现 | 含 `numpy,pandas,rasterio,wfmoments` 的 Python。 |
| `[parameters] rscript` | `Rscript` | 含 `ggplot2` 的绘图解释器。 |
| `[parameters] nx` | `20` | 东西方向 deme 数，1 至栅格列数。 |
| `[parameters] ny` | `20` | 南北方向 deme 数，1 至栅格行数。 |
| `[parameters] threshold` | `0.25` | block 内有效像元中正值比例达到该值才视为占据，范围 0-1。 |
| `[parameters] min_valid` | `0` | 一个 block 至少需要的有效像元数，非负整数。 |
| `[parameters] loss_mode` | `edge` | `edge` 为连续边缘收缩；`random` 为碎片化敏感性分析。 |
| `[parameters] direction` | `east_to_west` | edge 删除方向；允许 `east_to_west`、`west_to_east`、`north_to_south`、`south_to_north` 及 `eastwest/e2w/westeast/w2e/northsouth/n2s/southnorth/s2n` 别名。random 模式忽略方向。 |
| `[parameters] migration` | `25` | 相邻 deme 的模型迁移率，正数或 `auto`；不是直接的实测个体比例。 |
| `[parameters] migration_grid` | `0.1,0.3,1,3,10,25,50` | `migration=auto` 时的正数候选；选择模型东西半区 FST 与 `fst_global_est` 最接近者。 |
| `[parameters] fst_metric` | `hudson` | 自动校准使用 `hudson` 或 `nei` FST。 |
| `[parameters] theta` | `auto` | 正数或 `auto`；自动时按当前模型 equilibrium pi 线性缩放到 `pi_obs`。 |
| `[parameters] z_gdar` | 空 | 显式正数 GDAR 指数；非空时优先于 `param_file`。 |
| `[parameters] theta_probe` | `1e-4` | 自动迁移率/theta 标定的正初始 theta。 |
| `[parameters] time3` | `3` | 第一条未来 WF 曲线演化时间，通常解释为 3 generations。 |
| `[parameters] time5` | `5` | 第二条未来 WF 曲线时间，通常解释为 5 generations。 |
| `[parameters] mu` | `3.75e-8` | 正突变率；用于 `Ne = pi_obs/(4*mu)` 和自动 midterm 记录。 |
| `[parameters] midterm_generations` | `auto` | `auto` 记录 `Ne/2`；当前图仍只计算 `time3/time5`，该值不替代它们。 |
| `[parameters] replicates` | `1` | 删除顺序重复数；random 正式分析应大于 1，edge 同方向重复相同。 |
| `[parameters] seed` | `12345` | random 删除顺序种子。 |
| `[parameters] include_equilibrium` | `false` | 是否计算每个剩余 deme 集合的新 equilibrium，成本较高。 |
| `[parameters] plot_equilibrium` | `false` | 是否绘制 equilibrium；需同时计算。 |
| `[parameters] plot_raw` | `false` | 是否叠加原始点。 |

```bash
bioma wfmoment workflow.wfmoment.example.ini --dry-run
bioma wfmoment workflow.wfmoment.example.ini
```

`deme_table.tsv` 和 `deme_order.tsv` 记录栅格聚合与删除顺序；
`curve_replicates.tsv` 是不可改写的逐重复原始计算；`curve_summary.tsv` 是算术
均值；`migration_scan.tsv` 仅在自动迁移率时产生；`scenario_points.tsv` 把面积
情景插值到曲线；`metadata.*` 记录标定误差和参数。正式图的连续尾部可包含明确
标记的 display-only 平滑桥接，但 `Figure_wfmoment_2D_deme_curve_data.tsv`、
原始 summary 和 replicate 不被修改。100% habitat loss 的零 deme 是边界状态，
不是一个有观测生物学信息的普通 deme 状态。

## 10. 综合脆弱性

### 10.1 功能、原理和指标

模块在群体坐标上组合 niche continuous maladaptation、GF offset、RONA 汇总、
relaxed loadM 和 relaxed loadD。每个分组先取群体中位数，然后在本次分析所有
分组间对每项指标做 min-max 相对映射：

`score = round(1 + 2*(value-min)/(max-min))`，并限制为 1-3。

所有有限值相同时记 2，全部缺失时为 NA。高值固定解释为较高脆弱性，因此在
替换字段时必须确认方向一致。方块不是显著性、概率、权重或独立证据数量。

### 10.2 全部配置字段

| 字段 | 默认值 | 含义与约束 |
| --- | --- | --- |
| `[inputs] gf_offsets` | 必填/项目可 `auto` | 含 `lon,lat` 和 `gf_field` 的 GF ensemble 表。 |
| `[inputs] rona_ensemble` | 必填/项目可 `auto` | 含 `ID` 或 `pop` 及选定 `*_RONA` 列。 |
| `[inputs] load_predictors` | 必填/项目可 `auto` | 含 `pop,lon,lat`、分组、loadM/loadD 字段。 |
| `[inputs] niche_raster` | 必填/项目可 `auto` | 连续 maladaptation raster。 |
| `[analysis] output_dir` | 必填 | 输出目录。 |
| `[analysis] scenario_label` | `2061-2080 SSP245` | 图标题，不参与计算。 |
| `[parameters] gf_field` | `local_offset_mean` | GF 数值列；可选择 forward/reverse 均值，但必须在表中存在并在方法中说明。 |
| `[parameters] rona_variables` | `auto` | `auto/all` 发现所有 `*_RONA`，或显式逗号列表。 |
| `[parameters] rona_summary` | `mean` | 将多个 BIO RONA 合成一列：`mean`、`max` 或 `weighted`；逐群体忽略缺失。 |
| `[parameters] rona_weights` | 空 | `BIO1_RONA:0.3,BIO12_RONA:0.7`；非负，缺失值时对可用正权重重新归一化。 |
| `[parameters] rona_bio3_field` | 空、已弃用 | 旧配置兼容字段；新分析使用 `rona_variables`，软件不默认绑定 BIO3。 |
| `[parameters] rona_bio15_field` | 空、已弃用 | 旧配置兼容字段；软件不默认绑定 BIO15。 |
| `[parameters] loadM_field` | `mean_loadM_relax` | loadM 列。 |
| `[parameters] loadD_field` | `mean_loadD_relax` | loadD 列。 |
| `[parameters] group_column` | `cluster` | load predictor 表的分组列。 |
| `[parameters] group_order` | 空 | 显式逗号顺序；空/`auto` 保留首次出现顺序，未列出的实际分组追加并警告。 |
| `[parameters] exclude_groups` | 空 | 排除列表；空/`none` 不排除任何群体。Admixed 等规则只能写在物种配置中。 |
| `[parameters] nearest_tolerance` | `0.25` | 群体与最近 GF 网格点允许的最大经纬度欧氏距离，必须大于 0；不是 km。 |
| `[parameters] rscript` | `Rscript` | 含 `data.table,ggplot2,raster` 的解释器。 |

```bash
bioma vulnerability workflow.vulnerability.example.ini --dry-run
bioma vulnerability workflow.vulnerability.example.ini
```

`vulnerability_population_values.tsv` 保存未缩放群体值和 GF 匹配距离；
`vulnerability_rona_selection.tsv` 固化自动发现字段和权重；
`vulnerability_group_summary.tsv` 保存分组中位数和 1-3 分；
`Figure_multidimensional_vulnerability.png/.pdf` 是展示图。解释图时必须同时提供
原始表，避免把同一批数据衍生的指标当作完全独立证据。

## 11. 运行记录、质量控制和复现

每个独立模块的 `run_manifest.json` 记录有效配置、输入内容指纹、脚本指纹、
版本、场景和输出。项目模式另外写入 `00_project/input_contract.tsv`、
`input_manifest.json`、有效 INI、状态和 HTML 报告。目录哈希按稳定相对路径递归
计算；shapefile sidecar 作为一个逻辑数据集。路径相同但内容变化会使签名变化。

发表前至少保留：实际 INI、manifest、软件版本/tag、Conda lock、所有 QC/模型
性能表、GCM 列表、随机种子、位点/样本过滤规则和原始未插值结果。空间插值图
只用于展示，不增加独立观测，也不应用其平滑像元替代群体值做统计检验。

## 12. 方法文献

1. Ellis N, Smith SJ, Pitcher CR. Gradient forests: calculating importance gradients on physical predictors. *Ecology* (2012). https://doi.org/10.1890/11-0252.1
2. Fitzpatrick MC, Keller SR. Ecological genomics meets community-level modelling of biodiversity. *Ecology Letters* (2015). https://doi.org/10.1111/ele.12376
3. Rellstab C et al. Signatures of local adaptation in candidate genes of oaks with respect to present and future climatic conditions. *Molecular Ecology* (2016). https://doi.org/10.1111/mec.13889
4. Rellstab C, Dauphin B, Exposito-Alonso M. Prospects and limitations of genomic offset in conservation management. *Evolutionary Applications* (2021). https://doi.org/10.1111/eva.13205
5. Pina-Martins F et al. New insights into adaptation and population structure of cork oak using genotyping by sequencing. *Global Change Biology* (2019). https://doi.org/10.1111/gcb.14497
6. Sang Y et al. Genomic insights into local adaptation and future climate-induced vulnerability of a keystone forest tree in East Asia. *Nature Communications* (2022). https://doi.org/10.1038/s41467-022-34206-8
7. Exposito-Alonso M et al. Genetic diversity loss in the Anthropocene. *Science* (2022). https://doi.org/10.1126/science.abn5642
8. Lin M et al. marApp: An R package and web portal to calculate mutations- and genetic diversity-area relationship for conservation. bioRxiv (2025). https://doi.org/10.1101/2025.09.09.675155
9. Breiman L. Random Forests. *Machine Learning* (2001). https://doi.org/10.1023/A:1010933404324
10. Vaser R et al. SIFT missense predictions for genomes. *Nature Protocols* (2016). https://doi.org/10.1038/nprot.2015.123
11. Phillips SJ, Anderson RP, Schapire RE. Maximum entropy modeling of species geographic distributions. *Ecological Modelling* (2006). https://doi.org/10.1016/j.ecolmodel.2005.03.026
12. Elith J et al. A statistical explanation of MaxEnt for ecologists. *Diversity and Distributions* (2011). https://doi.org/10.1111/j.1472-4642.2010.00725.x
13. Mualim KS et al. Large future genetic diversity losses are predicted from conservation indicators even with habitat protection. *PNAS* (2026). https://doi.org/10.1073/pnas.2514371123

这些文献解释方法来源和限制。引用 BioMA 分析时还应引用实际使用的外部包及其
版本；正式 `CITATION.cff` 将在作者和论文信息确定后补充。
