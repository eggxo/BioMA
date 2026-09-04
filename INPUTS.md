# BioMA 输入文件清单与合并建议

本文档以当前服务器上的 7 个模块为准，区分用户原始输入、模块中间结果
和仅用于运行环境的配置。路径只是当前项目示例；正式项目应在
`workflow.project.ini` 中替换为自己的路径。

项目模式的 `[inputs]` 是统一输入层。用户填写的文件会保留在
`00_project/input_contract.tsv`；GF samples、MAR lonlat、load keep-list、
load predictors 和 MaxEnt occurrence 等兼容文件由 BioMA 写入
`00_project/generated_inputs/`，不修改用户原始文件。

项目同时写入 `00_project/input_manifest.json` 保存输入内容指纹：普通文件按
块流式计算 SHA-256，目录按稳定排序的递归文件清单计算聚合 SHA-256，并记录
每个子文件的 SHA-256。`.shp` 会自动纳入同名 `.dbf`、`.shx`、`.prj`、`.cpg`
等 sidecar。项目签名使用这些内容哈希，因此仅修改文件内容（即使路径不变）
也会触发重新确认；未使用 `--overwrite` 时，检查发生在生成兼容文件之前。

`input_contract.tsv` 的列为 `key`, `path`, `description`, `origin`, `kind`,
`size_bytes`, `sha256` 和 `entries_or_components`；目录和 shapefile 的完整
逐文件明细请查阅 `input_manifest.json`。

面向读者的输入只保留三条基因组数据流：GF/RONA 共用适应性位点 VCF
及其自动生成的群体 ALT 频率表；MAR 使用全基因组过滤 VCF；loadM/loadD
使用 derived + SIFT VCF。RONA 专用的 `alt.frq` 不再作为独立用户输入。

## 总览

| 模块 | 当前主要输入 | 是否可与其他模块共用 |
| --- | --- | --- |
| GF | 适应性位点 VCF、样本-群体表、群体环境表、当代/未来 BIO 栅格、mask | 环境栅格和 mask 与生态位/RONA 共用；适应性 VCF 不与 MAR 或 load VCF 合并 |
| RONA | GF 生成的适应性位点群体频率、群体环境表、未来 cut 栅格、mask | 与 GF 共用适应性位点频率、环境栅格和 mask；LD 文件属于当前实现细节 |
| MAR | Beagle 后全基因组 `mis0.9 + maf0.00001` VCF、逐样本经纬度表、可选 scenario 表 | 坐标表可由统一人口表派生；VCF 与 GF/RONA 适应性位点集不同 |
| loadM/loadD | derived 推断后、SIFT 注释的 4 个 VCF、群体 keep-list 目录、当前 BIO 预测表、未来预测 CSV、mask | keep-list、坐标和 BIO 变量可由统一人口/环境表派生；4 个 VCF 不能与其他 VCF 合并 |
| 生态位 MaxEnt | occurrence 点、当代/未来 BIO 栅格、mask、合法 `maxent.jar` | occurrence 可由群体坐标去重生成；栅格和 mask 可与 GF/RONA 共用 |
| WFmoments | 当前二值 habitat raster、物种 `pi_obs`、结构/FST、`z_gdar` 参数、面积情景表 | 物种参数可以合并为一个 species-parameter 表；面积情景可由未来 mask 预处理得到 |
| 综合脆弱性 | GF ensemble、RONA ensemble、load predictors、niche maladaptation raster | 全部应由项目运行自动连接，用户不应重复填写 |

## 逐模块输入

### 1. GF

当前模板：`workflow.example.ini`；完整 GF 工作流内部包含频率计算、气候
准备、GF 训练、三种 offset、模型均值和距离图。

- `vcf`：自适应位点 VCF，必须有 `GT`。请填写本地或集群上的实际路径。
- `samples` 或 `sample_groups_dir`：样本到群体的对应关系。推荐最终统一为
  两列表 `sample_id`, `population_id`；群体组别和坐标放入环境表。
- `coordinates`：`ID`, `pop`, `lon`, `lat`, `bio1`--`bio19` 的群体环境表。
- `present_climate`：当代 BIO1--BIO19 TIFF 目录。
- `future_climate`：未来场景目录，目录名为
  `<period>-<ssp>-<model>`，每个目录含 BIO1--BIO19。
- `current_mask`：当前物种分布 mask shapefile。
- `future_mask`：可选未来搜索 mask；留空时使用 `current_mask`。

### 2. RONA

当前模板：`workflow.rona.example.ini`。

- `alt_frequency`：由 GF 的适应性位点频率步骤自动生成，默认是 ALT allele
  frequency。项目模式不再要求用户单独提供 `alt.frq`；旧版
  `adaptive_frequency` 键仍仅作为兼容项接受，不应在新配置中填写。
- `unld_dir`：19 个 `LD_BIO1.prune.in`--`LD_BIO19.prune.in` 文件目录。
  这些文件通常由适应性位点 VCF 经 PLINK `--indep-pairwise` 生成，RONA
  计算会逐 BIO 读取并将位点与频率表列取交集。`population_environment.tsv`
  只提供 BIO 环境预测值，不能替代基因型 LD 筛选；每个列表的大小和
  SHA-256 会写入 RONA manifest。后续若自动从 VCF 生成 LD，窗口、步长和
  r² 必须作为显式参数记录，不能从环境表猜测。
- `environment`：与 GF 相同的群体环境表。服务器中 GF 与 RONA 两份文件
  SHA-256 都是
  `cc7520e70e86334bc70d844c0aca3c3586a6661044ab4f833d82d7f8a50d42c0`，
  下一步可以只保留一份。
- `future_climate`：未来 cut 栅格根目录；场景目录和 BIO 文件名需符合
  RONA 脚本的规则。
- `mask`：绘图 mask，可与 GF/生态位共用。

### 3. MAR

当前模板：`workflow.mar.example.ini`。

- `vcf`：MAR 输入 VCF，`maxsnps = auto` 会自动统计位点数。
- `lonlat`：逐样本 `ID`, `LONGITUDE`, `LATITUDE` 表。
- `scenario_file`：可选的面积/情景表，仅用于曲线上的情景点。

MAR 的逐样本坐标可由 GF 的样本设计和群体坐标自动生成。MAR VCF 与 GF
VCF 只有在位点集合、样本名和基因型编码完全一致时才可共用；当前 6410
位点 MAR VCF 和 6140 位点 GF VCF 应继续分开。

### 4. loadM/loadD

当前模板：`workflow.load.example.ini`。

- `vcf_dir`：必须恰好包含：
  `strict.synonymous.vcf`、`strict.nonsynonymous.vcf`、
  `strict.deleterious_plain.vcf`、`relaxed.deleterious_with_warning.vcf`。
  这些是已经完成 derived 判断和 SIFT4G 注释的输入；本模块不重新推断
  ancestral/derived 状态。
- `population_dir`：一个 keep-list 文件对应一个群体，文件名为群体名。
  可由统一样本设计自动生成。
- `predictors`：当前群体 BIO、坐标和其他 RF 预测变量表；load 结果会写入
  副本，不修改原表。
- `future_dir`：未来预测环境 CSV 所在目录；`future_pattern` 默认匹配
  `future_env_ssp*_mean.csv`。
- `mask`：load 地图 mask，可与其他模块共用。

4 个 derived + SIFT VCF 不能与 GF/RONA 的适应性 VCF 或 MAR 的全基因组
VCF 直接合并；它们对应不同的注释类别和 load 计算语义。未来 CSV 可以由
统一气候情景目录自动生成，
  但当前版本仍要求 CSV 格式的独立输入。

### 5. 生态位 MaxEnt

当前模板：`workflow.niche.example.ini`。

- `occurrence_csv`：`species`, `lon`, `lat` 点表。当前点来自同一批群体坐标，
  可由环境表按坐标去重生成。
- `current_env_dir`：当代 BIO1--BIO19 TIFF 目录。
- `future_root`：未来 `<period>-<ssp>-<gcm>` 目录集合。
- `mask_shp`：物种 mask；可与 GF/RONA 共用。当前 GF 和 MAXENT 的 mask
  `.shp` 文件 SHA-256 相同。
- `maxent_jar`：用户依法取得的 MaxEnt jar。它不随 BioMA 分发，运行时会
  使用该文件并记录校验值。

相关性阈值、候选变量子集、feature class、beta multiplier、CV 折数和
`mean/median` ensemble 都是参数，不是输入数据文件。

### 6. WFmoments 2-D deme

当前模板：`workflow.wfmoment.example.ini`。

- `current_raster`：当前二值 habitat/suitability raster。
- `pi_file`：物种级 `species`, `pi_obs` 表。
- `structure_file`：可选物种级 FST 表；`migration = auto` 时需要。
- `param_file`：可选物种级 `z_gdar`/MAR/FST 参数表。
- `area_file`：`Scenario`, `Area_km2` 表；或使用 `future_masks_json` 提供
  未来二值 mask。

建议将 `pi_obs`, `fst_global_est`, `z_gdar`, `z_mar` 整理为一个物种参数
表，而不是每次分别维护 3 个文件。当前分析使用 `area_file`，原始 raster
和表都只读，不会被绘图过程改写。

### 7. 综合脆弱性

当前模板：`workflow.vulnerability.project.example.ini`。项目模式下四项
都填 `auto`：

- GF：`ensemble_mean_offsets.tsv.gz`；
- RONA：项目内 `ensemble_mean` 表；
- load：项目内带 loadM/loadD 的 predictor 副本；
- niche：项目内连续 maladaptation ensemble raster。

用户只需在项目配置中指定 integration period、SSP 和 ensemble method，
不再重复填写这 4 个路径。测试配置可显式设置 `exclude_groups`（例如
排除混合群体），正式分析由用户按研究设计配置；软件默认不排除任何群体。

综合脆弱性参数也不绑定某个物种：`rona_variables` 可列出要汇总的
`*_RONA` 字段，或设置为 `auto` 自动发现；`rona_summary` 支持 `mean`、
`max` 和 `weighted`（加权时填写 `rona_weights`）。`group_order` 为空时按
输入表首次出现顺序排列，`exclude_groups` 为空时不排除任何群体。

## 当前最适合合并的输入

下一阶段可以先把下面 5 类文件收缩成统一输入层：

1. `population_samples.tsv`：个体样本到群体的主表；自动派生 GF sample
   表、MAR lonlat、load keep-list 和 MaxEnt occurrence 点。每个个体一行。
2. `population_environment.tsv`：唯一的群体 BIO1--BIO19 环境表；自动适配
   GF、RONA、load 和 MaxEnt 的表格要求。每个群体一行，不要重复成个体表。
3. `climate_catalog/`：唯一的当代/未来 BIO 栅格及 period/SSP/GCM 清单；
   自动建立 GF、RONA 和 MaxEnt 的目录视图，并可派生 load 的未来 CSV。
4. `species_parameters.tsv`：`pi_obs`, `fst_global_est`, `z_gdar` 等物种
   参数；自动生成 WFmoments 所需的兼容表。
5. `species_mask`：只保留一份 mask 及其 `.dbf/.shx/.prj` 配套文件。

不建议现在强行合并三类 VCF。GF/RONA 的群体 ALT 频率由适应性 VCF 的
GF 频率结果自动适配；MAR 全基因组 VCF 和 load 的 derived/SIFT VCF 仍应
保留独立入口，并在配置中明确标注位点处理历史。RONA 的 LD 选择不能用
环境表替代，需等候独立的 VCF/LD 规则设计。
