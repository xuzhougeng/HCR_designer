# 多倍体背景下的探针设计问题与后续任务

更新日期：2026-10-15。代码基线：`main` 的 `5bdc4b7`；服务器已部署同一代码版本并通过真实请求验证，尚未推送 GitHub。部署记录与回滚材料见服务器 `/home/xzg/web_app/HCR_designer_deployments/20261014-subgenome/`。

本轮完成了**第二阶段的完整近匹配检索与目标策略**，以及**第三阶段的候选生成边界修正**：双臂配对改为完整无缺口替换扫描，28 对历史漏检全部补回；目标判定可按物种、严格转录本或显式允许集合配置，非允许配对的排除阈值单独保存和展示；候选生成不再按三臂长度预留空间，也不再提前占用位置。

**仍未完成：** 参考与同源关系元数据（亚基因组、同源组注释）、围绕差异区域的候选生成、插入缺失全库复核、逐对实验对应表，以及三探针迁移。当前结果属于计算候选，不能直接视为单拷贝专一或可下单清单。

## 已完成的基础工作

| 工作 | 当前状态 |
| --- | --- |
| 数据与历史记录 | 已对应 13 个目标、397 对历史候选、51 对已购买探针和完整合成序列，保留历史编号 |
| 参考库 | Bo/Bju cDNA 共 155,577 条，Bo 55,125 条、Bju 100,452 条；文件 SHA256 已固定 |
| 双臂判定 | 核对同一参考转录本上的 L/R 全长命中、方向和间隔，保留坐标及错配位置 |
| 筛选修复 | 使用参考库时，缺少评估、无目标完全配对、跨物种或未知物种非目标完全配对均排除；同物种其他转录本仍只参与排序 |
| 结果追溯 | 保存实际评分、筛选原因、参考信息和参数；JSON 坐标统一为 1-based inclusive |
| 网页与部署 | 序列与参考 → 参数与检查 → 结果复核；支持单条、批量、证据展开及 ZIP 下载，已在服务器验证 |
| 完整近匹配检索 | `src/common/paired_scan.cpp` 按鸽笼种子穷尽扫描，每臂最多 3 个替换；397 对历史候选的 2,157 个配对位点与独立基线逐位点一致 |
| 目标策略 | 物种模式、严格转录本模式、显式允许集合、亚基因组模式；非允许配对排除阈值可配置并写入结果 |
| 亚基因组归属 | 可按亚基因组分多个 FASTA 上传（每个文件一个亚基因组，服务器合并并生成归属文件），或为单一参考附带两列归属文件；归属只来自文件，不由编号推断 |
| 候选生成边界 | 枚举全部可放下两臂的起点，隔离测试中 40 nt、两臂各 15 nt、间隔 0 枚举出 11 个几何位置 |

已提交的 `1316cc0` 保存服务器原有桥接库更新，`f716f4a` 包含上一轮功能修复。三探针尚未迁移到新的双臂特异性判定。

## 当前问题和影响

以下 P0 表示下一阶段算法开发的首要依赖，P1 表示随后需要补齐的能力，P2 表示工程维护。优先级不代表已经完成。

| 优先级 | 问题 | 已有证据或实现 | 影响 |
| --- | --- | --- | --- |
| 已完成 | 检测目标层级可配置 | `TargetPolicy` 提供 `species` / `transcript` / `allowlist` / `subgenome` 四种模式及独立排除阈值；单条表单与批量 CSV 第 5 列可指定允许共同检测的转录本，亚基因组模式按归属文件确定允许集合，二者都与参考库核对 | 同源拷贝命中是预期信号还是脱靶，由每个任务的策略决定并写入结果 |
| 已完成 | 同物种内拷贝可直接排除 | 严格模式及允许集合下，非允许转录本的完全配对直接排除，原因码区分 `same_species_exact_pair` 与 `same_species_near_pair` | 目标为某个 `BjuVA…` 时，完全配对到另一个 `BjuVB…` 可按策略排除；但单拷贝专一的结论仍受参考库收录范围限制 |
| 已完成 | 带错配检索已穷尽 | `paired_scan.cpp` 把较长一臂切成 `K+1` 段种子，逐位点核验两条整臂；397 对历史候选在 0/1/2/3 替换下分别检出 155/326/372/388 对，与独立基线一致，28 对漏检全部补回 | 不再有已知的替换错配漏检；范围仍限于正向转录本、无插入缺失 |
| P0 | 参考与同源关系元数据不完整 | 亚基因组归属已可按文件接入（分文件上传或两列文件），但**归属文件本身尚未提供**；正式参考版本、品种、基因与转录本对应、同源组注释仍未确认 | 亚基因组模式随时可用，前提是给出归属文件；同源组（而非亚基因组）层面的允许集合仍需逐个列出编号；检索范围限于参考库已收录的转录本 |
| P1 | 候选生成尚未针对差异区域 | 边界与占位问题已修复（枚举全部起点，不提前占用），但仍是沿序列均匀枚举 | 候选覆盖完整的可设计位置，但未优先搜索物种或拷贝间的差异区域 |
| P1 | 插入缺失检查未进入线上流程 | 独立补查仅覆盖 9 对候选与 619 条已知跨物种序列的 665 次比较；扫描结果中 `indels_checked` 固定为 `false` | “未检出替换错配”不能覆盖插入缺失；现有补查也不是全参考库完整检索 |
| P1 | 排序与实验表现缺少对应验证 | 当前错配权重用于排序，未建立经过实验校准的交叉反应模型；实验标签主要是基因级 | 不能把分数当作特异性概率，也不能把基因级标签分配给每一对探针 |
| P2 | 原有测试及运行维护待整理 | 全套 51 项测试中有 5 项既有失败、4 项既有错误，与原始基线相同；同步 HTTP 计算、tmux 服务，无后台队列 | 后续需要固定测试语义；完整扫描的耗时与内存已测量（见下），但仍无后台队列与任务恢复 |

### 检索差距已消除

同一批 397 对历史候选、同一参考库；下表统计“至少存在一个跨物种配对的候选对数”，不是命中位点数，也不是实验失败数。

| 每臂最多替换数 | 独立无缺口扫描 | 原 BLAST 流程 | 当前完整扫描 |
| --- | ---: | ---: | ---: |
| 0 | 155 | 155 | 155 |
| 1 | 326 | 326 | 326 |
| 2 | 372 | 351 | 372 |
| 3 | 388 | 360 | 388 |

逐位点核对同样一致：2,157 个配对位点的转录本编号、坐标、每臂错配数、错配位置和匹配到的参考片段与独立基线完全相同，没有遗漏也没有多出。检索范围为正向转录本、全长、无插入缺失、按设计间隔配对。该范围内的完整检出，不等于覆盖所有实验结合机制。

155 对跨物种完全配对中，包括 7 对已购买探针及 48 对原程序自动保留探针，现已从新筛选中排除。按默认物种模式重新筛选仍得到 64 对，不是新构建或已验证的探针清单。原先因漏检而进入保留集合的 5 对，现已带正确证据重新判定：它们的跨物种配对都出现在每臂 2–3 个替换处，默认阈值（只排除完全配对）下仍然保留，但特异性排序分不再为 0。

| 历史编号 | 跨物种配对转录本数（0/1/2/3 替换） | 默认阈值下的判定 | 排序分 |
| --- | --- | --- | ---: |
| `113CLV1.C06-10` | 0 / 0 / 0 / 1 | 保留 | 2.0 |
| `113CLV1.C06-14` | 0 / 0 / 2 / 2 | 保留 | 12.0 |
| `113CLV1.C06-3` | 0 / 0 / 0 / 1 | 保留 | 2.0 |
| `49ABCG11.A-1` | 0 / 0 / 1 / 2 | 保留 | 10.5 |
| `79ABCG13.A-12` | 0 / 0 / 0 / 1 | 保留 | 3.0 |

把排除阈值提高到每臂 2 个替换会排除其中 2 对，提高到 3 个会排除全部 5 对；该阈值未经实验校准，默认仍为 0。

耗时与资源（本机 Apple Silicon、209 Mbase 参考）：397 对按基因分 13 次扫描合计 71.8 s，扫描进程峰值常驻内存 2.4 MB（流式读取参考，不整库载入）；397 对一次性扫描为 24.0 s。真实单条设计（2,320 nt 目标、314 对候选）扫描 17.1 s，含单臂 BLAST 与输出共 36.8 s。

数据来源：[本轮回归摘要](../outputs/phase2_regression/summary.json)、[原 BLAST 流程回归](../outputs/phase1_regression/summary.json)、[独立复核汇总](../../analysis/20261014_probe_audit/results/summary.json)。

## 相关文件

路径以本文件所在目录为基准。代码和文档属于 Git 仓库；参考、历史输入和 `outputs` 是本地或服务器数据，不随普通代码克隆提供。`../../analysis` 等链接指向仓库旁的本地复核资料。

### 算法与网页代码

| 文件 | 负责内容 | 下一阶段关注点 |
| --- | --- | --- |
| [paired_scan.cpp](../src/common/paired_scan.cpp)、[paired_scan.py](../src/common/paired_scan.py) | 完整无缺口替换扫描及其按需编译、超时、截断检测和运行统计 | 插入缺失检索；按需的并行或分块策略；服务器编译器与二进制管理 |
| [reference_annotation.py](../src/common/reference_annotation.py) | 亚基因组归属文件的解析、校验、缓存和覆盖度统计 | 扩展到 `gene_id`、同源组等其他注释列 |
| [paired_specificity.py](../src/common/paired_specificity.py) | 目标策略、配对分类、排除条件、特异性分数及单臂 BLAST 证据 | 接入亚基因组/同源组注释；未知注释处理；阈值的实验校准 |
| [designer.py](../src/dual_probe/designer.py)、[config.py](../src/dual_probe/config.py) | 候选生成和设计约束 | 差异区域优先；长度枚举策略；候选数量与耗时的平衡 |
| [probeset_select.py](../src/common/probeset_select.py) | 硬性排除、排序及间距选择 | 覆盖度评估；重叠候选增多后的选择策略 |
| [output.py](../src/dual_probe/output.py)、[create_split_probe.py](../scripts/create_split_probe.py) | 证据输出及双探针流程串联 | 三探针迁移；证据体积（单臂明细仍较大） |
| [blast_utils.py](../src/common/blast_utils.py)、[utils.py](../scripts/utils.py) | 原有单臂分析和索引构建 | 单臂 BLAST 现仅用于排序，可评估是否继续保留 |
| [split_workflow.py](../scripts/split_workflow.py)、[split_web.py](../scripts/split_web.py) | 输入校验、任务保存、目标策略、单条/批量 HTTP 入口 | 后台队列与进度；按同源组批量指定允许集合 |
| [split_wizard.html](../templates/split_wizard.html)、[split-wizard.js](../static/split-wizard.js) | 分步表单、目标策略输入与结果证据 | 展示亚基因组/同源组；证据表的近似配对排序与筛选 |
| [recheck_historical_pairs.py](../scripts/recheck_historical_pairs.py) | 不改变编号和序列的历史回归，含逐位点基线比较 | 新策略下的对照运行；插入缺失复核接入后的回归 |

### 数据与验证材料

| 文件 | 内容与用途 |
| --- | --- |
| [bol_bju_cdna.fasta](../uploads/bol_bju_cdna.fasta) | 当前线上与本地使用的 Bo/Bju cDNA 参考；SHA256 为 `5de12a5800a09889b73b901b12787c16fb6461d1c5a2de42e1864109726921d0` |
| [original_397_pairs.json](../data/probe_audit_20261014/original_397_pairs.json)、[purchased_51_pairs.json](../data/probe_audit_20261014/purchased_51_pairs.json) | 历史候选和已购买子集；保留原编号、靶向臂及完整合成序列 |
| [target_13_genes.fasta](../data/probe_audit_20261014/target_13_genes.fasta)、[provenance.json](../data/probe_audit_20261014/provenance.json) | 13 个历史目标及来源、校验值 |
| [probe_table.txt](../resources/probe_table.txt)、[Bridge_Probes1-2578.csv](../Bridge_Probes1-2578.csv) | 当前桥接探针运行表与来源数据，共 2,478 条；不能只凭编号替换历史桥接序列 |
| [historical_bridge_sequences.json](../data/probe_audit_20261014/historical_bridge_sequences.json) | 历史设计所用桥接序列；用于复现及完整合成序列核对 |
| [原合成表](../../物种特异性测试二-260914.xlsx)、[复核结果工作簿](../../outputs/20261014_probe_audit/移动mRNA探针复核_20261014.xlsx) | 实验原始清单和已生成的复核结果；无需重新提供已有序列 |
| [audit.json](../../analysis/20261014_probe_audit/results/audit.json)、[paired_hits_max3.tsv](../../analysis/20261014_probe_audit/results/paired_hits_max3.tsv) | 独立扫描的机器可读结果及逐位点明细，作为回归对照 |
| [scan_pairs.cpp](../../analysis/20261014_probe_audit/scripts/scan_pairs.cpp)、[verify_scan.py](../../analysis/20261014_probe_audit/scripts/verify_scan.py) | 独立无缺口扫描原型及暴力枚举对照测试；已作为线上扫描器的设计来源和回归对照，原型本身保持不变 |
| [prepare_indels.py](../../analysis/20261014_probe_audit/scripts/prepare_indels.py)、[check_indels.cpp](../../analysis/20261014_probe_audit/scripts/check_indels.cpp)、[indel_supplement.tsv](../../analysis/20261014_probe_audit/results/indel_supplement.tsv) | 有限范围的插入缺失补查原型和结果，不能直接宣称为全库复核 |
| [phase1_regression](../outputs/phase1_regression) | 原 BLAST 算法对历史候选的结果、逐臂证据及漏检清单，保留不动 |
| [phase2_regression](../outputs/phase2_regression) | 本轮完整扫描对同一批历史候选的结果、逐位点基线比较及耗时内存 |
| [phase2_e2e.log](../outputs/phase2_e2e.log)、[phase2_http_e2e.log](../outputs/phase2_http_e2e.log) | 真实参考上的单条设计耗时，以及物种/严格/允许集合三种策略的真实 HTTP 运行与 ZIP 一致性检查 |
| [phase2_subgenome_e2e.log](../outputs/phase2_subgenome_e2e.log) | 亚基因组模式的真实 HTTP 运行：两列归属文件、按文件上传合并、关系分类与排除数对比；所用示例归属文件按编号拼写生成，仅用于检验机制 |
| [test_paired_scan.py](../tests/test_paired_scan.py)、[test_paired_specificity.py](../tests/test_paired_specificity.py)、[test_reference_annotation.py](../tests/test_reference_annotation.py)、[test_split_workflow.py](../tests/test_split_workflow.py)、[test_probe_table.py](../tests/test_probe_table.py) | 扫描器暴力对照、配对判定与目标策略、归属文件解析、网页流程及桥接库测试，当前对应 3、17、5、18、1 项测试 |
| [phase1_test_comparison.json](../outputs/phase1_test_comparison.json)、[staging-verification.json](../outputs/deployment_20261014/staging-verification.json) | 原始代码的 9 项测试问题对照，以及上一轮部署版本的 42 项测试结果 |

### 运行与交接

- 项目概览：[README](../README.md)；历史问题及实验资料：[工作目录 README](../../README.md)。
- 线上入口：[双探针页面](http://10.10.87.96:6789/split)；运行命令与回滚：[deployment.md](deployment.md)。
- 部署记录：[deployment.json](../outputs/deployment_20261014/deployment.json)；真实 HTTP 验证：[live-verification.json](../outputs/deployment_20261014/live-verification.json)。
- 服务器项目目录：`/home/xzg/web_app/HCR_designer`；备份目录：`/home/xzg/web_app/HCR_designer_deployments/20261014-guided-ui/`。

## 后续任务与验收条件

### 第一阶段 固定目标规则和回归基线

- [x] **明确检测目标。** 已实现按任务配置：物种模式、严格转录本模式、亚基因组模式，或显式列出允许共同检测的转录本（单条表单、批量 CSV 第 5 列）。允许集合必须逐个列明或来自归属文件，不会默认放行同物种所有基因。**仍需实验方确认每个目标实际要区分到哪一层**，代码只提供机制。
- [ ] **整理参考元数据。** 亚基因组一列已可接入（`transcript_id`+`subgenome` 两列文件，或按亚基因组分文件上传）；仍需 `gene_id、species、homolog_group、reference_version` 的对应，记录来源和未知项，不以名称后缀代替注释证据。
- [x] **固定不变的对照集。** 397 对历史输入、51 对已购买子集、参考校验值、独立扫描结果和原 64 对重筛结果均保留；本轮结果写入 `outputs/phase2_regression`，未覆盖 `outputs/phase1_regression`，编号不变。回归脚本现在逐位点比较独立基线（`--audit-sites`），任一差异即断言失败。

验收：机制部分已满足，每个任务的目标、允许命中和必须排除集合可从结果中的 `target_policy` 复现；亚基因组/同源组注释仍未接入，属未确认项。

### 第二阶段 补足近匹配检索并接入目标策略 — 已完成

- [x] **完整替换检索已工程化。** `src/common/paired_scan.cpp` 为项目内新实现，不修改独立原型：种子长度按上限 10 nt 截断（长臂不再报错，50 nt 臂长可用），校验查询列、非 ACGT 碱基、FASTA 结构，输出末尾写 `#DONE` 统计行。`paired_scan.py` 按源码哈希自动编译到 `.tools/paired_scan/`，支持 `HCR_SCAN_BIN`、超时、非零退出和截断检测，并把结果转成应用统一的配对格式。
- [x] **可配置排除规则已建立。** `TargetPolicy` 的四种模式（物种、严格转录本、允许集合、亚基因组）各自定义允许集合；非允许配对的排除阈值 `exclusion_max_mismatches` 单独定义、随结果保存、在网页第二步可设置并在运行信息中展示。排序权重与该阈值分离，原因码区分完全配对与近似配对。
- [x] **检查范围与失败状态已保存。** 每个任务的 `search` 记录方法、参考、方向、每臂阈值、`indels_checked`、含模糊碱基窗口的处理、扫描序列数与碱基数、耗时和峰值内存。扫描失败、超时、截断、输出格式异常都抛出异常；网页仅在全部任务完成扫描时标记 `search_exhaustive`，否则按未评估提示。

验收结果：397 对历史候选在每臂 0/1/2/3 替换下的跨物种候选集合为 155/326/372/388，与独立基线一致；2,157 个配对位点的坐标、错配数、错配位置与匹配片段逐条相同。28 对漏检全部补回；原先进入保留集合的 5 对带正确证据重新判定，默认阈值下仍保留（见上表），未预设必须排除。22 条随机查询 × 多种边界/重复/模糊碱基构造的合成参考与暴力枚举完全一致（`tests/test_paired_scan.py`，含 0 和 3 两个阈值）。真实参考耗时与峰值内存见上。

### 第三阶段 改进候选生成与插入缺失复核

- [x] **两臂边界和搜索范围已修正。** 上界改为 `len(序列) - (2 × 最小臂长 + 臂间间隔)`，枚举每个起点且不再占用位置，右臂起点允许到序列末端。隔离测试中 40 nt、两臂各 15 nt、间隔 0 枚举出 11 个几何位置，末尾候选的右臂结束于第 40 位。真实目标（2,320 nt）候选数由 33 对增加到 314 对，全部仍需通过原有物理与复杂度条件。
- [ ] **围绕差异区域生成候选。** 将目标与明确的非目标同源序列比较，枚举覆盖差异位置的 L/R 长度和窗口，再通过完整特异性检查及质量过滤。保留无法找到足够候选的原因，不为凑数放宽规则。
- [ ] **定义插入缺失复核范围。** 明确允许的编辑数、每臂与连接处的插入缺失处理、方向和配对几何约束；扩展现有有限补查，并保留“未完成全库检查”的状态直到范围确实覆盖。

验收：边界部分已通过（11 个几何位置的隔离测试在 `tests/test_paired_specificity.py`）。差异区域优先与插入缺失复核尚未实施，`indels_checked` 仍为 `false`。新设计继续使用本次任务的编号空间，不复用历史编号。

### 第四阶段 实验对照与工程收尾

- [ ] **形成逐对或逐组合的实验对应表。** 将样本、实际混池、浓度、条件、原始图像和对照关联到完整探针对编号，再分析误检与无信号；2 个特异、6 个非特异、5 个无信号的现有基因级标签不能替代逐对验证。
- [ ] **整理 9 项旧测试问题。** 处理缺失的 `filter_probe` / `tcr` 模块测试、`check_poly_n` 接口、GC 单位、连接位点和小写序列语义；先确定预期行为，再修实现或更新过时测试。
- [ ] **评估运行方式。** 耗时已测量（单条设计约 37 s，其中扫描 17 s；按基因重复扫描同一参考是主要成本，批量多任务会线性累加），仍需决定后台任务、进度、失败恢复与清理策略；可考虑同一参考上合并多任务扫描。三探针另行迁移并验证，避免直接套用双臂规则。
- [ ] **部署本轮改动。** 服务器需可用的 C++ 编译器或预编译 `HCR_SCAN_BIN`，并在服务器 BLAST 2.14.1 与真实参考上重跑回归与真实请求后再切换。

验收：形成有证据和实验状态的候选报告；实验结果可对应到具体单对或组合。发布前完成新旧回归、真实单条/批量请求和导出一致性检查，并保留回滚版本。

## 需要补充或确认的内容

| 内容 | 最少需要的信息 | 用途 |
| --- | --- | --- |
| 每个目标的检测层级 | 要区分 Bo/Bju、某个亚基因组/拷贝，还是允许哪些转录本共同检测 | 定义允许命中与排除集合，是筛选策略的首要输入 |
| 亚基因组归属 | 两列文本（第一列转录本编号、第二列亚基因组），或按亚基因组分好的多个 FASTA 文件；两种方式都可直接在网页第一步提交 | 启用亚基因组模式，确定每个目标的允许集合 |
| 其余参考来源及同源注释 | 正式版本、品种、转录本与基因关系、同源组映射；现有 GFF/GTF 或对应表均可用于整理 | 固定参考范围，避免只靠编号推断关系 |
| 实验材料序列 | 若已有，提供实际品种的目标及相关同源拷贝序列或已知差异 | 评估参考与实验材料是否一致 |
| 实验组合与图像 | 实验编号、样本、完整探针对编号、混池组成、浓度、条件、图像及对照对应 | 将计算判定与实验表现对齐 |

已有探针序列、桥接表和 Bo/Bju cDNA 库可以继续使用，不需要重新提供。完整近匹配检索和目标策略机制已就绪（含亚基因组模式），**当前最关键的输入是每个目标的检测层级和亚基因组归属文件**：前者决定选择哪种模式，后者让允许集合不必逐个列出编号。同源组层面的区分仍需额外注释。

## 当前基线的复现入口

从 `HCR_designer` 项目根目录执行：

```bash
source .venv/bin/activate
python -m unittest discover -s tests -p test_paired_scan.py -v
python -m unittest discover -s tests -p test_paired_specificity.py -v
python -m unittest discover -s tests -p test_split_workflow.py -v
python -m unittest discover -s tests -p test_probe_table.py -v

# 使用新的输出目录，保留已有回归结果。逐位点比较独立基线。
python -m scripts.recheck_historical_pairs \
  --audit-file ../analysis/20261014_probe_audit/results/audit.json \
  --audit-sites ../analysis/20261014_probe_audit/results/paired_hits_max3.tsv \
  --output outputs/phase2_baseline_review

# 严格模式下的对照运行（不覆盖默认物种模式的结果）。
python -m scripts.recheck_historical_pairs \
  --policy-mode transcript \
  --output outputs/phase2_strict_review
```

回归脚本现在对三件事设硬性断言：每个阈值上的跨物种候选集合与独立基线一致、逐位点证据一致、保留集合中没有跨物种完全配对。执行成功表示替换错配已在"正向转录本、每臂 ≤3 个替换"的范围内完整检出；**插入缺失仍未复核**，`summary.json` 的 `indel_checked` 保持 `false`。

亚基因组模式的对照运行：

```bash
python -m scripts.recheck_historical_pairs \
  --policy-mode subgenome --subgenome-map /path/to/subgenomes.tsv \
  --output outputs/phase2_subgenome_review
```

归属文件缺少某个目标时脚本直接报错，不会退回其他判定。

全套测试：`python -m unittest discover -s tests` 共 62 项，53 项通过，5 项失败、4 项错误与原始基线（`outputs/phase1_test_comparison.json`）完全相同，没有新增失败。
