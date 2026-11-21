# scATAC 硬红线 × 常见误读速查表

运行时对照用。完整决策过程（调用 → 读 → 决定 → 红线）以 [`analysis-playbook.md`](analysis-playbook.md) 为准；
本表把每条红线压成"红线 / 为什么 / 常见误读 / 正确做法"四列，便于在每一步快速自检。

| # | 硬红线 | 为什么 | 常见误读 | 正确做法 |
| --- | --- | --- | --- | --- |
| 1 | **原始可及性出发**：peak calling、差异可及性、scrublet 都用 fragment/counts，不用校正/缩放/积分值 | 统计检验与 peak caller 假设原始计数分布；校正破坏均值-方差关系 | 拿 integration 后的嵌入或 scaled 矩阵去 call peaks / diff | 始终从 counts/fragment 层出发；校正值只用于可视化与聚类 |
| 2 | **backed AnnData 退出前必须 close()** | snapatac2 对象背在 HDF5 上自动同步，不 close 会损坏文件 | 以为像 scanpy 内存对象一样不用管；进程结束文件坏了 | 脚本链每步落盘并 close；notebook 末尾 close 格；内存态对象显式 write |
| 3 | **分析单元是细胞（barcode），不是 spot/bin** | scATAC 一个 barcode = 一份染色质；无空间坐标 | 把 scATAC 当空间数据找 obsm['spatial']；把 tile 当细胞 | 细胞级 QC/聚类；tile/peak 是基因组特征维度，不是分析单元 |
| 4 | **doublet 会抬高可及性，必须 scrublet** | 两个细胞共一个 barcode → 片段数与可及性虚高，伪装成"高活性"新类型 | 只看 n_fragment 高就当好细胞；忽略 doublet | scrublet 算 doublet_probability，filter_doublets 去掉；报告 doublet 率 |
| 5 | **TSSe 是 ATAC 特有信噪比，阈值来自分布** | 低 TSSe = 片段散布无核小体/TSS 结构，即便片段数够也是噪声 | 只用 n_fragment 过滤；TSSe 阈值凭默认 | 读 tsse 分布再定 min_tsse；kept_fraction < 一半停下重审 |
| 6 | **peak/tile 策略是生物学决策** | tile 大小、peak 宽度、是否 per-cluster call 决定分辨率与细胞类型特异调控 | 用默认 bin/全局 peak 一套到底 | 粗粒度组织域用大 tile；细胞类型特异调控用 per-cluster peak；说清选择与理由 |
| 7 | **gene activity 是代理，不是实测 RNA** | 由启动子/基因体可及性推断表达，稀疏且偏倚 | 把 gene activity 当表达量做差异表达声称 | gene activity 仅用于注释/marker 可视化；差异表达需 scRNA；标注是推断 |
| 8 | **label transfer 需要 scRNA 参考** | 没有参考就没有标签来源 | 无参考却声称做了 label transfer / 细胞类型注释 | 参考来自 SingleCell agent（counts+obs/pseudobulk 共享出口）；无参考只做 gene-activity marker 推断并说明 |
| 9 | **motif 富集需 FASTA + motif 库；结果是序列预测** | motif 富集检验 marker 区域里 TF motif 是否过代表，是序列层面的预测 | 把 motif 命中当实测 TF 结合；把 motif 富集冒充 chromVAR 逐细胞偏差 | 说清 motif 库（JASPAR/cis_bp）与背景；motif 是预测需实验验证；chromVAR 逐细胞 TF 偏差是 R 路缺口、未实现，要它报缺口走项目环境 |
| 10 | **样本是实验单元，不是细胞** | 条件级声称的重复是生物学样本，细胞是伪重复 | 把成百上千细胞当重复做条件差异 | 每组 ≥2 生物学重复才做条件级声称；否则只描述样本内模式 |
| 11 | **整合不要过校正；校正值不进差异检验** | 批次校正可能抹掉真实生物学差异 | 拿 integration 后的值做差异可及性 | 校正只用于联合聚类/可视化；差异可及性回到原始 counts，按样本建模 |
| 12 | **不做比对/Cell Ranger；缺基因组资源报缺口** | 原始 reads 处理是 nf-core wrapper；插件环境只读 | 在插件里跑 Cell Ranger / pip 装 episcanpy / 装 chromVAR / 下载 845MB FASTA 进插件环境 | 从 fragment/BAM/peak 起步（macs3 已随环境）；缺包/缺参考 → env_probe 报事实 + 项目环境 env_request 或指引下载到项目目录 |

## 使用方式

- 每跑一步前，扫一遍与本步相关的行（QC→1/4/5；降维聚类→1/3；peaks→1/6；注释→7/8；motif→9；整合/跨样本→10/11；任何重依赖→12）。
- 交付前对照第 2 行（对象都 close/落盘了吗）与第 12 行（"没跑"有没有被说成"没问题"）。
- 阈值与策略永远"先读分布/报告，再决定"，不是先验默认。引用工具 JSON 里的数字，不要自己重算。