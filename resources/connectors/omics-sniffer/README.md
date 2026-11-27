Omics Sniffer 是跨源公开组学数据检索引擎，统一汇总 GEO、SRA、ArrayExpress、BioStudies、PRIDE、Massive/ProteomeXchange、CZ CELLxGENE、ENCODE、GDC/TCGA、GTEx、MetaboLights、Metabolomics Workbench、Zenodo、Figshare 等 30+ 公共科学数据仓库的标准化数据集、样本与文件索引。

### 能做什么

- 查看全局跨库目录规模与各数据源/组学模态覆盖度（`get_catalog_stats`）。
- 按关键词、疾病、组织、细胞类型、物种、检测技术（Assay）、组学模态、仓库来源、开放权限及最低样本量多维组合检索数据集（`search_omics_datasets`）。
- 根据统一数据集 ID 或原始入藏号（如 `GSE...`、`E-MTAB-...`、`PXD...`、`TCGA-...`）获取完整研究元数据、样本统计、文件统计与关联文献（`get_dataset`）。
- 浏览数据集下的样本元数据（组织、疾病、性别、年龄等临床/生物学属性，`list_samples`）及可下载数据产物列表（矩阵、原始测序或受控状态标记，`list_artifacts`）。

### 科研中的用法

- 在课题立项或验证阶段一次性跨库盘点特定癌种、器官或罕见病的公开转录组、单细胞、空间组学、表观组、蛋白组与代谢组队列。
- 快速筛选达到统计效能门槛（如 `min_samples >= 50`）且具备开放处理矩阵的数据集，直接衔接下游生物信息学分析流。
- 通过 DOI 或 PubMed ID 反查论文对应的多模态公开原始数据与补充矩阵入藏号。

### 使用说明

- 无需登录，也无需 API key，连接后即可检索公开元数据快照。
- 后端数据由 `data-pipelines/omics-sniffer` 增量同步管线定时拉取并原子发布为只读 Parquet 快照；受控访问数据（Controlled Access）仅返回元数据且 `url` 置空。
- 服务主页：[Phi Packages](https://github.com/dxsbiocc/phi-packages)。复用各仓库原始数据时，请遵守对应数据仓库及原始作者的数据使用许可与引用规范。
