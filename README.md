# RNA 双端基因计数后端

基于 FastAPI + pysam 的 RNA-seq paired-end 片段级基因计数服务。提交 GTF 注释与
按 QNAME 分组（queryname-sorted）的双端 BAM，得到每个 `gene_id` 的片段计数。

## 计数规则

**GTF（`app/gtf.py`）**
- 只读取 `exon` 特征；`gene_id` 属性必填，缺失整份拒绝。
- 同一 `gene_id` 的外显子必须位于同一染色体、同一正负链，否则拒绝。
- 坐标由 GTF 的 1 起始闭区间转换为 0 起始半开区间；同基因重叠外显子合并。
- 上限：每个上传文件 20 MiB；外显子 50,000 条。

**BAM 片段（`app/bam.py`）**
- 以 QNAME 连续成组；忽略 secondary（0x100）与 supplementary（0x800）比对。
- 成组后该 QNAME 又在其他组之后重现（非连续重现）时整份拒绝。
- 片段必须恰好有 1 条 R1 与 1 条 R2 主比对。缺端、任一端未比对、标记
  duplicate（0x400）、QC fail（0x200）、MAPQ < `min_mapq` 或 NH > 1 均过滤；
  NH 缺失视为 1。
- 上限：100,000 个片段（QNAME 组）。BAM 非法或截断整份拒绝。

**归属（`app/counter.py`）**
- CIGAR 中仅 M、=、X 产生参考覆盖；D、N 只推进参考坐标（跨越内含子不计覆盖）。
- 链模式：`forward` 取 R1 同链、`reverse` 取 R1 反链、`unstranded` 忽略链。
- 双端命中的 `gene_id` 取并集：恰好 1 个计 1；多于 1 个记歧义（ambiguous）；
  空集记未分配（unassigned）。配对端互相重叠、命中同一外显子多次均不重复计数。
- 先过滤、后归属。返回全部基因（含零计数）及 `total_fragments`、`filtered`、
  `ambiguous`、`unassigned`、`assigned`，满足
  `filtered + ambiguous + unassigned + assigned = total_fragments`。

任何错误都不会发布部分计数；失败任务不落盘。

## 持久化与并发

- 结果写入 `$RNACOUNT_DATA_DIR`（默认 `./data`），每个任务一个 UUID 目录，
  含 `meta.json` 与 `counts.tsv`；通过临时文件 + `os.replace` 原子发布，
  服务重启后自动扫描恢复，历史任务仍可查询与下载。
- 任务 ID 使用 UUID4，内存注册表加锁；并发提交互不覆盖。

## 代码结构

```
app/
  config.py    限制与数据目录
  models.py    枚举与数据类、CountError
  gtf.py       GTF 解析、外显子合并、分链区间索引
  bam.py       BAM 分组、CIGAR 覆盖、读段级过滤
  counter.py   片段过滤与基因归属
  storage.py   原子持久化仓储与重启恢复
  main.py      FastAPI 路由
scripts/
  make_examples.py  生成示例 GTF/BAM
tests/         pytest 自测（含 HTTP 往返、持久化、并发）
examples/      示例数据
```

## 快速开始

```bash
# 1. 生成示例数据
.venv/bin/python -m scripts.make_examples examples

# 2. 启动
RNACOUNT_DATA_DIR=./data .venv/bin/python -m uvicorn app.main:app --host 127.0.0.1 --port 8000

# 3. 提交任务（multipart 表单）
curl -s -X POST http://127.0.0.1:8000/jobs \
  -F gtf=@examples/annotation.gtf \
  -F bam=@examples/reads.bam \
  -F min_mapq=10 \
  -F strand=forward

# 4. 查询与下载（用返回的 job_id）
curl -s http://127.0.0.1:8000/jobs/<job_id>
curl -s http://127.0.0.1:8000/jobs/<job_id>/counts.tsv
```

表单字段：
- `gtf`：GTF 文件（必填）
- `bam`：QNAME 分组的双端 BAM（必填，无需 `.bai` 索引）
- `min_mapq`：整数 0–255（必填）
- `strand`：`forward` | `reverse` | `unstranded`（必填）

## HTTP 接口

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| POST | `/jobs` | 提交计数任务，返回完整计数 JSON（201） |
| GET | `/jobs` | 列出已完成任务 |
| GET | `/jobs/{job_id}` | 查询任务与计数 |
| GET | `/jobs/{job_id}/counts.tsv` | 下载 TSV（`gene_id<TAB>count`，全基因含零） |
| GET | `/health` | 健康检查 |

错误码：输入非法/截断/超限返回 `400`（超 20 MiB 为 `413`）；参数越界
`min_mapq` 返回 `422`；任务不存在返回 `404`。

## 示例预期结果（min_mapq=10）

| 链模式 | GENE_A | GENE_B | GENE_C | filtered | ambiguous | unassigned | assigned |
| --- | --- | --- | --- | --- | --- | --- | --- |
| forward | 2 | 0 | 0 | 6 | 1 | 3 | 2 |
| reverse | 0 | 1 | 1 | 6 | 0 | 4 | 2 |
| unstranded | 2 | 1 | 1 | 6 | 1 | 2 | 3 |

总片段数恒为 12，四个分项之和等于总数。

## 测试

```bash
.venv/bin/python -m pytest -q          # 13 项自测
.venv/bin/python -m compileall app scripts tests
```
