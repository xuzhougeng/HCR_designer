# Probe Generator

该脚本用于生成探针，并根据一些过滤条件进行探针筛选。最终输出满足条件的探针列表。

## 用法

1. 安装依赖

确保已安装以下依赖：
- Python 3.x
- 相关Python库（具体依赖见requirements.txt）
- BLAST软件（如果使用了BLAST功能）

2. 准备输入文件

- `seq.txt`：包含待处理序列的文本文件。

3. 运行脚本

```bash
python probe_generator.py --input seq.txt --probe-size 50 --initiator-type B1 --polyN 2 --min-gc 0.3 --max-gc 0.5 --output probes.txt --blastdb blastdb.fasta
```

参数说明：

--input：输入序列文件的路径。
--probe-size：探针的大小。
--initiator-type：选择的初始化器类型。
--polyN：Poly N过滤的阈值。
--min-gc：GC含量过滤的最小阈值。
--max-gc：GC含量过滤的最大阈值。
--output：输出文件的路径。
--blastdb：可选参数，BLAST数据库的路径。如果提供了该参数，将基于BLAST结果进行探针过滤。

查看结果

脚本会生成一个包含满足条件的探针列表的输出文件（probes.txt）。你可以查看该文件来获取生成的探针信息。

> 目前代码中BLAST的地址如果不在PATH中，需要手动设置HCR_prober_generator.py 中的 cmd = ""

## Flask应用

如果需要运行网页端，需要额外安装一个flask

```bash
pip install flask
```

运行方法

```
python app.py
```