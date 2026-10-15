# 服务器部署与回滚

服务器为 `xzg@10.10.87.96`，项目位于 `/home/xzg/web_app/HCR_designer`。现有服务使用 `tmux` 会话 `HCR_designer` 的 `0.0` 窗格、TCP 6789 端口；保持该运行方式。更新后用专用启动脚本关闭调试器及自动重载，避免文件更新过程中触发半完成版本启动。

```bash
cd /home/xzg/web_app/HCR_designer
HCR_PYTHON=/home/xzg/mambaforge/bin/python \
HCR_BLAST_BIN=/home/xzg/mambaforge/bin \
HCR_HOST=0.0.0.0 bash scripts/run_web.sh
```

服务器现有 Python 3.10.12、Biopython 1.81、NumPy 1.26.4、Flask 3.0.0、pandas 2.1.1、Requests 2.31.0 满足当前代码依赖。BLAST+ 为 2.14.1；本地开发的 BLAST+ 为 2.17.0，部署前应在服务器版本上完成真实参考库测试，不直接复制 macOS 可执行文件。无需修改共享 Python 环境。

2026-10-15 的部署（`5bdc4b7`）用 `git bundle` 把发布提交带到服务器，在 `HCR_designer_deployments/20261014-subgenome/staging` 的独立克隆里，用服务器的 Python 3.10、BLAST 2.14.1 和真实参考完成测试、历史回归和 HTTP 验证，再对运行目录做 `merge --ff-only`；扫描器在切换时预先编译，首个用户请求不承担编译耗时。该目录的 `before.bundle` 与 `before-code.tar.gz` 保存了部署前状态，回滚时 `git reset --keep f716f4a` 后按下面的命令重启即可，参考索引未变动，无需恢复。服务器上单次设计约 2 分钟，其中配对扫描约 32 s，仍是同步请求。

双臂配对扫描器需要 C++17 编译器：首次运行会用 `c++ -std=c++17 -O2` 把 `src/common/paired_scan.cpp` 编译到 `.tools/paired_scan/<源码哈希>/`，该目录不随代码提交，应用进程需要对它有写权限。服务器没有编译器时，在同架构机器上编译后放到服务器并设置 `HCR_SCAN_BIN` 指向该文件；macOS 上编译的二进制不能复制到 Linux 服务器。两者都不可用时应用会直接报错，不会退回近似检索。部署前先在服务器上执行一次历史回归（`python -m scripts.recheck_historical_pairs --audit-file … --audit-sites …`），确认 155/326/372/388 与逐位点结果一致。

本轮升级需要带原始序列编号的 BLAST 索引。`uploads/bol_bju_cdna.fasta` 保持原文件，SHA256 应为 `5de12a5800a09889b73b901b12787c16fb6461d1c5a2de42e1864109726921d0`。在独立部署目录用以下命令生成新索引，验证后再替换索引文件；旧索引与代码一起备份：

```bash
/home/xzg/mambaforge/bin/makeblastdb \
  -in /home/xzg/web_app/HCR_designer/uploads/bol_bju_cdna.fasta \
  -dbtype nucl -parse_seqids -out /path/to/staging/bol_bju_cdna.fasta
```

部署过程：

1. 保存服务器 HEAD、工作目录差异、已有未跟踪文件及旧索引，并与本地同步清单核对。不要覆盖部署期间产生的新改动。
2. 在独立目录检出待部署提交，使用服务器依赖和参考数据执行新增测试、完整旧测试对照及真实设计请求。只在验证通过后更新运行目录。
3. 停止 `HCR_designer:0.0` 窗格内的旧应用。用 Git stash 保存原有桥接库修改，再将 `main` 快进至待部署提交；不要对其他分支或服务使用强制重置。桥接库更新应已有独立提交，因而不再重复应用 stash。`test_debug.py` 保留为原有未跟踪调试文件。
4. 放入验证过的参考索引，以上述命令启动同一窗格。检查 `/split`、`/split/batch`、其他原有入口及真实设计、结果和下载接口。
5. 保存部署版本、备份路径、验证结果及重启命令。参考、历史实验数据、生成结果和本地环境不放入代码提交。

备份存放在项目同级的 `HCR_designer_deployments/<部署标识>/`，避免与网页上传和输出混在一起。回滚时先停止该应用，确认没有部署后的用户改动，再恢复备份代码和旧索引，或检出保存的旧版本并恢复备份中的原有未提交文件；然后按记录的旧命令启动。不要删除 `uploads/` 中的参考、用户上传或新产生的结果。

双探针设计在应用进程内用单个后台线程排队执行（`scripts/split_web.py` 的 `JobRunner`），一次只算一个任务；状态与进度写在 `uploads/split_results/<任务>/request.json` 和 `progress.json`。旧的一次性 `POST /split`、`/split/batch` 仍同步返回 ZIP。

**重启前先确认没有进行中的任务**：`curl -s -H 'Accept: application/json' http://127.0.0.1:6789/jobs` 中没有 `queued` 或 `running`。进程退出时正在计算或排队的任务不会继续；下次启动会把它们标记为 `interrupted`，用户可在任务页按原参数重新运行。仍是现有 tmux 服务，没有开机自启或其他进程管理；多进程 WSGI 部署前需要把队列换成独立进程或外部队列，否则每个进程各有一个队列。近匹配不穷尽、插入缺失未复核的算法限制见主 README。
