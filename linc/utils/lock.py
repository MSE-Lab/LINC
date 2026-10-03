"""进程级文件锁, 防止并发对同一 `-o` 输出目录执行 LINC pipeline。

工程审查问题7(2026-10-03 与用户确认): `search()` 里的 `shutil.rmtree`+重建、
`FileRepository._save()` 的 `.tmp`+`os.replace` 均未加任何并发保护, 两个
`linc run -o 同一目录` 同时跑会互相覆盖彼此的中间文件/state, 产生不可预期
的损坏结果。

方案: 基于标准库 `fcntl.flock`(POSIX, macOS/Linux 均支持, 不引入第三方
依赖如 `filelock`), 以**非阻塞方式**获取独占锁; 若已被另一进程持有,
立即报错(不排队等待, 这是"互斥"场景而非"排队"场景)。
"""
import fcntl
import os
from contextlib import contextmanager
from typing import Iterator


class PipelineLockedError(RuntimeError):
    """`-o` 输出目录已被另一个正在运行的 LINC 进程占用。"""


@contextmanager
def acquire_pipeline_lock(output_dir: str) -> Iterator[None]:
    """在 `<output_dir>/.linc_work/pipeline.lock` 上获取非阻塞独占锁。

    整个 pipeline 执行过程(run/resume/append)应持有该锁直到命令结束,
    正常或异常退出都会在 `with` 块结束时自动释放。若锁已被占用,
    立即抛出 `PipelineLockedError`, 不阻塞等待。
    """
    work_dir = os.path.join(output_dir, ".linc_work")
    os.makedirs(work_dir, exist_ok=True)
    lock_path = os.path.join(work_dir, "pipeline.lock")
    fd = os.open(lock_path, os.O_CREAT | os.O_RDWR)
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise PipelineLockedError(
                f"输出目录 {output_dir} 当前已被另一个正在运行的 linc 进程"
                "占用(无法获取 .linc_work/pipeline.lock 的独占锁), 请等待"
                "该进程结束后重试, 不支持并发对同一 -o 目录执行 pipeline。"
            ) from exc
        yield
    finally:
        try:
            fcntl.flock(fd, fcntl.LOCK_UN)
        except OSError:
            pass
        os.close(fd)
