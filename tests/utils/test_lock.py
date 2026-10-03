"""linc.utils.lock 的单元测试(工程审查问题7, 2026-10-03 与用户确认):

进程级文件锁, 防止并发对同一 `-o` 输出目录执行 LINC pipeline。
"""
import fcntl
import os

import pytest

from linc.utils.lock import PipelineLockedError, acquire_pipeline_lock


def test_acquire_lock_creates_lock_file_under_linc_work(tmp_path):
    """获取锁应在 `<output_dir>/.linc_work/pipeline.lock` 创建锁文件。"""
    output_dir = str(tmp_path / "out")
    with acquire_pipeline_lock(output_dir):
        lock_path = os.path.join(output_dir, ".linc_work", "pipeline.lock")
        assert os.path.exists(lock_path)


def test_acquire_lock_twice_in_same_process_does_not_deadlock_is_reentrant_safe(tmp_path):
    """同一进程内先获取锁、正常释放后, 应能再次获取(验证 finally 释放逻辑
    正确, 不会残留锁导致后续调用永久失败)。
    """
    output_dir = str(tmp_path / "out")
    with acquire_pipeline_lock(output_dir):
        pass
    with acquire_pipeline_lock(output_dir):
        pass


def test_acquire_lock_raises_when_already_held_by_another_fd(tmp_path):
    """若锁已被(模拟的)另一进程持有(通过独立 fd 手动 flock), 非阻塞获取
    应立即抛出 `PipelineLockedError`, 而不是阻塞等待。
    """
    output_dir = str(tmp_path / "out")
    work_dir = os.path.join(output_dir, ".linc_work")
    os.makedirs(work_dir, exist_ok=True)
    lock_path = os.path.join(work_dir, "pipeline.lock")

    other_fd = os.open(lock_path, os.O_CREAT | os.O_RDWR)
    fcntl.flock(other_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    try:
        with pytest.raises(PipelineLockedError):
            with acquire_pipeline_lock(output_dir):
                pass
    finally:
        fcntl.flock(other_fd, fcntl.LOCK_UN)
        os.close(other_fd)


def test_lock_released_after_context_exits_even_on_exception(tmp_path):
    """即使 `with` 块内抛出异常, 锁也应在退出时正确释放, 后续调用能重新
    获取(验证 finally 释放逻辑在异常路径下同样生效)。
    """
    output_dir = str(tmp_path / "out")

    class _Boom(Exception):
        pass

    with pytest.raises(_Boom):
        with acquire_pipeline_lock(output_dir):
            raise _Boom()

    # 锁已释放, 应能再次正常获取。
    with acquire_pipeline_lock(output_dir):
        pass
