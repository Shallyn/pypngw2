import ctypes
import sys
from pathlib import Path


if sys.platform == "darwin":
    SHARED_LIBRARY_SUFFIX = ".dylib"
elif sys.platform == "win32":
    SHARED_LIBRARY_SUFFIX = ".dll"
else:
    SHARED_LIBRARY_SUFFIX = ".so"


def shared_library_name(name):
    return f"{name}{SHARED_LIBRARY_SUFFIX}"


def shared_library_path(package_dir, subdir, name):
    path = Path(package_dir) / subdir / shared_library_name(name)
    if path.exists():
        return path

    available = sorted(
        str(p.relative_to(package_dir))
        for p in (Path(package_dir) / subdir).glob(f"{name}.*")
    )
    raise FileNotFoundError(
        f"Cannot find native library {shared_library_name(name)} in {path.parent}. "
        f"Available candidates: {available or 'none'}"
    )


def load_shared_library(package_dir, subdir, name, global_symbols=True):
    mode = ctypes.DEFAULT_MODE
    if global_symbols:
        mode |= getattr(ctypes, "RTLD_GLOBAL", 0)
    return ctypes.CDLL(str(shared_library_path(package_dir, subdir, name)), mode=mode)
