"""Python shim for the BUSHWHACKER Rust CLI.

The ``bushwhacker`` wheels bundle the compiled Rust binary. This package
exists so that ``python -m bushwhacker`` works and so downstream Python
code can locate the binary via :func:`binary_path`. All inference happens
in the Rust engine; there is no Python implementation.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import sysconfig
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("bushwhacker")
except PackageNotFoundError:  # running from a source checkout
    __version__ = "0.0.0+unknown"


def binary_path() -> str:
    """Return the path to the installed ``bushwhacker`` binary."""
    exe_name = "bushwhacker.exe" if os.name == "nt" else "bushwhacker"
    # Prefer the binary bundled into the active Python environment (venv/conda)
    # so `python -m bushwhacker` cannot accidentally pick up an unrelated binary
    # earlier on PATH.
    env_dirs: list[str] = []
    scripts_dir = sysconfig.get_path("scripts")
    if scripts_dir:
        env_dirs.append(scripts_dir)
    env_dirs.extend(
        (
            os.path.join(sys.exec_prefix, "bin"),
            os.path.join(sys.exec_prefix, "Scripts"),
        )
    )
    try:
        user_scripts = sysconfig.get_path("scripts", f"{os.name}_user")
        if user_scripts:
            env_dirs.append(user_scripts)
    except KeyError:
        pass

    seen: set[str] = set()
    for d in env_dirs:
        norm = os.path.normcase(os.path.abspath(d))
        if norm in seen:
            continue
        seen.add(norm)
        cand = os.path.join(d, exe_name)
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand

    found = shutil.which("bushwhacker")
    if found is not None:
        return found

    raise FileNotFoundError(
        "Could not locate the 'bushwhacker' binary. If you installed "
        "from source, run `cargo build --release` or reinstall with "
        "`pip install --force-reinstall bushwhacker`."
    )


def main(argv: list[str] | None = None) -> int:
    """Run the ``bushwhacker`` CLI, forwarding all arguments.

    Returns the CLI's exit code.
    """
    args = sys.argv[1:] if argv is None else argv
    return subprocess.run([binary_path(), *args], check=False).returncode
