#!/usr/bin/env python3
"""Smoke-test the installed CMake package with a clean consumer project."""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run(command: list[str], cwd: Path, env: dict[str, str]) -> None:
    subprocess.run(command, cwd=cwd, env=env, check=True)


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="pb11-install-") as td:
        base = Path(td)
        build = base / "build"
        prefix = base / "prefix"
        consumer = base / "consumer"
        consumer.mkdir()
        (consumer / "CMakeLists.txt").write_text(
            "cmake_minimum_required(VERSION 3.16)\n"
            "project(pb11_consumer LANGUAGES CXX)\n"
            "find_package(pb11_reaction 1.0 CONFIG REQUIRED)\n"
            "add_executable(consumer main.cpp)\n"
            "target_link_libraries(consumer PRIVATE pb11::pb11)\n"
        )
        (consumer / "main.cpp").write_text(
            '#include "pb11_c.h"\n'
            '#include <cmath>\n'
            'int main() { double r = 0.0;\n'
            '  if (pb11_c_abi_version() != 1) return 1;\n'
            '  if (pb11_c_reactivity_fast(100.0, &r) != 0 || !std::isfinite(r) || r < 0.0) return 2;\n'
            '  return 0; }\n'
        )
        env = os.environ.copy()
        env.setdefault("CMAKE_BUILD_PARALLEL_LEVEL", "1")
        run(["cmake", "-S", str(ROOT), "-B", str(build), "-DPB11_BUILD_FORTRAN=OFF"], ROOT, env)
        run(["cmake", "--build", str(build), "--target", "pb11", "-j1"], ROOT, env)
        run(["cmake", "--install", str(build), "--prefix", str(prefix)], ROOT, env)
        consumer_build = base / "consumer-build"
        run([
            "cmake", "-S", str(consumer), "-B", str(consumer_build),
            f"-DCMAKE_PREFIX_PATH={prefix}", "-DCMAKE_BUILD_TYPE=Release",
        ], consumer, env)
        run(["cmake", "--build", str(consumer_build), "-j1"], consumer, env)
        run([str(consumer_build / "consumer")], consumer, env)
        print("installed CMake package consumer: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
