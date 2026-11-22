# -*- mode: python ; coding: utf-8 -*-
#
# PyInstaller build file for Graph State Orbit Explorer.
#
# Build from the project root (the folder that contains src/) with
#
#     pip install pyinstaller
#     pyinstaller graphCreator.spec
#
# The application is written to dist/GraphStateOrbitExplorer/. Build on each
# target system separately: PyInstaller does not cross-compile.

import os

ROOT = os.path.abspath(SPECPATH)

a = Analysis(
    [os.path.join(ROOT, "src", "graphCreator.py")],
    # The project root must be on the search path so that the `src.` imports resolve.
    pathex=[ROOT],
    binaries=[],
    datas=[],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    # Large optional packages that the application does not use.
    excludes=["matplotlib", "scipy", "pandas", "IPython", "pytest"],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="GraphStateOrbitExplorer",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    # Keep the console so that error messages remain visible.
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name="GraphStateOrbitExplorer",
)
