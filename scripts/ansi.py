#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ansi.py
=======

Minimal colour support for the console output. No dependencies.

Colour is switched off automatically when the output is not a terminal - that
is, when it is redirected into a file or passed through a pipe. Otherwise
control sequences such as ``\\x1b[32m`` end up in the log files.

The module also respects the ``NO_COLOR`` convention (see
https://no-color.org): if that environment variable is set, everything stays
plain. ``FORCE_COLOR`` forces the opposite.

On Windows the processing of ANSI sequences has to be enabled once
(``ENABLE_VIRTUAL_TERMINAL_PROCESSING``); ``_enable_windows_vt`` does that
through the Win32 API, without colorama as an extra package.
"""

from __future__ import annotations

import os
import re
import sys

RESET = "\033[0m"
BOLD = "\033[1m"
GREEN = "\033[32m"
CYAN = "\033[36m"
YELLOW = "\033[33m"

def _enable_windows_vt() -> bool:
    """Enables ANSI processing in the Windows console."""
    try:
        import ctypes

        kernel32 = ctypes.windll.kernel32
        handle = kernel32.GetStdHandle(-11)          # STD_OUTPUT_HANDLE
        mode = ctypes.c_uint32()
        if not kernel32.GetConsoleMode(handle, ctypes.byref(mode)):
            return False
        ENABLE_VIRTUAL_TERMINAL_PROCESSING = 0x0004
        return bool(kernel32.SetConsoleMode(
            handle, mode.value | ENABLE_VIRTUAL_TERMINAL_PROCESSING))
    except Exception:
        return False

def _supported() -> bool:
    if os.environ.get("NO_COLOR"):
        return False
    if os.environ.get("FORCE_COLOR"):
        return True
    if not (hasattr(sys.stdout, "isatty") and sys.stdout.isatty()):
        return False
    if os.name == "nt":
        return _enable_windows_vt()
    return True

ENABLED = _supported()

def disable():
    """Switches colour output off at run time (for --no-color)."""
    global ENABLED
    ENABLED = False

def paint(text, code) -> str:
    """Colours ``text``, or returns it unchanged."""
    if not ENABLED or not text:
        return text
    return f"{code}{text}{RESET}"

def green(text):
    return paint(text, GREEN)

def cyan(text):
    return paint(text, CYAN)

def yellow(text):
    return paint(text, YELLOW)

def bold(text):
    return paint(text, BOLD)

# ----------------------------------------------------------------------------
# Chemistry-specific
# ----------------------------------------------------------------------------

HALOGEN_SYMBOLS = ("F", "Cl", "Br", "I", "At")

_LABEL = re.compile(r"^([A-Za-z]+)(\d*)$")

def atom_label(label: str) -> str:
    """Colours the element symbol of an atom label if it is a halogen.

    ``"Cl12"`` -> cyan ``Cl`` plus plain ``12``. The running number stays
    uncoloured, so that the symbol catches the eye and not the index.
    Non-halogens such as ``"H5"`` are returned unchanged.
    """
    if not label:
        return label
    m = _LABEL.match(label)
    if not m:
        return label
    symbol, number = m.group(1), m.group(2)
    if symbol in HALOGEN_SYMBOLS:
        return cyan(symbol) + number
    return label

def element(symbol: str) -> str:
    """Colours a bare element symbol if it is a halogen."""
    return cyan(symbol) if symbol in HALOGEN_SYMBOLS else symbol
