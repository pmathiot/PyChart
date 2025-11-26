"""
log.py

This module provides logging utilities for the PyChart application. It includes functions for enabling debug mode, logging informational messages, and logging debug messages.

Functions:
- set_debug(enabled: bool): Enables or disables debug mode.
- info(msg: str): Logs an informational message.
- debug(msg: str): Logs a debug message if debug mode is enabled.
"""

# pychart/log.py
DEBUG_ENABLED = False

def set_debug(enabled: bool):
    global DEBUG_ENABLED
    DEBUG_ENABLED = enabled
    if DEBUG_ENABLED:
        print("[DEBUG] Debug mode enabled.")

def info(msg: str):
    print(f"[PyChart] {msg}")

def debug(msg: str):
    if DEBUG_ENABLED:
        print(f"[DEBUG] {msg}")