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