import ctypes

# List of DLLs to check
dll_names = ["IlmThread-3_1_d.dll", "Imath-3_1_d.dll", "OpenEXR-3_1_d.dll"]

def check_dll_loading(dll_name):
    try:
        ctypes.CDLL(dll_name)
        return True
    except OSError:
        return False

for dll_name in dll_names:
    loaded = check_dll_loading(dll_name)
    if loaded:
        print(f"{dll_name} is loaded successfully.")
    else:
        print(f"{dll_name} failed to load.")