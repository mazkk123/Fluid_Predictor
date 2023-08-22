import json
import shutil
import os
import sys

sys.path.append()

try:
    from Fluid_Utilities.system import FluidSystem
except ImportError:
    sys.stderr(f"Module {path} not found, try again")

class ExportUtils:
    
    def __init__(self,
                 file_path:str = None,
                 file_type:str = None,
                 system_obj:FluidSystem = None):
        
        if file_path is not None:
            self.file_path = file_path
        if file_type is not None:
            self.file_type = file_type
        if system_obj is not None:
            self.system_obj = system_obj
            
        self.JSON_export = {}
        self.JSON_import = None

    def export_to_json(self):
        pass

    def import_from_json(self):
        pass
    
    def import_from_houdini(self):
        pass
        
    def export_to_houdini(self):
        pass
    
    def import_from_maya(self):
        pass
        
    def export_to_maya(self):
        pass
    
    def export_to_txt(self):
        pass

    def import_from_txt(self):
        pass