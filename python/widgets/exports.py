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
        
        self.geo_import = ""
        self.geo_export = ""
        
        self.txt_import = ""
        self.txt_export = ""

    def export_to_json(self):
        self.write_to_json_dict()
        if self.file_path is not None:
            with open(file_path, "w") as write_js:
                json.dumps(self.JSON_export, write.js)
            
    def write_to_json_dict(self):
        for particle in self.system_obj.particle_list:
            try:
                self.JSON_export["id"].append(particle.id)
                self.JSON_export["position"].append(particle.initial_pos)
                self.JSON_export["velocity"].append(particle.velocity)
                self.JSON_export["acceleration"].append(particle.acceleration)
                self.JSON_export["pressure"].append(particle.pressure_force)
                self.JSON_export["surface tension"].append(particle.surface_tension)
                self.JSON_export["buoyancy"].append(particle.buoyancy)
                self.JSON_export["viscosity"].append(particle.viscosity)
                self.JSON_export["mass"].append(particle.mass)
            except KeyError:
                self.JSON_export["id"] = [particle.id]
                self.JSON_export["position"] = [particle.initial_pos]
                self.JSON_export["velocity"] = [particle.velocity]
                self.JSON_export["acceleration"] = [particle.acceleration]
                self.JSON_export["pressure"] = [particle.pressure_force]
                self.JSON_export["surface tension"] = [particle.surface_tension]
                self.JSON_export["buoyancy"] = [particle.buoyancy]
                self.JSON_export["viscosity"] = [particle.viscosity]
                self.JSON_export["mass"] = [particle.mass]
                
    def import_from_json(self):
        if self.JSON_import is not None:
            try:
                with open(self.file_path, "r") as read_js:
                    self.JSON_import = json.loads(read_js)
            except FileNotFoundError:
                self.export_to_JSON()
        
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