import json
import os
import sys
import pyfbx

class ExportUtils:
    
    def __init__(self,
                 parent_dir:str = None,
                 file_path:str = "Particle",
                 export_prefix:str = "particle",
                 frame_padding:int = 4,
                 file_type:str = "geo",
                 frame_number:int = 0,
                 system_obj = None):
        
        if parent_dir is not None:
            self.parent_dir = parent_dir
        if file_path is not None:
            self.file_path = file_path
        if file_type is not None:
            self.file_type = file_type
        if frame_number is not None:
            self.frame_number = frame_number
        if export_prefix is not None:
            self.export_prefix = export_prefix
        if frame_padding is not None:
            self.frame_padding = frame_padding
        if system_obj is not None:
            self.system_obj = system_obj
            
        self.JSON_export = {}
        self.JSON_import = {}
        
        self.geo_import = ""
        self.geo_export = ""
        self.geo_data_frame = []

        self.output_suffixes = ["geo", "json", "abc", "fbx", "txt"]
        
        self.txt_export = ""
        self.particle_data = """
            Frame {0}
            Position {1}
            Velocity {2}
            Acceleration {3}
            Pressure {4}
            Id {5}
            \n
        """
        self.txt_import = ""

        self.create_directories()

    def clean_geo_data(self, frame_number):
        try:
            pass
        except OSError:
            pass

    def export_data(self):
        if self.file_type is not None:

            if self.file_type == "json":
                self.export_to_json()
            elif self.file_type == "geo":
                self.export_to_houdini(self.frame_number)
                print(f"exported {self.frame_number} to geo ...")
            elif self.file_type == "txt":
                self.export_to_txt()
            elif self.file_type == "fbx":
                self.export_as_fbx()

    def export_to_json(self):
        self.write_to_json_dict()
        if self.file_path is not None:
            with open(self.file_path, "w") as write_js:
                write_js.write(json.dumps(self.JSON_export))
            
    def write_to_json_dict(self):
        for frame in range(self.system_obj.num_frames):
            self.JSON_export.setdefault(frame, {})
            for particle in self.system_obj.particle_list:
                self.JSON_export[frame].setdefault("id", []).append(particle.id.copy())
                self.JSON_export[frame].setdefault("position", []).append(particle.initial_pos.copy())
                self.JSON_export[frame].setdefault("velocity", []).append(particle.velocity.copy())
                self.JSON_export[frame].setdefault("acceleration", []).append(particle.acceleration.copy())
                self.JSON_export[frame].setdefault("pressure force", []).append(particle.pressure_force.copy())
                self.JSON_export[frame].setdefault("surface tension", []).append(particle.surface_tension.copy())
                self.JSON_export[frame].setdefault("buoyancy", []).append(particle.buoyancy.copy())
                self.JSON_export[frame].setdefault("viscosity", []).append(particle.viscosity.copy())
                self.JSON_export[frame].setdefault("mass", []).append(particle.mass.copy())
                    
    def import_from_json(self):
        if self.JSON_import is not None:
            try:
                with open(self.file_path, "r") as read_js:
                    JSON_import = json.loads(read_js)
                
                for frame in range(self.system_obj.num_frames):
                    for i in range(len(self.system_obj.num_particles)):
                        try:
                            self.JSON_import.setdefault(frame, []).append(JSON_import[frame]["position"][i])
                        except KeyError:
                            print(f"Key {JSON_import[frame]} not found, try again")

            except FileNotFoundError:
                self.export_to_json()
        
    def import_from_houdini(self):
        pass

    def create_directories(self):
        
        self.parent_dir = os.getcwd()
        self.dir_title = self.file_path

        try:
            self.new_dir = os.path.join(self.parent_dir, self.dir_title)

            try:    
                os.mkdir(self.new_dir)

                print(f"created directory: {self.new_dir}")

                os.chdir(self.new_dir)

            except FileExistsError:
                print(f"Directory {self.new_dir} already exists!")
                os.chdir(self.new_dir)
                return

        except ValueError:
            print(f"Directory name is not of type {type(str)}, try again")
            return

    def write_to_geo(self, frame_number):
        if self.geo_export is not None:
            
            self.geo_data_frame = []
            self.geo_export = ""  # Initialize GEO data for this frame

            self.geo_export += "PGEOMETRY V5\n"
            self.geo_export += "NPoints " + str(self.system_obj.num_particles) + " NPrims 1\n"
            self.geo_export += "NPointGroups 0 NPrimGroups 0\n"
            self.geo_export += "NPointAttrib 2 NVertexAttrib 0 NPrimAttrib 1 NAttrib 0\n"
            self.geo_export += "PointAttrib\n"
            self.geo_export += "Cd 3 float 1 1 1\n"
            self.geo_export += "pscale 1 float 1\n"

            for particle in self.system_obj.particle_list:
                self.geo_export +=  f"{particle.initial_pos[0]} {particle.initial_pos[1]} {particle.initial_pos[2]} 1 ("
                self.geo_export += f"{particle.colour[0]} {particle.colour[1]} {particle.colour[2]} {particle.size})\n"

            self.geo_export += "PrimitiveAttrib\n"
            self.geo_export += "generator 1 index 1 papi\n"
            
            self.geo_export += f"Part {self.system_obj.num_particles}\n"
            for i in range(self.system_obj.num_particles):
                self.geo_export += f"{i} "
            self.geo_export += f"[{frame_number}]\n"

            self.geo_export += "beginExtra\n"
            self.geo_export += "endExtra\n"

            self.geo_data_frame.append(self.geo_export)

    def export_to_houdini(self, frame_number):

        self.write_to_geo(frame_number)
            
        all_geo_data = "\n".join(self.geo_data_frame)

        self.output_geo_file = f"{self.export_prefix}.{frame_number:0{self.frame_padding}d}.{self.output_suffixes[0]}"

        with open(self.output_geo_file, 'w') as geo_file:
            geo_file.write(all_geo_data)

    def export_as_fbx(self):
        pass

    def export_as_alembic(self):
        pass
    
    def write_to_txt(self):

        for frame in self.system_obj.num_frames:
            for particle in self.system_obj.particle_list:
                self.particle_data.format(frame, 
                                          particle.initial_pos,
                                          particle.velocity,
                                          particle.acceleration,
                                          particle.pressure,
                                          particle.id
                                        )
                self.txt_export += self.particle_data

    def export_to_txt(self):
        
        self.write_to_txt()

        with open(self.output_txt_name, "w") as output_file:
            output_file.write(self.txt_export)

    def import_from_txt(self):
        
        try:
            with open(self.output_txt_name, "r") as read_file:
                self.txt_import = read_file.read()

            self.parse_txt()
                
        except FileNotFoundError:
            print(f"file {self.output_txt_name} doesn't exist")

    def parse_txt(self):

        if len(self.txt_import) > 0:
            for attribs in self.txt_import.split("\n"):
                self.positions, self.frames = [], []
                for attrib in attribs.split(" "):
                    if attrib.startswith("F"):
                        self.frames.append(attrib)
                    elif attrib.startswith("P"):
                        self.positions.append(attrib)