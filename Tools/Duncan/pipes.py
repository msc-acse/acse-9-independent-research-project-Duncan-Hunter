import gmsh
import numpy as np
import os

model = gmsh.model
factory = model.occ
mesh = model.mesh

gmsh.initialize()
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.2)
def vec_angle(vec1, vec2):
    return np.arccos((np.dot(vec1, vec2)/
                      (np.linalg.norm(vec1)*np.linalg.norm(vec2))))

def rotate_vector(vec, angle, axis_of_rotation):
    vec = np.cross

class Network():
    def __init__(self, length, radius, lcar):
        self.lcar = lcar
        self.in_centre = None
        self.out_centre = None
        self.length = length
        self.direction = np.array([1, 0, 0])  # direction of pipe
        self.up_vector =  np.array([0, 0, 1])
        self.radius = radius
        #self.lcar = lcar
        
        self.inflow_tag, self.in_centre, self.out_centre = self._create_rotate_cylinder(length, self.direction)
        self._set_mesh_size(self.inflow_tag, lcar)
        self.surfaces = self._get_surfaces(self.inflow_tag)
        print(self.out_centre)
        #self.edge_surface = self.surfaces[0]
        #self.out_surface = self.surfaces[1]
        #self.in_surface = self.surfaces[2]

        self.entities = [self.inflow_tag]  # list of dimtags of entities in pipe
        self.no_entities = 1
    
    def _create_rotate_cylinder(self, length, direction):
        cyl_tag = (3, factory.addCylinder(0, 0, 0, 0, 0, length, self.radius))
        factory.synchronize()
        surfaces = model.getBoundary([cyl_tag], False, False, False)
        
        rotation_axis = np.cross(direction, self.up_vector)
        rotation_angle = np.arccos((np.dot(direction, self.up_vector)/  # can replace
                                (np.linalg.norm(direction)*np.linalg.norm(self.up_vector))))
        
        factory.rotate([cyl_tag], 0, 0, length/2, 
                    rotation_axis[0], rotation_axis[1], rotation_axis[2],
                    -rotation_angle)
        
        out_centre = factory.getCenterOfMass(2, surfaces[1][1])
        in_centre = factory.getCenterOfMass(2, surfaces[2][1])        
        factory.synchronize()
        return cyl_tag, np.array(in_centre), np.array(out_centre)
    
    def _set_mesh_size(self, dimtag, lcar):
        ov = model.getBoundary(dimtag, False, False, True)
        mesh.setSize(ov, lcar)
        
    def _get_surfaces(self, dimtag):
        """
        dimtag: (tuple) dimtag of object you're trying to get the boundary of.
        """
        surfaces = model.getBoundary(dimtag, False, False, False)
        return surfaces
        #if model.getType(dim, tag) != 'Plane':
            
    def add_pipe(self, length, lcar):
        self.no_entities += 1
        temp_tag, in_centre, out_centre = self._create_rotate_cylinder(length, self.direction)
        
        # This line may need fixing
        translate_vector = self.out_centre - in_centre #- np.array([-length/2, 0, length/2])  # 
        #print(translate_vector)    
        factory.translate([temp_tag], translate_vector[0], translate_vector[1], translate_vector[2])
        factory.synchronize()
        self._set_mesh_size(temp_tag, lcar)
        
        ## May need a bug fix
        self.out_centre = out_centre+translate_vector
        print(self.out_centre)
        
        self.entities.append(temp_tag)
        self.no_entities += 1

    def add_curve(self, new_direction, centre_of_rotation, lcar):
        """
        new_direction: (list) Direction pipe will be facing in x, y, z vector format.
        e.g. [0, 1, 0] faces positive y.
        centre_of_rotation: (list) Relative centre of rotation to the out_flow_surface in x, y, z vector format.
         e.g. [0, 1, 0] is +1 in y from the current end of pipe. Ensure this is far enough from the edge of the pipe.
        Signs of values in new_direction should match signs in centre_of_rotation.
        lcar: (float) Size of mesh in this piece.
        """
        # Check input
        if (len(new_direction) == len(centre_of_rotation) > 3):
            raise ValueError("Array given is too long (should be length 3, for x, y, z)")
        sign1 = np.sign(centre_of_rotation)
        sign2 = np.sign(new_direction)
        for i in range(3):
            if sign1[i] * sign2[i] < 0:
                raise ValueError("Sign of array values should match")
        if (np.linalg.norm(centre_of_rotation) < 2*self.radius):
            raise ValueError("Bend of radius is not large enough, place centre_of_rotation further away")
        # Create disk
        disk_tag = (2, factory.addDisk(self.out_centre[0], self.out_centre[1], self.out_centre[2], 
                                       self.radius, self.radius))
        # Align disk
        rotate_axis = np.cross(self.up_vector, self.direction)
        rotate_angle = np.arccos((np.dot(self.direction, self.up_vector)/  # can replace
                                (np.linalg.norm(self.direction)*np.linalg.norm(self.up_vector))))
        factory.rotate([disk_tag],self.out_centre[0],
                        self.out_centre[1],
                        self.out_centre[2],*list(rotate_axis), rotate_angle)
        
        # Revolve disk (Extrude through revolution)
        centre = self.out_centre + np.array(centre_of_rotation)
        new_direction = np.array(new_direction)
        angle = np.arccos((np.dot(self.direction, new_direction)/
                                (np.linalg.norm(self.direction)*np.linalg.norm(new_direction))))
        axis_of_rotation = np.cross(self.direction, new_direction)
        revolve_tags = factory.revolve([disk_tag], *list(centre), *list(axis_of_rotation), angle)
        factory.synchronize()

        #self._set_mesh_size(revolve_tags[1], lcar)
        self.direction = new_direction
        self.out_centre = np.array(factory.getCenterOfMass(2, revolve_tags[0][1]))
        self.entities.append((revolve_tags[1]))
    
    def fuse_objects(self):
        fuse_tool = self.entities[0]
        for i in range(1, len(self.entities)):
            out_dim_tags, out_dim_tags_map = factory.fuse([fuse_tool], [self.entities[i]], removeObject=True)
            print(out_dim_tags)
            fuse_tool = out_dim_tags[-1]
        #factory.fuse([self.entities[0]], self.entities[1:])
        factory.synchronize()
        surfaces = model.getBoundary(model.getEntities(3))
        for dim, tag in surfaces:
            print(dim, tag, model.getType(dim, tag))
            if model.getType(dim, tag) == "Surface of Revolution":
                self._set_mesh_size((dim, tag), 0.05)
            else:
                self._set_mesh_size((dim, tag), 0.2)

network = Network(0.1, 0.5, 0.05)
network.add_pipe(1, 0.1)
#network.add_pipe(1, 0.2)
#network.add_pipe(3, 0.07)
network.add_curve([0, -1, 0], [0, -1, 0], 0.11)
network.add_curve([-1, -1, 0], [-6, 0, 0], 0.2)
network.add_curve([0, -1, 0], [1, -1, 0], 0.2)
#network.add_pipe(1, 0.1)
#network.add_curve([1, 0, 0], [1, 0, 0], 0.1)
#network.add_pipe(1, 0.1)
#
#network.add_pipe(2, 0.2)
#network.add_curve([0, 0, 1], [0, 0, 1], 0.1)
#network.add_curve([-1, 0, 0], [-1, 0, 0], 0.1)
network.fuse_objects()

mesh.generate(3)

gmsh.write("testmsh1.msh2")

os.rename("testmsh1.msh2", "testmsh1.msh")

gmsh.option.setNumber("Mesh.SurfaceFaces", 1)
gmsh.option.setNumber("Geometry.LabelType", 2)
gmsh.option.setNumber("Geometry.SurfaceNumbers", 1)
#gmsh.option.setColor("Geometry.Color.Text", 0, 0, 0)
model.setVisibility(model.getEntities(3),0)
gmsh.fltk.run()

gmsh.finalize()