import gmsh
import numpy as np
import os
import itertools
from scipy.spatial.transform import Rotation as R
import time

model = gmsh.model
factory = model.occ
mesh = model.mesh

#gmsh.option.setNumber("General.Terminal", 1)

def vec_angle(vec1, vec2):
    return np.arccos((np.dot(vec1, vec2)/
                      (np.linalg.norm(vec1)*np.linalg.norm(vec2))))

def proj(vec1, vec2):
    """
    Component of vec1 along vec2
    vec1, vec2: (np.array shape 3) xyz vector
    
    returns: (float) component of vec1 along vec2
    """
    return np.dot(vec1, vec2)/np.linalg.norm(vec2)

def check_intersect(objects, tools):
    """
    Check if entities in tools intersect entities in objects.
    objects: (list) 1st group of entities
    tools: (list) 2nd group of entities
    
    returns: (bool) True if there is an intersection.
                    False if there is no intersection.
    """
    intersect = factory.intersect(objects, tools, removeObject=False, removeTool=False)
    if len(intersect):
        factory.remove(intersect, True)
        return True
    else:
        return False

class Network():
    """
    Overall BUGS
        First Rotation is wrong for mitred and junction
    Todo
        - Junctions
            Implemented as add_junction
            Needs to be ready for when we add pipes to it
        - Keep track of surface tags in general
            Update - made easier by classes
            Need list of centres
        - Keep track of orientations ??
            Update - made easier by classes
        - Keep track of locations
            Update - made easier by classes
            append to vol_centres
            Use dictionary to map tags to centres?
        - Mesh size of pieces after fuse
            Centre to mesh size mapping
            Fields with locations
        - Locations for postprocessing
        - Physical groups
            Need to deal with excess outflows or inflows
        - Testing
    """
    def __init__(self, length, radius, direction, lcar):
        """
        length: (float) length of piece.
        radius: (float) radius of pipe.
        direction: (list) Direction pipe will be facing in x, y, z vector format.
        lcar: (float) Mesh size of this piece
        """
        
        self.lcar = lcar
        self.length = length
        self.direction = np.array(direction)  # direction of pipe
        self.up_vector =  np.array([0, 0, 1])
        self.radius = radius
        #self.lcar = lcar
        
        piece = Cylinder(length, radius, self.direction, lcar)        
        
        self.in_surface = piece.in_surface  # unsure of use
        self.in_centre = piece.in_centre
        self.out_centre = piece.out_centre
        self.vol_tag = None
        
        # self._set_mesh_size(self.in_surface, lcar)

        self.entities = [piece.vol_tag]  # list of dimtags of entities in pipe
        self.lcar_list = [lcar]
        self.vol_centres = [piece.vol_centre]
        
        self._has_physical_groups = False
        self.physical_in_surface = None
        self.physical_out_surface = None
        self.physical_no_slip = None
        self.physical_volume = None
    
    def _set_mesh_size(self, dimtag, lcar):
        """
        Use centres or faces of objects to set desired mesh size.
        
        Box field?
        
        
        """
        
        ov = model.getBoundary(dimtag, False, False, True)
        mesh.setSize(ov, lcar)
        
    def _get_surfaces(self, dimtag):
        """
        dimtag: (tuple) dimtag of object you're trying to get the boundary of.
        Get in/out surfaces of object
        """
        surfaces = model.getBoundary(dimtag, False, False, False)
        return surfaces
        #if model.getType(dim, tag) != 'Plane':
            
    def add_pipe(self, length, lcar):
        """
        length: (float) length of pipe.
        lcar: (float) mesh size of piece.
        """
        piece = Cylinder(length, self.radius, self.direction, lcar)

        translate_vector = self.out_centre - piece.in_centre
        factory.translate([piece.vol_tag], *list(translate_vector))
        factory.synchronize()
        piece.update_centres()

        # self._set_mesh_size(piece.vol_tag, lcar)
        self.out_centre = piece.out_centre      
        self.entities.append(piece.vol_tag)

    def add_curve(self, new_direction, centre_of_rotation, lcar):
        """
        new_direction: (list) Direction pipe will be facing in x, y, z vector format.
        e.g. [0, 1, 0] faces positive y.
        centre_of_rotation: (list) Relative centre of rotation to the out_flow_surface in x, y, z vector format.
         e.g. [0, 1, 0] is +1 in y from the current end of pipe. Ensure this is far enough from the edge of the pipe.
        Signs of values in new_direction should match signs in centre_of_rotation.
        lcar: (float) Size of mesh in this piece.
        
        Todo:
            - add check for angle sign (dot product thing in add_mitered)
        """
        # Check input
        if (len(new_direction) == len(centre_of_rotation) > 3):
            raise ValueError("Array given is too long (should be length 3, for x, y, z)")
        sign1 = np.sign(centre_of_rotation)
        sign2 = np.sign(new_direction)
        for i in range(3):
            if sign1[i] * sign2[i] < 0:
                raise ValueError("Sign of array values should match")
        if (np.linalg.norm(centre_of_rotation) < 1.1* self.radius):
            raise ValueError("Bend of radius is not large enough, place centre_of_rotation further away")
        
        piece = Curve(self.radius, self.direction, new_direction, centre_of_rotation, lcar)
        translate_vector = self.out_centre - piece.in_centre
        factory.translate([piece.vol_tag], *list(translate_vector))
        #piece.out_centre += translate_vector
        factory.synchronize()
        piece.update_centres()

        #self._set_mesh_size(revolve_tags[1], lcar)
        self.direction = piece.out_direction
        self.out_centre = piece.out_centre  # np.array(factory.getCenterOfMass(2, revolve_tags[0][1]))
        self.entities.append(piece.vol_tag)  #(revolve_tags[1]))
        
        
    def add_mitered(self, new_direction, lcar):
        # Create the piece
        """
        Add a mitered bend to the pipe
        new_direction: (list, length 3) xyz vector representing the new direction of the pipe
        lcar: (float) size of mesh of this piece
        
        BUGS:
            - Two in a row non-90 degree bends doesn't behave as expected.
                Suspected cause: Applying rotation to v2
            - Contains planes when fused - is translation at the right place?
            
        """
        piece = Mitered(self.radius, self.direction, new_direction, lcar)
        translate_vector = self.out_centre - piece.in_centre
        factory.translate([piece.vol_tag], *list(translate_vector))
        piece.update_centres()
        factory.synchronize()
        self.out_centre = piece.out_centre
        self.direction = piece.out_direction
        self.entities.append(piece.vol_tag)
        
    def change_radius(self, length, new_radius, change_length, lcar):
        """
        Cylinder piece with change of radius of the pipe to new_radius.
        length: (float) Length of the piece.
        new_radius: (float) radius to change to.
        change_length: (float) Length that the change takes place over. Must be less than
                        length. Longer means a more gentle change.
        lcar: (float) Mesh size for this piece.
        """
        assert(change_length < length)
        assert(change_length > 0)
        if new_radius > self.radius:
            piece = Cylinder(length, new_radius, self.direction, lcar)
        else:
            piece = Cylinder(length, self.radius, self.direction, lcar)
        translate_vector = self.out_centre - piece.in_centre
        factory.translate([piece.vol_tag], *list(translate_vector))
        piece.update_centres()
        factory.synchronize()
        if new_radius > self.radius:
            lines = model.getBoundary([piece.in_surface], False, False)
            factory.chamfer([piece.vol_tag[1]], [lines[0][1]], [piece.in_surface[1]], [new_radius-self.radius, change_length])
            factory.synchronize()
        elif new_radius < self.radius:
            lines = model.getBoundary(piece.out_surface, False, False)
            factory.chamfer([piece.vol_tag[1]], [lines[0][1]], [piece.out_surface[1]], [self.radius-new_radius, change_length])
            factory.synchronize()
        self.out_centre = piece.out_centre
        self.radius = new_radius
        self.entities.append(piece.vol_tag)
        
    def add_T_junction(self, T_direction):
        """
        Add a T junction. T_Direction cannot be the same as self.direction 
        """
        piece = T_junction(self.radius, self.direction, T_direction)
        translate_vector = self.out_centre - piece.in_centre
        factory.translate([piece.vol_tag], *list(translate_vector))
        piece.update_centres()
        factory.synchronize()
        self.out_centre = piece.out_centre
        self.direction = piece.out_direction
        self.entities.append(piece.vol_tag)
    
    def fuse_objects(self):
        #fuse_tool = self.entities[0]
        """
        Fuse objects together.
        Find out which surface is inflow/outflow using norms.
        Update inflow_tag, outflow_tag
        Update add functions to find out which is outflow tag.
        Set mesh sizes
        
        Todo
            - Decide on what attributes we want to use
        """
        out_dim_tags, out_dim_tag_map = factory.fuse([self.entities[0]],self.entities[1:])
        factory.synchronize()
        self.vol_tag = out_dim_tags[0]
        surfaces = model.getBoundary([self.vol_tag], False)
        n_surfaces = len(surfaces)
        planes = [surfaces[i] for i in range(n_surfaces) if model.getType(*surfaces[i]) == "Plane"]
        no_slip = [surfaces[i] for i in range(n_surfaces) if model.getType(*surfaces[i]) != "Plane"]
        no_slip_tags = [dimtag[1] for dimtag in no_slip]
        found_in = False  # being extra safe
        found_out = False
        # Find in surface. use np.isclose
        for dim, tag in planes:
            loc = np.array(factory.getCenterOfMass(dim, tag))
            if np.allclose(loc, self.in_centre) and not found_in:
                self.in_surface = (dim, tag)
                found_in = True
            elif np.allclose(loc, self.out_centre) and not found_out:
                self.out_surface = (dim, tag)
                self.out_centre = loc
                found_out = True
        
        if self._has_physical_groups:
            model.removePhysicalGroups()
        
        self.physical_in_surface = model.addPhysicalGroup(2, [self.in_surface[1]])
        self.physical_out_surface = model.addPhysicalGroup(2, [self.out_surface[1]])
        self.physical_no_slip = model.addPhysicalGroup(2, no_slip_tags)
        self.physical_volume = model.addPhysicalGroup(3, [self.vol_tag[1]])
        self._has_physical_groups = True
        self.entities = out_dim_tags
        # Find out surface
                 
    def rotate_network(self, new_direction, out_flow_direction=True):
        """
        rotate pipe so outflow is pointed in new_direction
        update position of where in_centre and out_centre is
        """


class PipePiece():
    def __init__(self, radius, vol_tag, in_tag, out_tag,
                 in_direction, out_direction):
        """
        radius: (float) radius of piece
        vol_tag: (tuple, length 2) GMSH dimtag (dimension, tag) representing volume
        in_tag: (tuple, length 2) GMSH dimtag representing in surface
        out_tag: (tuple, length 2) GMSH dimtag representing out surface
        in_direction: (np array, shape 3) xyz vector representing direction going in
        out_direction: (np array, shape 3) xyz vector representing direction going out
        """
        self.radius = radius
        self.vol_tag  = vol_tag
        self.in_surface = in_tag
        self.out_surface = out_tag
        self.vol_centre = np.array(factory.getCenterOfMass(*vol_tag))
        self.in_centre = np.array(factory.getCenterOfMass(*in_tag))
        self.out_centre = factory.getCenterOfMass(*out_tag)
        self.in_direction = in_direction
        self.out_direction = out_direction
        
    def update_centres(self):
        self.vol_centre = np.array(factory.getCenterOfMass(*self.vol_tag))
        self.in_centre = np.array(factory.getCenterOfMass(*self.in_surface))
        self.out_centre = np.array(factory.getCenterOfMass(*self.out_surface))     


class Cylinder(PipePiece):
    """
    Class representing a GMSH cylinder with base at 0,0,0 facing upwards.
    """
    def __init__(self, length, radius, direction, lcar):
        """
        length: (float) length of cylinder
        radius: (float) radius of cylinder
        direction: (list, length 3) xyz vector representing direction cylinder is facing
        lcar: (float) Mesh size for this piece
        """
        self.length = length
        self.lcar = lcar
        vol_tag = (3, factory.addCylinder(0,0,0,0,0,length,radius))
        factory.synchronize()
        surfaces = model.getBoundary([vol_tag], False)
        in_surface = surfaces[2]  # these are switched
        out_surface = surfaces[1]
        direction = np.array(direction)
        up_vector = np.array([0, 0, 1])

        if np.allclose(direction, up_vector) is False:
            cross = np.cross(up_vector, direction)
            rotation_axis = np.cross(direction, up_vector)
            if np.allclose(rotation_axis, np.array([0, 0, 0])):  # if rotating 180 degrees
                rotation_axis = np.array([1, 0, 0])
            rotation_angle = np.arccos((np.dot(direction, up_vector)/  # can replace
                                    (np.linalg.norm(direction)*np.linalg.norm(up_vector))))
            
            if np.dot(rotation_axis, cross) > 0:
                rotation_angle *= -1
            factory.rotate([vol_tag], 0, 0, length/2, 
                        rotation_axis[0], rotation_axis[1], rotation_axis[2],
                        -rotation_angle)
            factory.synchronize()
            out_direction = np.copy(direction)
            in_direction = np.copy(direction)
        else:
            out_direction = np.copy(up_vector)
            in_direction = np.copy(up_vector)
            
        super(Cylinder, self).__init__(radius, vol_tag, in_surface, out_surface, in_direction, out_direction)

class Curve(PipePiece):
    """
    Class representing a GMSH curve by revolution
    """
    def __init__(self, radius, in_direction, out_direction, centre_of_rotation,
                 lcar):
        """
        radius: (float) radius of the pipe
        in_direction: (list, length 3) xyz vector representing direction going in
        out_direction: (list, length 3) xyz vector representing direction going out
        centre_of_rotation: (list, length 3) xyz coordinates of centre of rotation (distance from 0,0,0)
                            Signs of out_direction and centre_of rotation need to match to get sensible
                            answers. e.g. [0,1,0] and [0,1,0] get a 90 degree +ve bend to +y.
                            e.g. [0,0,1] and [0,1,1]
        lcar: (float) Mesh size for this piece
        """
        in_tag = (2, factory.addDisk(0,0,0, radius, radius))
        up_vector = np.array([0,0,1])
        disk_centre = factory.getCenterOfMass(*in_tag)
        in_direction = np.array(in_direction)
        out_direction = np.array(out_direction)
        if np.allclose(in_direction, np.array([0, 0, 1])) is False:  # only have to rotate if its not facing up
            rotate_axis = np.cross(up_vector, in_direction)
            if np.allclose(rotate_axis, np.zeros(3)):
                rotate_axis = np.array([1, 0, 0])
            rotate_angle = np.arccos((np.dot(in_direction, up_vector)/  # can replace
                                    (np.linalg.norm(in_direction)*np.linalg.norm(up_vector))))
            
            factory.rotate([in_tag],*disk_centre,*list(rotate_axis), rotate_angle)
        revolve_axis = np.cross(out_direction, in_direction)
        revolve_angle = vec_angle(in_direction, out_direction)
        revolve_tags = factory.revolve([in_tag], *centre_of_rotation, *list(revolve_axis), -revolve_angle)
        factory.synchronize()
        out_tag = revolve_tags[0]
        vol_tag = revolve_tags[1]
        super(Curve, self).__init__(radius, vol_tag, in_tag, out_tag, in_direction, out_direction)
        

class Mitered(PipePiece):
    """
    Piece creation is done by the painful process of masking (intersect)
        a cylinder with a chamfered box. The piece is then mirrored, rotated
        and fused.
        
        The piece is then rotated to face the direction of the outflow.
        It is then rotated about the direction of outflow to match the new direction
    """
    def __init__(self, radius, in_direction, out_direction, lcar):
        """
        radius: (float) radius of the pipe
        in_direction: (list, length 3) xyz vector representing direction going in
        out_direction: (list, length 3) xyz vector representing direction going out
        lcar: (float) Mesh size for this piece
        
        BUGS:
            - If in_direction is not in direction of x, y, or z, may get wrong out_direction
        """
        in_direction = np.array(in_direction)
        out_direction = np.array(out_direction)  # clean up v's
        up_direction = np.array([0, 0, 1])
        v1 = -up_direction  # original inlet direction
        v3 = np.array(in_direction)
        v4 = np.array(out_direction)
        height = 2*radius
        # Chamfer cylinder
        angle = vec_angle(v4, v3)
        v2 = np.array([np.sin(np.pi-angle), 0, -np.cos(np.pi-angle)])  # original outlet direction in xz plane
        cyl1 = (3, factory.addCylinder(0,0,0,0,0,height, radius))  # create cylinder
        box1 = (3, factory.addBox(-radius-1,-radius,-1, 2*radius+1, 2*radius, height+1))  # create box
        factory.synchronize()
        surface = model.getBoundary([box1], False, False)[5]  # Top surface to chamfer from
        line = model.getBoundary([surface], False, False)[2]  # -x line on top
        sdist = 2*radius*np.tan(angle/2)  # distances to chamfer to get correct angle
        #ldist = 2*radius*np.cos(angle/2)
        factory.chamfer([box1[1]], [line[1]], [surface[1]], [2*radius, sdist])  # chamfer
        int_tag = factory.intersect([box1], [cyl1])  # intersect (chamfer on cylinder)
        fuse = factory.fuse([int_tag[0]], int_tag[1:])[0]  # fuse
        fuse2 = factory.copy([fuse])  # create copy and mirror
        factory.symmetrize([fuse2], 1, 0, 0, 0)
        factory.synchronize()
        surface = model.getBoundary([fuse2], False, False)[1]  # get center of rotation
        com = factory.getCenterOfMass(*surface)
        factory.rotate([fuse2], *com, 0, 1, 0, -(np.pi-angle))  # rotate to create piece
        vol_tag = factory.fuse([fuse], [fuse2])[0][0]  # fuse
        factory.synchronize()

        surfaces = model.getBoundary([vol_tag], False, False)
        in_surface = surfaces[3]  # bottom
        out_surface = surfaces[6]  # angled
        
        # Rot1: rotate object so inlet faces correct direction
        if np.allclose(in_direction, np.array([0, 0, 1])) is False:
            rot_axis = np.cross(v1, v3)
            if np.allclose(rot_axis, np.zeros(3)):
                rot_axis = np.array([1, 0, 0])
            rot_angle = vec_angle(v1, v3)
            in_centre = factory.getCenterOfMass(*in_surface)
            factory.rotate([vol_tag], *in_centre, *list(rot_axis), -rot_angle)
            factory.synchronize()
            rot_vec = -rot_angle * rot_axis  # minus sign here
            rot1 = R.from_rotvec(rot_vec)
            v2 = rot1.apply(v2)  # find new outlet vector post rotation
        else:
            in_centre = factory.getCenterOfMass(*in_surface)

        #Rot2
        b1 = np.cross(v4, v3)  # basis vectors perpendicular to direction (rotation axis)
        b2 = np.cross(b1, v3)  # and perpendicular to other basis
        alpha = np.array([proj(v2, b1), proj(v2, b2)]) # v2 in basis
        beta = np.array([proj(v4, b1), proj(v4, b2)])  # v4 (new_direction) in basis
        # Find angle between two vectors in basis
        rot2_angle = vec_angle(alpha, beta)
        cross = np.cross(v2, v4)
        if np.dot(v3, cross) < 0:
            rot2_angle *= -1
        

        factory.rotate([vol_tag], *in_centre, *list(v3), rot2_angle)
        factory.synchronize()

        super(Mitered, self).__init__(radius, vol_tag, in_surface, out_surface, v3, v4)


class T_junction(PipePiece):
    def __init__(self, radius, direction, T_direction, T_type='in'):
        self.out_num=2
        direction = np.array(direction)
        self.T_direction = np.array(T_direction)
        

        
        T_angle = vec_angle(direction, self.T_direction)
        if T_angle < np.pi/2:
            rot_sign = -1
            beta = np.pi/2 - abs(T_angle)
        else:
            rot_sign = 1
            beta = abs(T_angle) - np.pi/2
            
        height = radius*np.tan(beta) + radius/np.cos(beta)
        # Calculating height needed to emerge from merge

        in_tag = (3, factory.addCylinder(0,0,0,0,0,1.1*height, radius))
        mid_tag = (3, factory.addCylinder(0,0,0,1.1*height,0,0,radius))
        out_tag = (3, factory.addCylinder(0,0,0,0,0,-1.1*height, radius))
        factory.rotate([mid_tag], 0,0,0, 0, 1, 0, rot_sign*beta)
        factory.synchronize()
        
        vol_tags, out_dim_tag_map = factory.fuse([in_tag], [mid_tag, out_tag])
        vol_tag = vol_tags[0]
        factory.synchronize()
        
        surfaces = model.getBoundary([vol_tag], False)
        in_surface = surfaces[5]  # 5 3
        out_surface = surfaces[3]
        self.mid_surface = surfaces[4]

        mid_direction = np.array([1, 0, 0])

        up_vector = np.array([0, 0, 1])
        if np.allclose(direction, up_vector) is False:
            cross = np.cross(up_vector, direction)
            rotation_axis = np.cross(direction, up_vector)
            if np.allclose(rotation_axis, np.array([0, 0, 0])):  # if rotating 180 degrees
                rotation_axis = np.array([1, 0, 0])
            rotation_angle = np.arccos((np.dot(direction, up_vector)/  # can replace
                                    (np.linalg.norm(direction)*np.linalg.norm(up_vector))))
            
            if np.dot(rotation_axis, cross) > 0:
                rotation_angle *= -1
            vol_centre = factory.getCenterOfMass(*vol_tag)
            factory.rotate([vol_tag], *vol_centre, 
                        rotation_axis[0], rotation_axis[1], rotation_axis[2],
                        -rotation_angle)
            factory.synchronize()
            
            rot_vec = -rotation_angle * rotation_axis
            rot1 = R.from_rotvec(rot_vec)
            mid_direction = rot1.apply(mid_direction)
            print(mid_direction)
            
            out_direction = np.copy(direction)
            in_direction = np.copy(direction)
        else:
            out_direction = np.copy(up_vector)
            in_direction = np.copy(up_vector)

        b1 = np.cross(T_direction, direction)  # basis vectors perpendicular to direction (rotation axis)
        b2 = np.cross(b1, direction)  # and perpendicular to other basis
        alpha = np.array([proj(mid_direction, b1), proj(mid_direction, b2)]) # v2 in basis
        beta = np.array([proj(T_direction, b1), proj(T_direction, b2)])  # v4 (new_direction) in basis
        # Find angle between two vectors in basis
        rot2_angle = vec_angle(alpha, beta)
        cross = np.cross(mid_direction, T_direction)
        if np.dot(direction, cross) < 0:
            rot2_angle *= -1
        
        in_centre = factory.getCenterOfMass(2, surfaces[3][1])
        self.mid_centre = factory.getCenterOfMass(2, surfaces[4][1])
        
        print(rot2_angle)
        factory.rotate([vol_tag], *in_centre, *list(direction), rot2_angle)
        factory.synchronize()
        
        super(T_junction, self).__init__(radius, vol_tag, in_surface, out_surface, out_direction, out_direction)
    
    def update_centres(self):
        """
        need to update T centre too
        """
        self.vol_centre = np.array(factory.getCenterOfMass(*self.vol_tag))
        self.in_centre = np.array(factory.getCenterOfMass(*self.in_surface))
        self.out_centre = np.array(factory.getCenterOfMass(*self.out_surface))
        self.mid_centre = np.array(factory.getCenterOfMass(*self.mid_surface))