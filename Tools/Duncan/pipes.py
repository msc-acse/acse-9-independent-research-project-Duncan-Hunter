"""
Create pipes and improve this docstring.

Bugs:

Todo
    - mesh size of pieces after fuse
        Update: Basic version working - lots of bugs still
        to fix.
        Overall pipe mesh size
    - Junctions
        Make side without junction smaller.
    - Docstrings
    - Physical groups
        Need to deal with excess inflows (use out_surfaces?)
    - Testing
        Test case with ICFERST
    - Standardize variable names
"""
import gmsh
import numpy as np
from scipy.spatial.transform import Rotation
import pieces

MODEL = gmsh.model
FACTORY = MODEL.occ
MESH = MODEL.mesh
#gmsh.option.setNumber("General.Terminal", 1)

def round_0(values):
    values = np.array(values)
    values[np.abs(values) < 1e-8] = 0
    return values

def vec_angle(vec1, vec2):
    """Returns the angle between two numpy array vectors"""
    return np.arccos(
        (np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))))


def proj(vec1, vec2):
    """
    Returns the component of vec1 along vec2

    Args:
        vec1, vec2: (np.array shape 3) xyz vector.
    """
    return np.dot(vec1, vec2) / np.linalg.norm(vec2)


def check_intersect(objects, tools):
    """
    Check if entities in tools intersect entities in objects.
    objects: (list) 1st group of entities
    tools: (list) 2nd group of entities

    returns: (bool) True if there is an intersection.
                    False if there is no intersection.
    """
    intersect = FACTORY.intersect(objects,
                                  tools,
                                  removeObject=False,
                                  removeTool=False)
    if intersect:
        FACTORY.remove(intersect, True)
        return True
    else:
        return False

class Network():
    """Represents a pipe or network of pipes.

    Pipes are built from an inlet in a sequential, modular fashion.

    Attributes:
        physical_out_surface:
        physical_in_surface:
        phyiscal_no_slip:
        physical_volume:
    """

    def __init__(self, length, radius, direction, lcar):
        """
        Creates the inlet cylinder of a pipe.

        This is the beginning of a pipe, from which you can add more
        pieces. This piece has the inlet surface, and stores the
        geometrical information of the pipe, such as what direction
        the outlet is facing, where it is, and what the radius of the
        pipe is.

        Args:
        length: (float) length of piece.
        radius: (float) radius of pipe.
        direction: (list) Direction pipe will be facing in x, y, z
            vector format.
        lcar: (float) Mesh size of this piece.
        """
        piece = pieces.Cylinder(length, radius, np.array(direction), lcar)
        self.in_direction = np.array(direction)  # not sure I'll need this
        self.radius = radius  # not sure i'll need this

        # self.out_centre = piece.out_centre
        self.piece_list = [piece]
        self.in_surfaces = [piece.in_surface]  # need to be careful with this
        self.out_surfaces = [piece.out_surface]

        self.networks = [self]  # for use when adding a pipe.
        # self.joined = False  # for use when adding a pipe.

        self.fused = False
        self.vol_tag = None  # Overall vol_tag, only created after fuse.
        self.physical_in_out_surfaces = {}  # dictionary of physical dim tags
        # to surface objects, which has information
        self.physical_no_slip = None
        self.physical_volume = None

    def _set_mesh_sizes(self):
        """Sets the mesh size for all pieces."""
        field_list = []
        for network in self.networks:
            for piece in network.piece_list:
                half_length = np.abs(np.linalg.norm(piece.in_surface.centre -
                                                    piece.vol_centre))
                field_length = np.linalg.norm(np.array([half_length,
                                                        piece.in_surface.radius]))
                point = FACTORY.addPoint(*list(piece.vol_centre),
                                meshSize=piece.lcar)
                FACTORY.synchronize()
                dist_field = MESH.field.add("Distance")
                MESH.field.setNumbers(dist_field, "NodesList", [point])
                thresh_field = MESH.field.add("Threshold")
                MESH.field.setNumber(thresh_field, "IField", dist_field)
                MESH.field.setNumber(thresh_field, "LcMin", piece.lcar)
                MESH.field.setNumber(thresh_field, "LcMax", 0.3)
                MESH.field.setNumber(thresh_field, "DistMin", field_length)
                MESH.field.setNumber(thresh_field, "DistMax", 1.1*field_length)
                field_list.append(thresh_field)
        min_field = MESH.field.add("Min")
        MESH.field.setNumbers(min_field, "FieldsList", field_list)
        MESH.field.setAsBackgroundMesh(min_field)

    def _out_number(self, out_number):
        if out_number > len(self.out_surfaces):
            raise ValueError("Out piece does not exist")
        if out_number <= 0:
            out_number = 0
        else:  # Number to index
            out_number -= 1
        return out_number

    def add_cylinder(self, length, lcar, out_number=0):
        """Adds a pipe to the Network at the outlet.

        Args:
            length: (float) length of pipe.
            lcar: (float) mesh size of piece.
        """
        out_number = self._out_number(out_number)
        out_surface = self.out_surfaces[out_number]
        piece = pieces.Cylinder(length, out_surface.radius, out_surface.direction, lcar)
        translate_vector = out_surface.centre - piece.in_surface.centre
        FACTORY.translate([piece.vol_tag], *list(translate_vector))
        FACTORY.synchronize()
        piece.update_centres()
        self.piece_list.append(piece)
        self.out_surfaces[out_number] = piece.out_surface

    def add_curve(self, new_direction, bend_radius, lcar, out_number=0):
        """Adds a curve to the Network at the outlet.

        Args:
            new_direction: (list) Direction pipe will be facing
                in x, y, z vector format.
                e.g. [0, 1, 0] faces positive y.
            bend_radius: (float) Radius of the bend.
            lcar: (float) Size of mesh in this piece.

        Raises:
            ValueError: new_direction vector isn't right size.
                Bend radius isn't big enough.
        """
        # Check input
        out_number = self._out_number(out_number)
        if len(new_direction) > 3:
            raise ValueError(
                """Array given is too long (should be length 3,
                for x, y, z)""")
        if bend_radius < 1.1 * self.radius:
            raise ValueError(
                """Bend of radius is not large enough, place
                centre_of_rotation further away"""
            )
        out_surface = self.out_surfaces[out_number]
        piece = pieces.Curve(out_surface.radius, out_surface.direction,
                             new_direction, bend_radius, lcar)
        translate_vector = out_surface.centre - piece.in_surface.centre
        FACTORY.translate([piece.vol_tag], *list(translate_vector))
        FACTORY.synchronize()
        piece.update_centres()
        self.piece_list.append(piece)
        self.out_surfaces[out_number] = piece.out_surface

    def add_mitered(self, new_direction, lcar, out_number=0):
        """Adds a mitered bend to the Network at the outlet.

        A mitered bend is a sharp change in direction. The piece
        contains planes that may have an effect when creating
        the physical surfaces at outlets.

        Args:
            new_direction: (list, length 3) xyz vector representing
                the new direction of the pipe.
            lcar: (float) size of mesh of this piece.
        """
        out_number = self._out_number(out_number)
        out_surface = self.out_surfaces[out_number]
        piece = pieces.Mitered(out_surface.radius, out_surface.direction,
                               new_direction, lcar)
        translate_vector = out_surface.centre - piece.in_surface.centre
        FACTORY.translate([piece.vol_tag], *list(translate_vector))
        FACTORY.synchronize()
        piece.update_centres()
        self.piece_list.append(piece)
        self.out_surfaces[out_number] = piece.out_surface

    def add_change_radius(self, length, new_radius, change_length, lcar,
                          out_number=0):
        """Adds a piece that changes the radius of the Network.

        The piece is length long, and changes the Network radius to
        new_radius, over change_length, which controls how gentle the
        change is.

        Args:
            length: (float) Length of the piece.
            new_radius: (float) radius to change to.
            change_length: (float) Length that the change takes
                place over. Must be less than length and > 0.
            lcar: (float) mesh size for this piece.

        Raises:
            ValueErrors: change_length is not between length and 0.
        """
        out_number = self._out_number(out_number)
        out_surface = self.out_surfaces[out_number]
        if change_length >= length:
            raise ValueError('change_length must be less than length')
        if change_length < 0:
            raise ValueError('change_length must be greater than 0')
        piece = pieces.ChangeRadius(length, change_length,
                                    out_surface.radius, new_radius,
                                    out_surface.direction, lcar)
        translate_vector = out_surface.centre - piece.in_surface.centre
        FACTORY.translate([piece.vol_tag], *list(translate_vector))
        FACTORY.synchronize()
        piece.update_centres()
        self.piece_list.append(piece)
        self.out_surfaces[out_number] = piece.out_surface

    def add_t_junction(self, t_direction, lcar, out_number=0):
        """Adds a T junction to the Network at the outlet.

        This represents a pipe joining this pipe, creating a place to
        add a Network to this Network.

        Args:
            t_direction: (list, length 3) representing the direction
                that the joining pipe's inlet is facing.
        """
        out_number = self._out_number(out_number)
        out_surface = self.out_surfaces[out_number]
        piece = pieces.TJunction(out_surface.radius, out_surface.direction, t_direction, lcar)
        translate_vector = out_surface.centre - piece.in_surface.centre
        FACTORY.translate([piece.vol_tag], *list(translate_vector))
        FACTORY.synchronize()
        piece.update_centres()
        self.piece_list.append(piece)
        self.out_surfaces.append(piece.t_surface)
        self.out_surfaces[out_number] = piece.out_surface

    def _fuse_objects(self):
        #fuse_tool = self.entities[0]
        """Fuses separate objects in Network together.

        Returns the dimtags of in and out surfaces.
        """
        if self.vol_tag:
            raise ValueError("Network already fused")
        vol_tags = [piece.vol_tag for piece in self.piece_list]
        out_dim_tags = FACTORY.fuse([vol_tags[0]], vol_tags[1:])[0]
        FACTORY.synchronize()
        self.vol_tag = out_dim_tags[0]
        for piece in self.piece_list:
            piece.vol_tag = None
        surfaces = MODEL.getBoundary([self.vol_tag], False)
        n_surfaces = len(surfaces)

        # Dimtags have changed
        planes = [
            surfaces[i] for i in range(n_surfaces)
            if MODEL.getType(*surfaces[i]) == "Plane"
        ]
        no_slip = [
            surfaces[i] for i in range(n_surfaces)
            if MODEL.getType(*surfaces[i]) != "Plane"
        ]
        # Iterate through planes
        for dim, tag in planes:
            added = False
            loc = np.array(FACTORY.getCenterOfMass(dim, tag))
            # Iterate over in surfaces
            for surf in self.in_surfaces:
                # If it is an in surface, add it.
                if np.allclose(loc, surf.centre):
                    surf.dimtag = (dim, tag)
                    added = True  # skip extra iterations
            # If it is an out surface, add it.
            if not added:
                for surf in self.out_surfaces:
                    if np.allclose(loc, surf.centre):
                        surf.dimtag = (dim, tag)
        return no_slip

    def set_physical_groups(self):
        """Use fused to create physical groups"""
        no_slip = self._fuse_objects()
        no_slip_tags = [tag[1] for tag in no_slip]
        self.physical_no_slip = MODEL.addPhysicalGroup(2, no_slip_tags)
        for surf in (self.out_surfaces + self.in_surfaces):
            phys_tag = MODEL.addPhysicalGroup(2, [surf.dimtag[1]])
            self.physical_in_out_surfaces[phys_tag] = surf
        self.physical_volume = MODEL.addPhysicalGroup(3, [self.vol_tag[1]])

    def rotate_network(self, axis, angle):
        """Rotates the network from old_direction to new_direction."""
        rot_axis = np.array(axis)
        dimtags = []
        for network in self.networks:
            for piece in network.piece_list:
                dimtags.append(piece.vol_tag)
        FACTORY.rotate(dimtags, 0, 0, 0, *list(rot_axis), angle)
        FACTORY.synchronize()
        for piece in self.piece_list:
            piece.update_centres()
            piece.update_directions(rot_axis, angle)

    def translate_network(self, vector):
        vector = np.array(vector)
        dimtags = []
        for network in self.networks:
            for piece in network.piece_list:
                dimtags.append(piece.vol_tag)
        FACTORY.translate(dimtags, *list(vector))
        FACTORY.synchronize()
        for piece in self.piece_list:
            piece.update_centres()

    def write_info(self, fname):
        with open(fname, 'w') as myfile:
            myfile.writelines("Physical Surface Number, centre, direction")
            for key in self.physical_in_out_surfaces:
                surf = self.physical_in_out_surfaces[key]
                myfile.writelines("\n{}\t{}\t{}".format(key, round_0(surf.centre), round_0(surf.direction)))
            myfile.close()