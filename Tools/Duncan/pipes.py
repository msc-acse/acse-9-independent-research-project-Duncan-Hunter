"""
Create pipes and improve this docstring.

Bugs:
    - There will be some.

Todo
    - Continue update of Surface class
    - Update how we add pieces
    - Aim is to be able to build from junctions.



    - mesh size of pieces after fuse
        Update: Basic version working - lots of bugs still
        to fix.
        Overall pipe mesh size
        Apply to all networks joined.
    - Junctions
        Make side without junction smaller.
        Be able to build from junctions, therefore knowing what direction
        you're going in.
    - Docstrings
    - Physical groups
        Need to deal with excess outflows or inflows
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
        self.in_direction = np.array(direction)
        # self.in_surface = piece.in_surface  # unsure of use
        # self.in_centre = piece.in_centre
        # self.direction = piece.out_direction  # direction of pipe
        self.radius = radius

        # self.out_centre = piece.out_centre
        self.piece_list = [piece]
        self.in_surfaces = [piece.in_surface]  # need to be careful with this
        self.out_surfaces = [piece.out_surface]

        self.t_pieces = []  # for use when adding a pipe.
        self.networks = [self]  # for use when adding a pipe.
        self.joined = False  # for use when adding a pipe.

        self.fused = False
        self.vol_tag = None  # Overall vol_tag, only created after fuse.
        self._has_physical_groups = False
        self.physical_in_surface = None
        self.physical_out_surface = None
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

    def _update_pieces(self):
        """
        Add a piece to the pipe
        Specify the type of piece and where to add it.
        Use piece type to determine what piece. Give individual
        functions the out_surface to determine where the out centre
        and direction is.
        Replace the out_surface with the new out piece.
        """
        # Update all pieces
        # Update out pieces
        # Update in pieces
        for piece in self.piece_list:
            piece.update()

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
        piece.update()
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
        piece = pieces.Mitered(out_surface.radius, out_surface.out_direction,
                               new_direction, lcar)
        translate_vector = out_surface.out_centre - piece.in_centre
        FACTORY.translate([piece.vol_tag], *list(translate_vector))
        FACTORY.synchronize()
        piece.update_centres()
        self.piece_list.append(piece)
        self.out_surfaces[out_number] = piece

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
        if new_radius > out_surface.radius:
            piece = pieces.Cylinder(length, new_radius, out_surface.out_direction, lcar)
        else:
            piece = pieces.Cylinder(length, out_surface.radius, out_surface.out_direction, lcar)
        translate_vector = out_surface.out_centre - piece.in_centre
        FACTORY.translate([piece.vol_tag], *list(translate_vector))
        piece.update_centres()
        FACTORY.synchronize()
        if new_radius > out_surface.radius:
            lines = MODEL.getBoundary([piece.in_surface], False, False)
            FACTORY.chamfer([piece.vol_tag[1]], [lines[0][1]],
                            [piece.in_surface[1]],
                            [new_radius - out_surface.radius, change_length])
            FACTORY.synchronize()
        elif new_radius < out_surface.radius:
            lines = MODEL.getBoundary(piece.out_surface, False, False)
            FACTORY.chamfer([piece.vol_tag[1]], [lines[0][1]],
                            [piece.out_surface[1]],
                            [out_surface.radius - new_radius, change_length])
            FACTORY.synchronize()
        piece.radius = new_radius
        self.piece_list.append(piece)
        self.out_surfaces[out_number] = piece

    def add_t_junction(self, t_direction, lcar, out_number=0):
        """Adds a T junction to the Network at the outlet.

        This represents a pipe joining this pipe, creating a place to
        add a Network to this Network.

        Args:
            t_direction: (list, length 3) representing the direction
                that the joining pipe's inlet is facing.
        """
        out_number = self._out_number(out_number)
        out_surface = self.out_surfaces
        piece = pieces.TJunction(self.radius, self.direction, t_direction, lcar)
        translate_vector = self.out_centre - piece.in_centre
        FACTORY.translate([piece.vol_tag], *list(translate_vector))
        piece.update_centres()
        FACTORY.synchronize()
        self.t_pieces.append(piece)
        self.piece_list.append(piece)
        self._update_properties()

    def fuse_objects(self):
        #fuse_tool = self.entities[0]
        """Fuses separate objects in Network together.

        Finds out which surface is inflow/outflow using norms.
        Updates inflow_tag, outflow_tag. Creates the physical
        surfaces needed for ICFERST.
        """
        vol_tags = [piece.vol_tag for piece in self.piece_list]
        # out_dim_tags = FACTORY.fuse([self.entities[0]], self.entities[1:])[0]
        out_dim_tags = FACTORY.fuse([vol_tags[0]], vol_tags[1:])[0]
        FACTORY.synchronize()
        self.vol_tag = out_dim_tags[0]
        for piece in self.piece_list:
            piece.vol_tag = None
        surfaces = MODEL.getBoundary([self.vol_tag], False)
        n_surfaces = len(surfaces)
        planes = [
            surfaces[i] for i in range(n_surfaces)
            if MODEL.getType(*surfaces[i]) == "Plane"
        ]
        no_slip = [
            surfaces[i] for i in range(n_surfaces)
            if MODEL.getType(*surfaces[i]) != "Plane"
        ]
        no_slip_tags = [dimtag[1] for dimtag in no_slip]
        found_in = False  # being extra safe
        found_out = False
        # Find in surface. use np.isclose
        for dim, tag in planes:
            loc = np.array(FACTORY.getCenterOfMass(dim, tag))
            if np.allclose(loc, self.in_centre) and not found_in:
                self.in_surface = (dim, tag)
                found_in = True
            elif np.allclose(loc, self.out_centre) and not found_out:
                self.out_surface = (dim, tag)
                self.out_centre = loc
                found_out = True

        if self._has_physical_groups:
            MODEL.removePhysicalGroups()

        self.physical_in_surface = MODEL.addPhysicalGroup(
            2, [self.in_surface[1]])
        self.physical_out_surface = MODEL.addPhysicalGroup(
            2, [self.out_surface[1]])
        self.physical_no_slip = MODEL.addPhysicalGroup(2, no_slip_tags)
        self.physical_volume = MODEL.addPhysicalGroup(3, [self.vol_tag[1]])
        self._has_physical_groups = True
        # Could init a PipePiece here.
        piece = PipePiece(self.radius, self.vol_tag, self.in_surface,
                          self.out_surface, self.in_direction,
                          self.direction, 0.2)
        # Uncertain this is the best thing to do
        # self.piece_list = [piece]

    def rotate_network(self, old_direction, new_direction):
        """Rotates the network from old_direction to new_direction."""
        old_direction = np.array(old_direction)
        new_direction = np.array(new_direction)
        rot_angle = vec_angle(old_direction, new_direction)
        rot_axis = np.cross(old_direction, new_direction)
        dimtags = []
        for network in self.networks:
            for piece in network.piece_list:
                dimtags.append(piece.vol_tag)
        FACTORY.rotate(dimtags, 0, 0, 0, *list(rot_axis), rot_angle)
        FACTORY.synchronize()
        self._update_networks()

    def translate_network(self, vector):
        vector = np.array(vector)
        dimtags = []
        for network in self.networks:
            for piece in network.piece_list:
                dimtags.append(piece.vol_tag)
        FACTORY.translate(dimtags, *list(vector))
        FACTORY.synchronize()
        self._update_networks()

    def _update_networks(self):
        MESH.generate(3)
        for network in self.networks:
            for piece in network.piece_list:
                piece.update()  # updates locations
            network._update_properties()  # updates out centre

    def add_network(self, network, junction_number=-1):
        """Adds a premade network to junction junction_number

        If junction number is -1, then it adds to the first junction
        face. Junctions are ordered in the order they are added, for
        example if it is added first, its junction number is 1.
        Find if junction number exists (raise ValueError)
        rotate the new pipe so its outlet is facing the opposite to
        the direction of the t_junction
        translate the pipe

        add networks to list of networks.
        say that the joining network can't be added to (optional but will help)
        add t_junctions in joining to this network.
        """
        print("We're working on it")
        if junction_number > len(self.t_pieces):
            raise ValueError("""Junction number is higher than number
                             of junctions""")
        if junction_number == 0:
            raise ValueError("Can't use 0th junction")
        if junction_number == -1:
            junction_number = 1
        t_piece = self.t_pieces[junction_number - 1]
        t_direction = t_piece.t_direction
        network.rotate_network(network.direction,
                               -t_direction)
        network._update_networks()
        translate_vector =  t_piece.mid_centre - network.out_centre
        network.translate_network(translate_vector)
        network._update_networks()
        self.networks += network.networks
