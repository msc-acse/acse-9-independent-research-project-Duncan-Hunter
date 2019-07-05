"""
Create pipes and improve this docstring.

Bugs:
    - There will be some.

Todo
    - mesh size of pieces after fuse
        Update: Basic version working - lots of bugs still
        to fix.
        Overall pipe mesh size
        Apply to all networks joined.
    - Junctions
        Make side without junction smaller.
        Update in direction of joining network
        Reliable way of knowing in_direction of joined network?
        After junction has been joined, remove junction from list
        Make joined pipe not able to add to
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
        piece = Cylinder(length, radius, np.array(direction), lcar)
        self.in_direction = np.array(direction)
        self.in_surface = piece.in_surface  # unsure of use
        self.in_centre = piece.in_centre
        self.direction = piece.out_direction  # direction of pipe
        self.radius = radius

        self.out_centre = piece.out_centre
        self.piece_list = [piece]

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
                half_length = np.abs(np.linalg.norm(piece.in_centre -
                                                    piece.vol_centre))
                field_length = np.linalg.norm(np.array([half_length,
                                                        piece.radius]))
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

    def _update_properties(self):
        """Updates properties of pipe after a piece is added."""
        self.out_centre = self.piece_list[-1].out_centre
        self.direction = self.piece_list[-1].out_direction
        self.in_centre = self.piece_list[0].in_centre

    def add_pipe(self, length, lcar):
        """Adds a pipe to the Network at the outlet.

        Args:
            length: (float) length of pipe.
            lcar: (float) mesh size of piece.
        """
        piece = Cylinder(length, self.radius, self.direction, lcar)

        translate_vector = self.out_centre - piece.in_centre
        FACTORY.translate([piece.vol_tag], *list(translate_vector))
        FACTORY.synchronize()
        piece.update_centres()
        self.piece_list.append(piece)
        self._update_properties()

    def add_curve(self, new_direction, bend_radius, lcar):
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
        if len(new_direction) > 3:
            raise ValueError(
                "Array given is too long (should be length 3, for x, y, z)")
        if bend_radius < 1.1 * self.radius:
            raise ValueError(
                "Bend of radius is not large enough, place centre_of_rotation further away"
            )

        piece = Curve(self.radius, self.direction, new_direction, bend_radius,
                      lcar)
        translate_vector = self.out_centre - piece.in_centre
        FACTORY.translate([piece.vol_tag], *list(translate_vector))
        FACTORY.synchronize()
        piece.update_centres()
        self.piece_list.append(piece)
        self._update_properties()

    def add_mitered(self, new_direction, lcar):
        """Adds a mitered bend to the Network at the outlet.

        A mitered bend is a sharp change in direction. The piece
        contains planes that may have an effect when creating
        the physical surfaces at outlets.

        Args:
            new_direction: (list, length 3) xyz vector representing
                the new direction of the pipe.
            lcar: (float) size of mesh of this piece.
        """
        piece = Mitered(self.radius, self.direction, new_direction, lcar)

        translate_vector = self.out_centre - piece.in_centre
        FACTORY.translate([piece.vol_tag], *list(translate_vector))
        FACTORY.synchronize()
        piece.update_centres()
        self.piece_list.append(piece)
        self._update_properties()

    def add_change_radius(self, length, new_radius, change_length, lcar):
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
        if change_length >= length:
            raise ValueError('change_length must be less than length')
        if change_length < 0:
            raise ValueError('change_length must be greater than 0')
        if new_radius > self.radius:
            piece = Cylinder(length, new_radius, self.direction, lcar)
        else:
            piece = Cylinder(length, self.radius, self.direction, lcar)
        translate_vector = self.out_centre - piece.in_centre
        FACTORY.translate([piece.vol_tag], *list(translate_vector))
        piece.update_centres()
        FACTORY.synchronize()
        if new_radius > self.radius:
            lines = MODEL.getBoundary([piece.in_surface], False, False)
            FACTORY.chamfer([piece.vol_tag[1]], [lines[0][1]],
                            [piece.in_surface[1]],
                            [new_radius - self.radius, change_length])
            FACTORY.synchronize()
        elif new_radius < self.radius:
            lines = MODEL.getBoundary(piece.out_surface, False, False)
            FACTORY.chamfer([piece.vol_tag[1]], [lines[0][1]],
                            [piece.out_surface[1]],
                            [self.radius - new_radius, change_length])
            FACTORY.synchronize()
        self.radius = new_radius
        self.piece_list.append(piece)
        self._update_properties()

    def add_t_junction(self, t_direction, lcar):
        """Adds a T junction to the Network at the outlet.

        This represents a pipe joining this pipe, creating a place to
        add a Network to this Network.

        Args:
            t_direction: (list, length 3) representing the direction
                that the joining pipe's inlet is facing.
        """
        piece = TJunction(self.radius, self.direction, t_direction, lcar)
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
        dimtags = []
        for network in self.networks:
            for piece in network.piece_list:
                dimtags.append(piece.vol_tag)
        FACTORY.translate(dimtags, *list(vector))
        FACTORY.synchronize()
        self._update_networks()

    def _update_networks(self):
        for network in self.networks:
            for piece in network.piece_list:
                piece.update_centres()  # updates locations
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

def _rotate_inlet(vol_tag, in_direction, out_direction):
    """Rotates up facing inlet to face in_direction.

    Calculates the new outlet direction after it has been
    transformed.

    Args:
        vol_tag: dimtag tuple of volume being rotated.
        in_direction: xyz array direction to rotate the inlet to.
        out_direction: Direction the outlet is facing before
            rotation.

    Returns:
        new_out_direction: Direction outlet is facing after
            rotation as xyz array.
    """
    up_vector = np.array([0, 0, 1])
    if np.allclose(in_direction, up_vector
                    ) is False:  # only have to rotate if its not facing up
        rotate_axis = np.cross(up_vector, in_direction)
        if np.allclose(rotate_axis, np.zeros(3)):
            rotate_axis = np.array([1, 0, 0])
        rotate_axis = rotate_axis / np.linalg.norm(rotate_axis)
        rotate_angle = np.arccos((
            np.dot(in_direction, up_vector) /  # can replace
            (np.linalg.norm(in_direction) * np.linalg.norm(up_vector))))

        FACTORY.rotate([vol_tag], *[0, 0, 0], *list(rotate_axis),
                        rotate_angle)
        FACTORY.synchronize()
        rot_vec = rotate_angle * rotate_axis
        rot1 = Rotation.from_rotvec(rot_vec)
        new_out_direction = rot1.apply(out_direction)
    return new_out_direction

def _rotate_outlet(vol_tag, out_direction, in_direction,
                       new_out_direction):
        """Rotates outlet about in_direction to face out_direction.

        Projects new_out_direction and out_direction onto basis axes
        that are perpendicular to each other and in_direction. The
        angle between the projections is found, and the object is
        rotated.

        Args:
            vol_tag: dimtag tuple of volume being rotated.
            out_direction: xyz array, the direction that the outlet
                will face after being rotated.
            in_direction: xyz array, the direction that the inlet is
                facing, and the axis that the object will be rotated
                about.
            new_out_direction: xyz array, the direction that the outlet
                faces before being rotated.
                Returned from _rotate_inlet.

        """
        basis_1 = np.cross(
            out_direction, in_direction
        )  # basis vectors perpendicular to direction (rotation axis)
        basis_2 = np.cross(basis_1,
                           in_direction)  # and perpendicular to other basis
        # Before rotation projection.
        alpha = np.array([
            proj(new_out_direction, basis_1),
            proj(new_out_direction, basis_2)
        ])
        # After rotation projection.
        beta = np.array(
            [proj(out_direction, basis_1),
             proj(out_direction, basis_2)])
        # Find angle between two vectors in bases.
        rot2_angle = vec_angle(alpha, beta)
        cross = np.cross(new_out_direction, out_direction)
        if np.dot(in_direction, cross) > 0:
            rot2_angle *= -1
        FACTORY.rotate([vol_tag], *[0, 0, 0], *list(in_direction), -rot2_angle)
        FACTORY.synchronize()

class PipePiece():
    """Parent class of pieces.

    Pieces are GMSH objects that can be used in creating pipes.
    This class has common information that all pieces have, such as
    radius. It also has functions that all the classes use, such as the
    need to update centres of pieces after they have been transformed.

    """

    def __init__(self, radius, vol_tag, in_tag, out_tag, in_direction,
                 out_direction, lcar):
        """Stores the information of a created piece.

        Args:
            radius: (float) radius of the piece.
            vol_tag: (tuple, length 2) GMSH dimtag (dimension, tag)
                representing volume.
            in_tag: (tuple, length 2) GMSH dimtag representing in surface
            out_tag: (tuple, length 2) GMSH dimtag representing out surface
            in_direction: (np array, shape 3) xyz vector representing
                direction going in.
            out_direction: (np array, shape 3) xyz vector representing
                direction going out.
        """
        self.radius = radius
        self.lcar = lcar
        self.vol_tag = vol_tag
        self.in_surface = in_tag
        self.out_surface = out_tag
        self.vol_centre = np.array(FACTORY.getCenterOfMass(*vol_tag))
        self.in_centre = np.array(FACTORY.getCenterOfMass(*in_tag))
        self.out_centre = FACTORY.getCenterOfMass(*out_tag)
        self.in_direction = in_direction
        self.out_direction = out_direction

    def get_surfaces(self):
        """This isnt used, but could be for the whole thing"""
        surfaces = MODEL.getBoundary([self.vol_tag], combined=False)
        self.in_surface = surfaces[2]
        self.out_surface = surfaces[1]

    def update_centres(self):
        """Updates centres of faces attributes of piece."""
        self.get_surfaces()  # overriden in most cases
        self.vol_centre = np.array(FACTORY.getCenterOfMass(*self.vol_tag))
        self.in_centre = np.array(FACTORY.getCenterOfMass(*self.in_surface))
        self.out_centre = np.array(FACTORY.getCenterOfMass(*self.out_surface))

class Cylinder(PipePiece):
    """
    Class representing a GMSH cylinder with base at 0,0,0 facing upwards.
    """

    def __init__(self, length, radius, direction, lcar):
        """Creates the cylinder with GMSH.

        Args:
            length: (float) length of cylinder.
            radius: (float) radius of cylinder.
            direction: (list, length 3) xyz vector representing direction
                cylinder is facing.
            lcar: (float) mesh size for this piece.
        """
        self.length = length
        self.lcar = lcar
        vol_tag = (3, FACTORY.addCylinder(0, 0, 0, 0, 0, length, radius))
        FACTORY.synchronize()
        surfaces = MODEL.getBoundary([vol_tag], False)
        in_surface = surfaces[2]
        out_surface = surfaces[1]
        direction = np.array(direction)
        # up_vector = np.array([0, 0, 1])

        _rotate_inlet(vol_tag, direction, direction)

        super(Cylinder,
              self).__init__(radius, vol_tag, in_surface, out_surface,
                             direction, direction, lcar)

    def get_surfaces(self):
        surfaces = MODEL.getBoundary([self.vol_tag], combined=False)
        self.in_surface = surfaces[2]
        self.out_surface = surfaces[1]

class Curve(PipePiece):
    """Class representing a GMSH curve by revolution."""

    def __init__(self, radius, in_direction, out_direction, bend_radius, lcar):
        """Inits the piece.

        Initially with inlet facing down, and outlet facing in x-plane.

        Args:
            radius: (float) radius of the pipe
            in_direction: (list, length 3) xyz vector representing
                direction going in.
            out_direction: (list, length 3) xyz vector representing
                direction going out.
            bend_radius: (float) radius of the bend of the curve.
            lcar: (float) mesh size for this piece.
        """
        in_tag = (2, FACTORY.addDisk(0, 0, 0, radius, radius))
        in_direction = np.array(in_direction)
        out_direction = np.array(out_direction)

        revolve_axis = [0, 1, 0]
        centre_of_rotation = [bend_radius, 0, 0]
        angle = vec_angle(in_direction, out_direction)
        # Revolve in x plane, bend with radius bend_radius
        revolve_tags = FACTORY.revolve([in_tag], *centre_of_rotation,
                                       *revolve_axis, angle)
        FACTORY.synchronize()

        vol_tag = revolve_tags[1]

        new_out_direction = np.array(
            [np.sin(np.pi - angle), 0,
             -np.cos(np.pi - angle)])  # direction out is currently facing
        # Rotate so in_direction faces right way "Rot1"
        new_out_direction = _rotate_inlet(vol_tag, in_direction,
                                               new_out_direction)
        # Rotate so out_direction faces right way "Rot2"
        _rotate_outlet(vol_tag, out_direction, in_direction,
                       new_out_direction)

        surfaces = MODEL.getBoundary([vol_tag], False, True)
        out_tag = surfaces[2]
        in_tag = surfaces[1]
        super(Curve, self).__init__(radius, vol_tag, in_tag, out_tag,
                                    in_direction, out_direction, lcar)

        # FACTORY.synchronize()
    def get_surfaces(self):
        surfaces = MODEL.getBoundary([self.vol_tag], combined=False)
        self.in_surface = surfaces[1]
        self.out_surface = surfaces[2]

class Mitered(PipePiece):
    """Class representing a mitered (sharp) pipe bend.

    Piece creation is done by masking (intersect) a cylinder with a
    chamfered box. The piece is then mirrored, rotated and fused.

    The piece is then rotated to face the direction of the outflow.
    It is then rotated about the direction of outflow to match the new direction
    """

    def __init__(self, radius, in_direction, out_direction, lcar):
        """Creates the GMSH piece.

        The inlet is facing up originally.

        Args:
            radius: (float) radius of the pipe.
            in_direction: (list, length 3) xyz vector representing
            direction going in.
            out_direction: (list, length 3) xyz vector representing
            direction going out.
            lcar: (float) mesh size for this piece.

        """
        in_direction = np.array(in_direction)
        out_direction = np.array(out_direction)  # clean up v's

        # Chamfer cylinder
        angle = vec_angle(out_direction, in_direction)
        height = 2.1 * radius * np.tan(angle / 2)
        new_out_direction = np.array([np.sin(np.pi - angle), 0, -np.cos(np.pi - angle)
                       ])  # original outlet direction in xz plane
        cyl1 = (3, FACTORY.addCylinder(0, 0, 0, 0, 0, height,
                                       radius))  # create cylinder
        box1 = (3,
                FACTORY.addBox(-radius - 1, -radius, -1, 2 * radius + 1,
                               2 * radius, height + 1))  # create box
        FACTORY.synchronize()
        surface = MODEL.getBoundary([box1], False,
                                    False)[5]  # Top surface to chamfer from
        line = MODEL.getBoundary([surface], False, False)[2]  # -x line on top
        sdist = 2 * radius * np.tan(
            angle / 2)  # distances to chamfer to get correct angle
        #ldist = 2*radius*np.cos(angle/2)
        FACTORY.chamfer([box1[1]], [line[1]], [surface[1]],
                        [2 * radius, sdist])  # chamfer
        int_tag = FACTORY.intersect([box1],
                                    [cyl1])  # intersect (chamfer on cylinder)
        fuse = FACTORY.fuse([int_tag[0]], int_tag[1:])[0]  # fuse
        fuse2 = FACTORY.copy([fuse])  # create copy and mirror
        FACTORY.symmetrize([fuse2], 1, 0, 0, 0)
        FACTORY.synchronize()
        surface = MODEL.getBoundary([fuse2], False,
                                    False)[1]  # get center of rotation
        com = FACTORY.getCenterOfMass(*surface)
        FACTORY.rotate([fuse2], *com, 0, 1, 0,
                       -(np.pi - angle))  # rotate to create piece
        vol_tag = FACTORY.fuse([fuse], [fuse2],
                               removeObject=True,
                               removeTool=True)[0][0]  # fuse
        FACTORY.synchronize()

        surfaces = MODEL.getBoundary([vol_tag], False)
        in_surface = surfaces[3]  # bottom
        out_surface = surfaces[6]  # angled

        # Rot1: rotate object so inlet faces correct direction
        new_out_direction = _rotate_inlet(vol_tag, in_direction, new_out_direction)
        #Rot2
        _rotate_outlet(vol_tag, out_direction, in_direction, new_out_direction)
        super(Mitered, self).__init__(radius, vol_tag, in_surface, out_surface,
                                      in_direction, out_direction, lcar)

    def get_surfaces(self):
        surfaces = MODEL.getBoundary([self.vol_tag], combined=False)
        self.in_surface = surfaces[3]
        self.out_surface = surfaces[6]

class TJunction(PipePiece):
    """Class representing a T-junction in GMSH"""

    def __init__(self, radius, direction, t_direction, lcar, T_type='in'):
        """Inits the piece

        Creates a piece with the t_direction facing in the x-plane and
        inlet facing up.

        Args:
            radius: (float) radius of the pipe.
            direction: (list, length 3) xyz vector representing
                direction going in.
            t_direction: (list, length 3) xyz vector represeting the
                direction that the junction inlet faces.
        """
        self.out_num = 2
        direction = np.array(direction)
        self.t_direction = np.array(t_direction)
        t_angle = vec_angle(direction, self.t_direction)
        if t_angle < np.pi / 2:
            rot_sign = -1
            beta = np.pi / 2 - abs(t_angle)
        else:
            rot_sign = 1
            beta = abs(t_angle) - np.pi / 2

        height = radius * np.tan(beta) + radius / np.cos(beta)
        # Calculating height needed to emerge from merge

        in_tag = (3, FACTORY.addCylinder(0, 0, 0, 0, 0, 1.1 * height, radius))
        mid_tag = (3, FACTORY.addCylinder(0, 0, 0, 1.1 * height, 0, 0, radius))
        out_tag = (3, FACTORY.addCylinder(0, 0, 0, 0, 0, -1.1 * height,
                                          radius))
        FACTORY.rotate([mid_tag], 0, 0, 0, 0, 1, 0, rot_sign * beta)
        FACTORY.synchronize()

        vol_tags, out_dim_tag_map = FACTORY.fuse([in_tag], [mid_tag, out_tag])
        vol_tag = vol_tags[0]
        FACTORY.synchronize()

        surfaces = MODEL.getBoundary([vol_tag], False)
        in_surface = surfaces[5]  # 5 3
        out_surface = surfaces[3]
        self.mid_surface = surfaces[4]

        mid_direction = np.array([1, 0, 0])
        mid_direction = _rotate_inlet(vol_tag, direction, mid_direction)
        out_direction = np.copy(direction)

        _rotate_outlet(vol_tag, t_direction, direction, mid_direction)

        self.mid_centre = FACTORY.getCenterOfMass(2, surfaces[4][1])

        super(TJunction,
              self).__init__(radius, vol_tag, in_surface, out_surface,
                             out_direction, out_direction, lcar)

    def update_centres(self):
        """See base class. Updates T_centre too."""
        self.get_surfaces()
        self.vol_centre = np.array(FACTORY.getCenterOfMass(*self.vol_tag))
        self.in_centre = np.array(FACTORY.getCenterOfMass(*self.in_surface))
        self.out_centre = np.array(FACTORY.getCenterOfMass(*self.out_surface))
        self.mid_centre = np.array(FACTORY.getCenterOfMass(*self.mid_surface))

    def get_surfaces(self):
        surfaces = MODEL.getBoundary([self.vol_tag], combined=False)
        self.in_surface = surfaces[5]
        self.out_surface = surfaces[3]
        self.mid_surface = surfaces[4]