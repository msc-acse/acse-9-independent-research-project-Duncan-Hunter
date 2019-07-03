"""
Create pipes and improve this docstring.
"""
import gmsh
import numpy as np
from scipy.spatial.transform import Rotation

MODEL = gmsh.model
FACTORY = MODEL.occ
MESH = MODEL.mesh
#gmsh.option.setNumber("General.Terminal", 1)
"""
Overall BUGS
    Todo
        - mesh size of pieces after fuse
            Centre to mesh size mapping
            Fields with locations
        - Junctions
            Ability to add to a junction
        - Keep track of surface tags in general
            Keep list of centres
        - Keep track of orientations ??
            Update - made easier by classes
        - Keep track of locations
            Update - made easier by classes
            append to vol_centres
            Use dictionary to map tags to centres?
            For post processing
        - Docstrings
        - Physical groups
            Need to deal with excess outflows or inflows
        - Testing
        - Standardize variable names
"""


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

        self.lcar = lcar
        self.length = length
        self.direction = np.array(direction)  # direction of pipe
        self.up_vector = np.array([0, 0, 1])
        self.radius = radius
        #self.lcar = lcar

        piece = Cylinder(length, radius, self.direction, lcar)

        self.in_surface = piece.in_surface  # unsure of use
        self.in_centre = piece.in_centre
        self.out_centre = piece.out_centre
        self.vol_tag = None

        # self._set_MESH_size(self.in_surface, lcar)

        self.entities = [piece.vol_tag]  # list of dimtags of entities in pipe
        self.lcar_list = [lcar]
        self.vol_centres = [piece.vol_centre]

        self._has_physical_groups = False
        self.physical_in_surface = None
        self.physical_out_surface = None
        self.physical_no_slip = None
        self.physical_volume = None

    def _set_mesh_sizes(self):
        """Sets the mesh size for all pieces."""
        # for point, lcar in zip(self.vol_centres, self.lcar_list):
        # self._set_point_size(point, lcar)

    def _set_point_size(self, point, field_length, lcar):
        """
        Sets the mesh size around a specific point.

        Args:
            point:
            field_length:
            lcar:
        """
        print(point, lcar)

    def _get_surfaces(self, dimtag):
        """Returns the dimtag of surfaces of a volume.

        NOT USED.
        dimtag: (tuple) dimtag of object you're trying to get the boundary of.
        Get in/out surfaces of object
        """
        surfaces = MODEL.getBoundary(dimtag, False, False, False)
        return surfaces
        #if model.getType(dim, tag) != 'Plane':

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

        # self._set_MESH_size(piece.vol_tag, lcar)
        self.out_centre = piece.out_centre
        self.entities.append(piece.vol_tag)
        self.lcar_list.append(lcar)
        self.vol_centres.append(piece.vol_centre)

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
        #piece.out_centre += translate_vector
        FACTORY.synchronize()
        piece.update_centres()

        #self._set_MESH_size(revolve_tags[1], lcar)
        self.direction = piece.out_direction
        self.out_centre = piece.out_centre
        self.entities.append(piece.vol_tag)  #(revolve_tags[1]))

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

        piece.update_centres()

        FACTORY.synchronize()

        self.out_centre = piece.out_centre
        self.direction = piece.out_direction
        self.entities.append(piece.vol_tag)

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
        self.out_centre = piece.out_centre
        self.radius = new_radius
        self.entities.append(piece.vol_tag)

    def add_t_junction(self, t_direction):
        """Adds a T junction to the Network at the outlet.

        This represents a pipe joining this pipe, creating a place to
        add a Network to this Network.

        Args:
            t_direction: (list, length 3) representing the direction
                that the joining pipe's inlet is facing.
        """
        piece = TJunction(self.radius, self.direction, t_direction)
        translate_vector = self.out_centre - piece.in_centre
        FACTORY.translate([piece.vol_tag], *list(translate_vector))
        piece.update_centres()
        FACTORY.synchronize()
        self.out_centre = piece.out_centre
        self.direction = piece.out_direction
        self.entities.append(piece.vol_tag)

    def fuse_objects(self):
        #fuse_tool = self.entities[0]
        """Fuses separate objects in Network together.

        Finds out which surface is inflow/outflow using norms.
        Updates inflow_tag, outflow_tag. Creates the physical
        surfaces needed for ICFERST.
        """
        out_dim_tags = FACTORY.fuse([self.entities[0]], self.entities[1:])[0]
        FACTORY.synchronize()
        self.vol_tag = out_dim_tags[0]
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
        self.entities = out_dim_tags
        # Find out surface

    def rotate_network(self, new_direction, out_flow_direction=True):
        """
        rotate pipe so outflow is pointed in new_direction
        update position of where in_centre and out_centre is
        """


class PipePiece():
    """Parent class of pieces.

    Pieces are GMSH objects that can be used in creating pipes.
    This class has common information that all pieces have, such as
    radius. It also has functions that all the classes use, such as the
    need to update centres of pieces after they have been transformed.

    """

    def __init__(self, radius, vol_tag, in_tag, out_tag, in_direction,
                 out_direction):
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
        self.vol_tag = vol_tag
        self.in_surface = in_tag
        self.out_surface = out_tag
        self.vol_centre = np.array(FACTORY.getCenterOfMass(*vol_tag))
        self.in_centre = np.array(FACTORY.getCenterOfMass(*in_tag))
        self.out_centre = FACTORY.getCenterOfMass(*out_tag)
        self.in_direction = in_direction
        self.out_direction = out_direction

    def update_centres(self):
        """Updates centres of faces attributes of piece."""
        self.vol_centre = np.array(FACTORY.getCenterOfMass(*self.vol_tag))
        self.in_centre = np.array(FACTORY.getCenterOfMass(*self.in_surface))
        self.out_centre = np.array(FACTORY.getCenterOfMass(*self.out_surface))

    def _rotate_inlet(self, vol_tag, in_direction, out_direction):
        """
        Rotates up facing inlet to face in_direction.

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

    def _rotate_outlet(self, vol_tag, out_direction, in_direction,
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
        in_surface = surfaces[2]  # these are switched
        out_surface = surfaces[1]
        direction = np.array(direction)
        up_vector = np.array([0, 0, 1])

        if np.allclose(direction, up_vector) is False:
            cross = np.cross(up_vector, direction)
            rotation_axis = np.cross(direction, up_vector)
            if np.allclose(rotation_axis,
                           np.array([0, 0, 0])):  # if rotating 180 degrees
                rotation_axis = np.array([1, 0, 0])
            rotation_angle = np.arccos((
                np.dot(direction, up_vector) /  # can replace
                (np.linalg.norm(direction) * np.linalg.norm(up_vector))))

            if np.dot(rotation_axis, cross) > 0:
                rotation_angle *= -1
            FACTORY.rotate([vol_tag], 0, 0, length / 2, rotation_axis[0],
                           rotation_axis[1], rotation_axis[2], -rotation_angle)
            FACTORY.synchronize()
            out_direction = np.copy(direction)
            in_direction = np.copy(direction)
        else:
            out_direction = np.copy(up_vector)
            in_direction = np.copy(up_vector)

        super(Cylinder,
              self).__init__(radius, vol_tag, in_surface, out_surface,
                             in_direction, out_direction)


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
        new_out_direction = self._rotate_inlet(vol_tag, in_direction,
                                               new_out_direction)
        # Rotate so out_direction faces right way "Rot2"
        self._rotate_outlet(vol_tag, out_direction, in_direction,
                            new_out_direction)

        surfaces = MODEL.getBoundary([vol_tag], False, True)
        out_tag = surfaces[2]
        in_tag = surfaces[1]
        super(Curve, self).__init__(radius, vol_tag, in_tag, out_tag,
                                    in_direction, out_direction)

        # FACTORY.synchronize()


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
        new_out_direction = self._rotate_inlet(vol_tag, in_direction, new_out_direction)
        #Rot2
        self._rotate_outlet(vol_tag, out_direction, in_direction, new_out_direction)
        super(Mitered, self).__init__(radius, vol_tag, in_surface, out_surface,
                                      in_direction, out_direction)


class TJunction(PipePiece):
    """Class representing a T-junction in GMSH"""

    def __init__(self, radius, direction, t_direction, T_type='in'):
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

        mid_direction = self._rotate_inlet(vol_tag, direction, mid_direction)
        out_direction = np.copy(direction)

        self._rotate_outlet(vol_tag, t_direction, direction, mid_direction)

        self.mid_centre = FACTORY.getCenterOfMass(2, surfaces[4][1])

        super(TJunction,
              self).__init__(radius, vol_tag, in_surface, out_surface,
                             out_direction, out_direction)

    def update_centres(self):
        """See base class. Updates T_centre too."""
        self.vol_centre = np.array(FACTORY.getCenterOfMass(*self.vol_tag))
        self.in_centre = np.array(FACTORY.getCenterOfMass(*self.in_surface))
        self.out_centre = np.array(FACTORY.getCenterOfMass(*self.out_surface))
        self.mid_centre = np.array(FACTORY.getCenterOfMass(*self.mid_surface))
