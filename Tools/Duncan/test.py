
import gmsh  # Download gmsh.py, and libgmsh files from gmsh-sdk
import os
import numpy as np
import pieces
import pipes

model = gmsh.model
factory = model.occ
mesh = model.mesh

gmsh.initialize()
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.1)
gmsh.model.add("Example")

"""
Uncomment for different example meshes
Change filename
"""
fname = "junc"


# Piece classes
# piece = pieces.Cylinder(1, 0.5, [1,0,0], 0.1)
# piece = pieces.ChangeRadius(2,1.8, 0.3, 0.2, [1,0,0], 0.1)
# piece = pieces.Mitered(0.5, [1,1,0], [1,0,0], 0.2)
# piece = pieces.Curve(0.5, [1,0,-1], [0,1,0], 1, 0.2)
piece = pieces.TJunction(0.5, [1,0,0], [-2,1,0], 0.1)

# Start a network
# network = pipes.Network(1, 0.3, [1,0,0], 0.1)

# Create pipe with junctions
# network.add_t_junction([-1,1,0], 0.1)
# network.add_t_junction([-1,-1,0], 0.1)
# network.add_cylinder(1, 0.1, out_number=2)
# network.add_curve([-1,0,0], 0.5, 0.1, out_number=3)
# network.add_cylinder(1.5, 0.1, out_number=3)

# Chicane
# network.add_curve([0,-1,0], 0.5, 0.2)
# network.add_pipe(5, 0.2)
# network.add_curve([0,0,-1],0.5,0.2)
# network.add_pipe(10,0.2)

# U bend positive x to -x, through -y
#network.add_pipe(1, 0.2)
#network.add_curve([0, -1, 0], [0, -1, 0], 0.1)
#network.add_curve([-1, 0, 0], [-1, 0, 0], 0.1)
#network.add_pipe(1, 0.2)

# U bend with changing radius
#network.change_radius(1, 0.4, 0.7)
#network.add_curve([0, -1, 0], [0, -1, 0], 0.1)
#network.add_curve([-1, 0, 0], [-1, 0, 0], 0.1)
#network.change_radius(1, 0.5, 0.2)

# Sharp Chicane
# network.add_cylinder(1, 0.1)
# network.add_mitered([0,1,0], 0.1)
# network.add_cylinder(1, 0.1)
# network.add_mitered([1,0,0], 0.1)
# network.add_cylinder(1, 0.1)

# Utility functions
# network.set_physical_groups()
# network._set_mesh_sizes()
# network.write_info("info.csv")


mesh.generate(3)
# os.rename(fname + ".msh2", fname + ".msh")
gmsh.option.setNumber("General.Axes", 2)
gmsh.option.setNumber("Mesh.SurfaceFaces", 1)

gmsh.fltk.run()

gmsh.finalize()