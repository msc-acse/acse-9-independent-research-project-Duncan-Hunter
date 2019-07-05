
import pipes

import gmsh  # Download gmsh.py, and libgmsh files from gmsh-sdk
import os
import numpy as np

model = gmsh.model
factory = model.occ
mesh = model.mesh

gmsh.initialize()
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.2)
gmsh.model.add("Example")

"""
Uncomment for different example meshes
Change filename
"""
fname = "pipe"

# Pieces available
# piece = pipes.Cylinder(1, 0.5, [1,0,0], 0.2)
# piece = pipes.Mitered(0.5, [-1,0,1], [-1,0,0], 0.2)
# piece = pipes.Curve(0.5, [1,0,0], [0,1,0], 1, 0.2)
# piece = pipes.T_junction(0.5, [1,0,0], [1,1,-1], 0.2)

# Create pipe with junctions
# network = pipes.Network(0.2, 0.25, [0,0,-1], 0.1)
# network.add_pipe(1, 0.1)
# network.add_t_junction([-1,0,1], 0.1)
# network.add_t_junction([1,0,1], 0.1)

# network2 = pipes.Network(0.2, 0.25, [1,0,0], 0.1)
# network2.add_pipe(3, 0.1)

# network3 = pipes.Network(0.2, 0.25, [1,0,0], 0.1)
# network3.add_pipe(3, 0.1)

# network.add_network(network2, junction_number=1)
# network.add_network(network3, junction_number=2)
# network.fuse_objects()
# network._set_mesh_sizes()

# Recreate pipe from mpml file
# network.add_curve([0,-1,0], 0.5, 0.2)
# network.add_pipe(5, 0.2)
# network.add_curve([0,0,-1],0.5,0.2)
# network.add_pipe(10,0.2)
# network.fuse_objects()

# Chicane
# network.add_pipe(1, 0.1)
# network.add_curve([0, -1, 0], 1, 0.05)
# network.add_pipe(0.5, 0.15)
# network.add_curve([1, 1, 0], 2, 0.05)
# network.add_pipe(4, 0.1)
# network.fuse_objects()
# network._set_mesh_sizes()

# Junction
# network.add_pipe(1, 0.2)
# network.add_T_junction([1,0,1])
# network.fuse_objects()

# Spiral - using loops
#for i in range(3):
#   network.add_pipe(1,0.2)
#   network.add_mitered([0, 0, 1], 0.2)
#   network.add_pipe(0.1, 0.2)
#   network.add_mitered([0, 1, 0], 0.2)
#   network.add_pipe(1, 0.2)
#   network.add_curve([-1,0,0],[-0.55,0,0],0.2)
#   network.add_pipe(1, 0.2)
#   network.add_mitered([0,0,1],0.2)
#   network.add_pipe(0.1, 0.2)
#   network.add_mitered([0,-1,0],0.2)
#   network.add_pipe(1, 0.2)
#   network.add_curve([1, 0, 0], [0.55, 0, 0], 0.2)
#network.fuse_objects()

# U bend positive x to -x, through -y
#network.add_pipe(1, 0.2)
#network.add_curve([0, -1, 0], [0, -1, 0], 0.1)
#network.add_curve([-1, 0, 0], [-1, 0, 0], 0.1)
#network.add_pipe(1, 0.2)
#network.fuse_objects()

# U bend with changing radius
#network.change_radius(1, 0.4, 0.7)
#network.add_curve([0, -1, 0], [0, -1, 0], 0.1)
#network.add_curve([-1, 0, 0], [-1, 0, 0], 0.1)
#network.change_radius(1, 0.5, 0.2)
#network.fuse_objects()

# Sharp Chicane
# network.add_pipe(1, 0.2)
# network.add_mitered([0,1,0], 0.2)
# network.add_pipe(1, 0.2)
# network.add_mitered([1,0,0], 0.2)
# network.add_pipe(1, 0.2)
# network.fuse_objects()

"""
Some notes
"""

mesh.generate(3)

# gmsh.option.setNumber("Mesh.Binary", 1)
# gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints", 1)
# gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)

# gmsh.write(fname + ".msh2")
# os.rename(fname + ".msh2", fname + ".msh")
gmsh.option.setNumber("General.Axes", 2)
gmsh.option.setNumber("Mesh.SurfaceFaces", 1)
# gmsh.option.setNumber("Geometry.LabelType", 2)
# gmsh.option.setNumber("Geometry.SurfaceNumbers", 1)
#model.setVisibility(model.getEntities(3),0)
gmsh.fltk.run()
gmsh.finalize()