import pipes
import gmsh  # Download gmsh.py, and libgmsh files from gmsh-sdk
import os
import numpy as np

model = gmsh.model
factory = model.occ
mesh = model.mesh

gmsh.initialize()
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.2)
# # gmsh.option.setNumber("General.Terminal", 1)
model.add("Example")

"""
Uncomment for different example meshes
Change filename
"""
fname = "chicane"


network = pipes.Network(0.1, 0.5, [1,0,0], 0.2)

# Junction
# network.add_pipe(1, 0.2)
# network.add_T_junction([0,0,1])
# network.fuse_objects()


# Spiral
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

# Chicane
# network.add_pipe(1, 0.1)
# network.add_curve([0, -1, 0], [0, -1, 0], 0.1)
# network.add_pipe(0.5, 0.1)
# network.add_curve([1, 0, 0], [2, 0, 0], 0.1)
# network.add_pipe(4, 0.1)
# network.fuse_objects()

#print("In", network.physical_in_surface)
#print("Out", network.physical_out_surface)
#print("Wall", network.physical_no_slip)
#print("Volume", network.physical_volume)

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



mesh.generate(3)

# gmsh.write(fname + ".msh2")
# os.rename(fname + ".msh2", fname + ".msh")
gmsh.option.setNumber("Mesh.Binary", 1)

gmsh.write(fname)

gmsh.option.setNumber("General.Axes", 2)
gmsh.option.setNumber("Mesh.SurfaceFaces", 1)
gmsh.option.setNumber("Geometry.LabelType", 2)
gmsh.option.setNumber("Geometry.SurfaceNumbers", 1)
#model.setVisibility(model.getEntities(3),0)
gmsh.fltk.run()
gmsh.finalize()

"""
Some notes
"""
# Setting mesh size

# stag = factory.addPoint(0, 0, 0, 0.2)
# ftag = factory.addPoint(1, 0, 0, 0.2)
# line = factory.addLine(stag, ftag)
# factory.synchronize()

# field_tag = mesh.field.add("Distance")

# Distance away from point(s)
# mesh.field.setNumbers(field_tag, "NodesList", [stag])

# Distance away from line
#mesh.field.setNumber(field_tag, "NNodesByEdge", 100)
#mesh.field.setNumbers(field_tag, "EdgesList", [line])

# Threshold field
# thresh_field = mesh.field.add("Threshold")
# mesh.field.setNumber(thresh_field, "IField", 1)
# mesh.field.setNumber(thresh_field, "LcMin", 0.05)
# mesh.field.setNumber(thresh_field, "LcMax", 0.3)
# mesh.field.setNumber(thresh_field, "DistMin", 0.3)
# mesh.field.setNumber(thresh_field, "DistMax", 0.5)

# mesh.field.setAsBackgroundMesh(thresh_field)
#model.mesh.field.setNumbers(1, "EdgesList", line)
#mesh.field.setNumbers

# Using factory.remove gets rid of the effect