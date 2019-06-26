import pipes
import gmsh  # Download gmsh.py, and libgmsh files from gmsh-sdk
import numpy as np

model = gmsh.model
factory = model.occ
mesh = model.mesh

gmsh.initialize()
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.2)
# gmsh.option.setNumber("General.Terminal", 1)
model.add("Example")

"""
Uncomment for different example meshes
Change filename
"""

fname = "test"

# Spiral
#for i in range(3):
#    network.add_pipe(1,0.2)
#    network.add_mitered([0, 0, 1])
#    network.add_mitered([0, 1, 0])
#    network.add_pipe(1, 0.2)
#    network.add_curve([-1,0,0],[-0.55,0,0],0.2)
#    network.add_pipe(1, 0.2)
#    network.add_mitered([0,0,1])
#    network.add_pipe(0.1, 0.2)
#    network.add_mitered([0,-1,0])
#    network.add_pipe(1, 0.2)
#    network.add_curve([1, 0, 0], [0.55, 0, 0], 0.2)
#network.fuse_objects()

#network.add_mitered([1,1,0])
#network.add_mitered([1,0,0])

#Chicane
#network.add_pipe(1, 0.1)
#network.add_curve([0, -1, 0], [0, -1, 0], 0.1)
#network.add_pipe(0.5, 0.1)
#network.add_curve([1, 0, 0], [2, 0, 0], 0.1)
#network.add_pipe(4, 0.1)
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
#network.add_pipe(1, 0.2)
#network.add_mitered([0,1,0])
#network.add_pipe(1, 0.2)
#network.add_mitered([1,0,0])
#network.add_pipe(1, 0.2)
#network.fuse_objects()

mesh.generate(3)

gmsh.option.setNumber("Mesh.Binary", 1)
gmsh.write(fname + ".msh2")
os.rename(fname + ".msh2", fname + ".msh")

gmsh.option.setNumber("General.Axes", 2)
gmsh.option.setNumber("Mesh.SurfaceFaces", 1)
gmsh.option.setNumber("Geometry.LabelType", 2)
gmsh.option.setNumber("Geometry.SurfaceNumbers", 1)
model.setVisibility(model.getEntities(3),0)
gmsh.fltk.run()
gmsh.finalize()