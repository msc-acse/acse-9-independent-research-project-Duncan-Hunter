# Tools that use the GMSH Python API

### Create classes that create GMSH objects and store information. 
e.g. cylinders, boxes. 
More complex objects such as extrusions, revolutions.

Classes store more information than GMSH provides.

*Easy to access information:*
- Location
  -   Centre of faces
  -   Centre of object
- Direction
  - Direction of faces
- GMSH tags
- GMSH dimensions	
- Radius*
- Length*

\* If applicable to that shape

*Easy to use functions:*
  - GMSH functions
    - Rotate, translate etc.
    - Fuse
		
*Custom function examples:*
- add an object to an object
- Find/define physical surfaces


Create a 'Network' (name still on for decision).
	Contains combinations of pieces (e.g. Cylinders, curves)
	Sequential pattern (add_straight -> add_curve)
	Automatically places new pieces in the right place
(at the end of the pipe)
	Fuses at the end to create a long pipe in desired shape.

Join Networks at junctions

### Requirements for pipes.py:
- libgmsh.so, from the GMSH SDK, available [here](http://gmsh.info/)
- numpy
- time to understand how it works, and deal with bugs, as it is incredibly raw
