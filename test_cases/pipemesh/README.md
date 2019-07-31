## pipemesh tests

This folder contains:
* Makefile: Generate a Python virtual environment which installs pipemesh and its dependencies. If you have the package installed already, this can be ignored. To create and activate the virtual environment:
```
$make
$source virtual_env/bin/activate
```
* 90_elbow_curve/: A test case of a pipe with a 90 degree elbow. Call make to generate a the mesh file required and test it with check_geometry.py. The .mpml file can then be used with IC-FERST to run a short simulation to check if it runs (serial)
```
$path/to/icferst 3d_pipe_FEM.mpml
```

* change_radius/: A test case of a pipe with a changing radius (getting smaller).

* cylinder/: A test case of a cylinder.

* t_junction/: A test case of a t_junction.