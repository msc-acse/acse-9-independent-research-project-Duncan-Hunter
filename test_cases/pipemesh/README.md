## pipemesh tests

This folder contains:
* Makefile: Generate a Python virtual environment which installs pipemesh and its dependencies. If you have the package installed already, this can be ignored. To create and activate the virtual environment:
```bash
$make
$source virtual_env/bin/activate
```

To remove the virtual environment:
```bash
$make clean_venv
```

To clean all subdirectories:
```bash
$make clean_all
```

* 90_elbow_curve/: A test case of a pipe with a 90 degree elbow. Call make to generate a the mesh file required and test it with check_geometry.py. The .mpml file can then be used with IC-FERST to run a short simulation to check if it runs.

```bash
$cd 90_elbow_curve
$make
$mpirun -n 8 path/to/icferst 3d_pipe_FEM.mpml
```

* change_radius/: A test case of a pipe with a changing radius (getting smaller). Call make to generate the mesh file. Decompes the mesh using fldecomp for 8 processors.

```bash
$cd change_radius
$make
$mpirun -n 8 /path/to/icferst 3d_pipe_FEM.mpml
```

* cylinder/: A test case of a cylinder. A test case of a cylinder, the same total length as the change_radius case (to look at pressure/velocity differences). Call make to generate the mesh and decompose for 8 processors.

```bash
$cd cylinder
$make
$mpirun -n 8 /path/to/icferst 3d_pipe_FEM.mpml
```

* t_junction/: A test case of a t_junction. Call make to generate the mesh and decompose for 8 processors.

```bash
$cd t_junction
$make
$mpirun -n 8 /path/to/icferst 3d_pipe_FEM.mpml
```

To remove the simulations in any folder, call
```bash
$make clean
```