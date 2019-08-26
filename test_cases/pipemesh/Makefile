SHELL = sh

venv:	clean
	mkdir virtual_env
	(	\
		python3 -m venv virtual_env;	\
		. virtual_env/bin/activate;	\
		python3 -m pip install --upgrade pip;	\
		python3 -m pip install scipy;	\
		python3 -m pip install -i https://test.pypi.org/simple/ pipemesh;	\
	)

clean_venv:
	rm -rf virtual_env

clean_all:
	cd cylinder && make clean
	cd 90_elbow_curve && make clean
	cd 90_elbow_miter && make clean
	cd change_radius && make clean
	cd t_junction && make clean