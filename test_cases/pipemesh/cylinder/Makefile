SHELL = sh

pipe:	src/generate.py
	cd src;	\
	python3 generate.py;	\
	fldecomp -n 8 pipe;
	python3 check_geometry.py;

clean:
	rm src/pipe*
	rm -rf 3d_pipe_cylinder_test_case*
	rm *.mpml
	rm -rf *.pvtu
	rm *.stat
	rm *.csv