SHELL = sh

pipe:	src/generate.py
	cd src;	\
	python3 generate.py;	\
	fldecomp -n 8 pipe;
	python3 check_geometry.py

clean:
	rm src/pipe*
	rm *.mpml
	rm -rf 3d_pipe_elbow_test_case*
	rm *.vtu
	rm *.stat
	rm *.csv
