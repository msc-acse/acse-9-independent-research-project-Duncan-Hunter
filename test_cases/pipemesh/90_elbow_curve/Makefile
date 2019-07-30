SHELL = sh

pipe:	src/generate.py
	cd src;	\
	python3 generate.py;
	python3 check_geometry.py

clean:
	rm src/pipe*
	rm -rf 3d_pipe_elbow_test_case*