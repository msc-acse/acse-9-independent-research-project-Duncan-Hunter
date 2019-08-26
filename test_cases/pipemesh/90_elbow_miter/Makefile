SHELL = sh

pipe:
	cd src;	\
	python3 generate.py;	\
	fldecomp -n 8 pipe;

clean:
	rm src/pipe*
	rm -rf 3d_pipe_miter_test_case*
	rm *.mpml
	rm *.vtu
	rm *.stat
	rm *.csv
