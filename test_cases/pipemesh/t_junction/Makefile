SHELL = sh

pipe:	src/generate.py
	cd src;	\
	python3 generate.py;	\
	fldecomp -n 8 pipe;
	python3 check_geometry.py;


clean:
	rm src/pipe*
	rm -rf junction_flow_test_case*
	rm generated_mpml.mpml
	rm *.vtu
	rm *.stat
	rm *.csv