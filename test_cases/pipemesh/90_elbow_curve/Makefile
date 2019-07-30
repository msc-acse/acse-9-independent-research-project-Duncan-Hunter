SHELL = sh

pipe:	src/generate.py
	cd src;	\
	python3 generate.py;
	python3 check_geometry.py

clean:
	rm src/pipe*