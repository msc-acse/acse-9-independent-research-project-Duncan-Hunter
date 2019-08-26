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

clean:
	rm -rf virtual_env