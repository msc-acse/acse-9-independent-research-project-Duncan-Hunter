import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pipe-gmsh",
    version="0.2.3",
    author="Duncan Hunter",
    author_email="dunchunter@hotmail.co.uk",
    description="A package for making pipe network meshes.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ImperialCollegeLondon/icferst_ACSE-IRP",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
    ],
)
