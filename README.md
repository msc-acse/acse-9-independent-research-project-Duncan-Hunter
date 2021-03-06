# ACSE-9 Independent Research Project
## Duncan Hunter (CID:01104228)

This is the submission repository, where the report can be found. Figures can be found in the Report/Figures directory. The code (with license, README, and pytest file) can be found in Tools/pipemesh. Some integration test cases for IC-FERST can be found in test_cases/pipemesh.

As the project is a pre-processing feature in the IC-FERST workflow, this repository is a mirror of a shared repository, where another person was working on post-processing. This repository only contains my work, and the commit history has been preserved. The shared repository can be found [here](https://github.com/ImperialCollegeLondon/icferst_ACSE-IRP) is a few commits behind this one. The python package created, *pipemesh* has a [public repository](https://github.com/Duncan-Hunter/pipemesh), which is a few pushes behind this one, and is used so readthedocs can generate and host documentation, as well as provide access to the source code for people using it.

This repository's structure:
```
+-- Tools/
|   +-- pipemesh/
|    // Python package directory, with pytest scripts and examples.py file.
|        +-- docs/
|        // Copy of what readthedocs uses in the public repository.
|        +-- pipemesh/
|        // Contains modules *pipes* and *pieces*.
|            +-- icferst/
|            // Contains *auto_mpml* module.
+-- test_cases/
|    +-- pipemesh/
|    // Contains examples of pipemesh's integration with IC-FERST
```
