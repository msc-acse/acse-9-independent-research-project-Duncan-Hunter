import os, sys

def test1():
    with open("src/pipe.msh", 'rb') as testFile:
        content = testFile.readlines()
        assert(content[1][:3] == b"2.2")
        assert(int(content[5][:3]) == 688)
        assert(len(content) == 386)
    print(".msh file correct")

test1()