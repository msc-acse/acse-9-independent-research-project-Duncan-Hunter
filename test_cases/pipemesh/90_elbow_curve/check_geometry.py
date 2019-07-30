import os, sys

def test1():
    with open("src/pipe.msh", 'r') as testFile:
        content = testFile.readlines()
        assert(content[1][:3] == "2.2")
        assert(int(content[4][:3]) == 688)
        assert(len(content) == 4213)
    print(".msh file correct")

test1()