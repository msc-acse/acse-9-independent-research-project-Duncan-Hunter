import os, sys

def test1():
    with open("src/pipe.msh", 'rb') as testFile:
        content = testFile.readlines()
        print(len(content))
        assert(content[1][:3] == b"2.2")
        assert(int(content[5][:3]) == 642)
        assert(len(content) == 387)
    print(".msh file correct")

test1()