#! /usr/bin/env python

# import packages from modules directory
import sys, os
pathname = os.path.dirname(sys.argv[0])
sys.path.append(os.path.join(os.path.abspath(pathname), 'modules'))

if len(sys.argv) < 3:
    print "Usage: %s <file1> <file2>" % sys.argv[0]
    sys.exit(1)

import doug.utils
import doug.formats

format = doug.formats.Fortran()

file1 = open(sys.argv[1])
try:
    vec1 = format.readVectorFromFile(file1)
    print "Solution 1 size: %d" % len(vec1)
finally:
    file1.close()

file2 = open(sys.argv[2])
try:
    vec2 = format.readVectorFromFile(file2)
    print "Solution 2 size: %d" % len(vec2)
finally:
    file1.close()

diff = doug.utils.diffVectors(vec1, vec2)
print "Difference %s" % diff

sys.exit(0)
