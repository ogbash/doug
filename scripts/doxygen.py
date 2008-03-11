#! /usr/bin/env python

#: Simple script that generates DOUG documentation and copies it to the
#: corresponding www directory.

import sys
import os

log = open("/srv/doug/output/doxygen.out", "w+")
log.write("ARGS: %s\n" % sys.argv)

rev = sys.argv[2]

orig_dir = os.getcwd()
os.chdir('/srv/doug/doug-doc')
try:
    ret = os.spawnlp(os.P_WAIT, 'svn', 'svn', 'up', '-r', rev)
    log.write("svn up returned %d\n" % ret)

    os.chdir('doc')
    ret = os.spawnl(os.P_WAIT, '/usr/local/bin/doxygen', 'doxygen')
    log.write("doxygen returned %d\n" % ret)

    ret = os.system('cp html/* /srv/www/dougdevel.org/doug/docs')
    log.write("copy returned %d\n" % ret)
finally:
    os.chdir(orig_dir)

