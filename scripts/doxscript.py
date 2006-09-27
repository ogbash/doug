#! /usr/bin/env python

import autotools
import logging
import getopt
import sys

logconfFileName = None
confFileName = None

def usage():
    sys.stderr.write("Usage:\n\t%s --logconf=<filename> --conf=<filename>\n" % sys.argv[0])

opts, extra = getopt.getopt(sys.argv[1:], "", ["logconf=", "conf="])

if len(extra)>0:
    usage()
    sys.exit(1)

for opt in opts:
    if opt[0]=="--logconf":
        logconfFileName = opt[1]
    if opt[0]=="--conf":
        confFileName = opt[1]

if not confFileName or not logconfFileName:
    usage()
    sys.exit(1)

import logging.config
logging.config.fileConfig(logconfFileName)
LOG = logging.getLogger('doxscript')

import time
LOG.info("Running doxscript at %s" % time.asctime())

# autotools
autotools.run([confFileName])

# doxygen 
LOG.info("Preparing to run doxygen")
import popen2
import os

from ConfigParser import SafeConfigParser
conf = SafeConfigParser()
conf.read([confFileName])

cmd = conf.get('doxscript', 'doxygen')
docdir = conf.get('doxscript', 'docdir')
docdirhtml = conf.get('doxscript', 'docdirhtml')
outfname = os.path.abspath(conf.get('doxscript', 'doxygen-outfilename'))
errfname = os.path.abspath(conf.get('doxscript', 'doxygen-errfilename'))

curdir = os.getcwd()
os.chdir(docdir)
LOG.debug('Changed directory to %s', docdir)

LOG.info("Running doxygen")
dx = popen2.Popen3("%s > %s 2> %s" % (cmd, outfname, errfname))
dx.wait()

dx_value = dx.poll()
if dx_value != 0:
    LOG.error("Error occured while running doxygen (value=%d), "
              "inspect output files (%s, %s) for error description." %
              (dx_value, outfname, errfname))
    sys.exit(1)

os.chdir(curdir)
LOG.debug('Changed directory to %s', curdir)

# copying to apache dir
apache_docdir = conf.get('apache', 'docdir')

LOG.info("Copying documentation to apache")
dx = popen2.Popen3("cp -R %s/* %s > %s 2> %s" % (docdirhtml, apache_docdir, outfname, errfname))
dx.wait()

dx_value = dx.poll()
if dx_value != 0:
    LOG.error("Error occured while copying documentation (value=%d), "
              "inspect output files (%s, %s) for error description." %
              (dx_value, outfname, errfname))
    sys.exit(1)
