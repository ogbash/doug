#! /usr/bin/env python

from scripts import ScriptException
import autotools
import logging
import getopt
import sys
from ConfigParser import SafeConfigParser
from StringIO import StringIO

defaultConfig = """
[doxscript]
doxygen: /usr/bin/doxygen
docdir: doc
docdirhtml: %(docdir)s/html
doxygen-outfilename: doxygen.out
doxygen-errfilename: doxygen.err

[apache]
docdir: /var/www/docs
"""

logconfFileNames = []
confFileNames = []

def usage():
    sys.stderr.write("Usage:\n"
                     "\t%s --conf=<filename> --logconf=<filename>\n"
                     "\t%s -g (to generate default config to stdout)\n"
                     % (sys.argv[0], sys.argv[0]))

try:
    opts, extra = getopt.getopt(sys.argv[1:], "g", ["logconf=", "conf="])
except getopt.GetoptError:
    usage()
    sys.exit(1)

if len(extra)>0:
    usage()
    sys.exit(1)

for opt in opts:
    if opt[0]=="--logconf":
        logconfFileNames.append(opt[1])
    if opt[0]=="--conf":
        confFileNames.append(opt[1])
    if opt[0]=="-g":
        sys.stdout.write(defaultConfig)
        sys.stdout.write(autotools.defaultConfig)        
        sys.exit(0)

import logging.config
for fname in logconfFileNames:
    logging.config.fileConfig(fname)
else:
    logging.basicConfig()
    
LOG = logging.getLogger('doxscript')

try:
    import time
    LOG.info("Running doxscript at %s" % time.asctime())

    # autotools
    autotools.run(confFileNames)

    # doxygen 
    LOG.info("Preparing to run doxygen")
    import popen2
    import os

    conf = SafeConfigParser()
    conf.readfp(StringIO(defaultConfig))
    conf.read(confFileNames)

    cmd = conf.get('doxscript', 'doxygen')
    docdir = conf.get('doxscript', 'docdir')
    docdirhtml = conf.get('doxscript', 'docdirhtml')
    outfname = conf.get('doxscript', 'doxygen-outfilename')
    errfname = conf.get('doxscript', 'doxygen-errfilename')

    curdir = os.getcwd()
    os.chdir(docdir)
    LOG.debug('Changed directory to %s', docdir)

    LOG.info("Running doxygen")
    dx = popen2.Popen3("%s > %s 2> %s" % (cmd, outfname, errfname))
    dx.wait()

    dx_value = dx.poll()
    if dx_value != 0:
        raise ScriptException("Error occured while running doxygen (value=%d), "
                  "inspect output files (%s, %s) for error description." %
                  (dx_value, outfname, errfname))

    os.chdir(curdir)
    LOG.debug('Changed directory to %s', curdir)

    # copying to apache dir
    apache_docdir = conf.get('apache', 'docdir')

    LOG.info("Copying documentation to apache")
    dx = popen2.Popen3("cp -R %s/* %s > %s 2> %s" % (docdirhtml, apache_docdir, outfname, errfname))
    dx.wait()

    dx_value = dx.poll()
    if dx_value != 0:
        raise ScriptException("Error occured while copying documentation (value=%d), "
                  "inspect output files (%s, %s) for error description." %
                  (dx_value, outfname, errfname))

except ScriptException, e:
    LOG.critical("Error while running script: %s" % e)
    sys.exit(1)
except Exception, e:
    LOG.critical("Unknown exception occured: %s" % e)
    raise
