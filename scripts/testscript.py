#! /usr/bin/env python

from scripts import ScriptException
import svnscripts
import autotools
import logging
import getopt
import sys
from ConfigParser import SafeConfigParser
from StringIO import StringIO
import os
import dougtest

defaultConfig = """
[testscript]
run-svn: yes
run-autotools: yes

[tests]
# set of properties in format
#  test<name>: <ctrlfile> <solutionfile>
# like
#  testN1: /tmp/test1/DOUG.dat /tmp/test1/correctanswer.dat
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

config = "\n".join([defaultConfig,
                    svnscripts.defaultConfig,
                    autotools.defaultConfig,
                    dougtest.defaultConfig])

for opt in opts:
    if opt[0]=="--logconf":
        logconfFileNames.append(opt[1])
    if opt[0]=="--conf":
        confFileNames.append(opt[1])
    if opt[0]=="-g":
        sys.stdout.write(config)
        sys.exit(0)

import logging.config
for fname in logconfFileNames:
    logging.config.fileConfig(fname)
else:
    logging.basicConfig()
    
LOG = logging.getLogger('testscript')

try:
    import time
    LOG.info("Running testscript at %s" % time.asctime())

    conf = SafeConfigParser()
    conf.readfp(StringIO(config))
    conf.read(confFileNames)

    # svn
    if conf.getboolean('testscript', 'run-svn'):
        svndoug = svnscripts.run(confFileNames)

    # autotools
    if conf.getboolean('testscript', 'run-autotools'):
        autotools.run(confFileNames)

    items = conf.items("tests")
    for name, value in items:
        if not name.startswith("test"): continue
        ctrlfname, solutionfname = tuple(value.strip().split(" ", 1))
        ctrlfname = os.path.abspath(ctrlfname)
        solutionfname = os.path.abspath(solutionfname)
        datadir = os.path.dirname(ctrlfname)
        LOG.info("Starting test '%s': %s" % (name, ctrlfname))
        test = dougtest.TestCase(datadir, ctrlfname, solutionfname, conf, 1)
        try:
            test.setUp()
            try:
                test.run()
            finally:
                test.tearDown()
        except ScriptException, e:
            LOG.info("Ended with error test '%s'" % (name, ))
            LOG.error("Error while running test: %s" % e)
        else:
            LOG.info("Ended successfully test '%s'" % (name, ))

    LOG.info("Ended testscript at %s" % time.asctime())

except ScriptException, e:
    LOG.critical("Error while running script: %s" % e)
    sys.exit(1)
except Exception, e:
    LOG.critical("Unknown exception occured: %s" % e)
    raise
