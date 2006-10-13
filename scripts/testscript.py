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
import unittest

defaultConfig = """
[testscript]
run-svn: yes
run-autotools: yes

[tests]
# set of properties in format
#  test<name>: <ctrlfile> <solutionfile> <commasep test confs>
# like
#  testN1: /tmp/test1/DOUG.dat /tmp/test1/correctanswer.dat conf1,conf2
# whereas configuration is specified as
#  [testconf_conf1]
#  solver=1,2
#  method=1
#  processors=1,9
"""

logconfFileNames = []
confFileNames = []
verbose = False

def usage():
    sys.stderr.write("Usage: %s\n"
                     "\t [--conf=<filename>]... [--logconf=<filename>]...\n"
                     "\t [-g] (to generate default config to stdout)\n"
                     "\t [-v] (verbose)\n"
                     % (sys.argv[0],))

try:
    opts, extra = getopt.getopt(sys.argv[1:], "gv", ["logconf=", "conf="])
except getopt.GetoptError:
    usage()
    sys.exit(1)

if len(extra)>0:
    usage()
    sys.exit(1)

config = "\n".join([defaultConfig,
                    svnscripts.defaultConfig,
                    autotools.defaultConfig,
                    dougtest.getDefaultConfig()])

for opt in opts:
    if opt[0]=="--logconf":
        logconfFileNames.append(opt[1])
    if opt[0]=="--conf":
        confFileNames.append(opt[1])
    if opt[0]=="-g":
        sys.stdout.write(config)
        sys.exit(0)
    if opt[0]=="-v":
        verbose = True

import logging.config
for fname in logconfFileNames:
    logging.config.fileConfig(fname)
else:
    logging.basicConfig()

if verbose:
    logging.getLogger().setLevel(logging.DEBUG)

LOG = logging.getLogger('testscript')

def generateTuples(*iters):
    elems = [[]]

    for i in range(0,len(iters)):
        appendelems = iters[i]
        elems = [x+[y] for x in elems for y in appendelems]
    
    return map(tuple, elems)

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

        # read test data
        ctrlfname, solutionfname, testconfs = tuple(value.strip().split(" ", 2))
        ctrlfname = os.path.abspath(ctrlfname)
        solutionfname = os.path.abspath(solutionfname)
        datadir = os.path.dirname(ctrlfname)
        LOG.info("Test '%s'" % (name, ))
        LOG.info("Control file: %s" % ctrlfname)
        LOG.info("Correct solution file: %s" % (solutionfname,))

        tarTestResult = dougtest.DougTarTestResult("/usr/svn/anytests.tar")
        import unittest
        testResults = [unittest._TextTestResult(unittest._WritelnDecorator(sys.stderr), False, 1),
                       tarTestResult]
        testRunner = dougtest.TestRunner(testResults)
        testSuite = unittest.TestSuite()

        # read test configurations
        testconfs = testconfs.split(",")
        for testconf in testconfs:
            testconfname = "testconf_%s" % testconf
            solvers = map(int, conf.get(testconfname, "solver").split(","))
            methods = map(int, conf.get(testconfname, "method").split(","))
            processors = map(int, conf.get(testconfname, "processors").split(","))

            testtuples = generateTuples(solvers, methods, processors)

            for testtuple in testtuples:
                test = dougtest.TestCase(name+"_"+testconf, datadir, ctrlfname, solutionfname, conf, *testtuple)
                testSuite.addTest(test)

    testRunner.run(testSuite)

    tarTestResult.close()
    LOG.debug("tar file closed")

    LOG.info("Ended testscript at %s" % time.asctime())

except ScriptException, e:
    LOG.critical("Error while running script: %s" % e)
    sys.exit(1)
except Exception, e:
    LOG.critical("Unknown exception occured: %s" % e)
    raise
