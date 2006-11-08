#! /usr/bin/env python

# DOUG - Domain decomposition On Unstructured Grids
# Copyright (C) 1998-2006 Faculty of Computer Science, University of Tartu and
# Department of Mathematics, University of Bath
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# or contact the authors (University of Tartu, Faculty of Computer Science, Chair
# of Distributed Systems, Liivi 2, 50409 Tartu, Estonia, http://dougdevel.org,
# mailto:info(at)dougdevel.org)

from scripts import ScriptException
import svnscripts
import autotools
import logging
import getopt
import sys
from ConfigParser import SafeConfigParser
from StringIO import StringIO
import os
import unittest
import dougtest

defaultConfig = """
[testscript]
run-svn: yes
run-autotools: yes
save-tar: no
save-mysql: no

tar-file: results.tar

mysql-host: localhost
mysql-user:
mysql-password:
mysql-database:

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

    
    # construct tests
    testSuite = unittest.TestSuite()
    
    items = conf.items("tests")
    for name, value in items:
        if not name.startswith("test"): continue

        # read test data
        ctrlfname, solutionfname, testconfs = tuple(value.strip().split(" ", 2))
        ctrlfname = os.path.abspath(ctrlfname)
        solutionfname = os.path.abspath(solutionfname)
        datadir = os.path.dirname(ctrlfname)
        LOG.debug("Constructing test '%s'" % (name, ))
        LOG.debug("Control file: %s" % ctrlfname)
        LOG.debug("Correct solution file: %s" % (solutionfname,))

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

    # run tests
    testResults = [unittest._TextTestResult(unittest._WritelnDecorator(sys.stderr), False, 1)]

    saveTar = conf.getboolean("testscript", "save-tar")
    saveMysql = conf.getboolean("testscript", "save-mysql")
    
    if saveTar:
        import dougtesttar
        tarFileName = os.path.abspath(conf.get("testscript", "tar-file"))
        tarTestResult = dougtesttar.DougTarTestResult(tarFileName)
        testResults.append(tarTestResult)

    if saveMysql:
        import dougtestmysql
        mysqlHost = conf.get("testscript", "mysql-host")
        mysqlUser = conf.get("testscript", "mysql-user")
        mysqlPassword = conf.get("testscript", "mysql-password")
        mysqlDatabase = conf.get("testscript", "mysql-database")
        mysqlTestResult = dougtestmysql.DougMySQLTestResult(mysqlHost, mysqlUser, mysqlPassword, mysqlDatabase, conf)
        testResults.append(mysqlTestResult)

    testRunner = dougtest.TestRunner(testResults)

    try:
        testRunner.run(testSuite)
    finally:
        if saveTar:
            try:
                tarTestResult.close()
                LOG.debug("tar file closed")
            except Exception, e:
                LOG.error("Could not close tar file: %s" % e)

        if saveMysql:
            try:
                mysqlTestResult.close()
                LOG.debug("mysql closed")
            except Exception, e:
                LOG.error("Could not close mysql: %s" % e)

    LOG.info("Ended testscript at %s" % time.asctime())

except ScriptException, e:
    LOG.critical("Error while running script: %s" % e)
    sys.exit(1)
except Exception, e:
    LOG.critical("Unknown exception occured: %s" % e)
    raise
