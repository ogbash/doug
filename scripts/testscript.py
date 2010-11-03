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
import doug.execution
from doug.config import DOUGConfigParser

defaultConfig = """
[testscript]
run-svn: yes
run-autotools: yes
run-tests: yes
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
#  testN1: tests/test1/DOUG.dat tests/test1/correctanswer.dat conf1,conf2
# whereas configuration is specified as
#  [testconf_conf1]
#  solver=1,2
#  method=1
#  levels=1,2
#  processors=1,4
#  executables=doug_main,doug_aggr
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

conf = SafeConfigParser()
conf.readfp(StringIO(config))
conf.read(confFileNames)

conf.set("DEFAULT", "cwd", os.getcwd())

def generateTuples(*iters):
    elems = [[]]

    for i in range(0,len(iters)):
        appendelems = iters[i]
        elems = [x+[y] for x in elems for y in appendelems]
    
    return map(tuple, elems)

def main(testResults):    
    # svn
    if conf.getboolean('testscript', 'run-svn'):
        svndoug = svnscripts.run(confFileNames, conf.defaults())

    # autotools
    if conf.getboolean('testscript', 'run-autotools'):
        autotools.run(confFileNames, conf.defaults())

    
    # construct tests

    if conf.getboolean('testscript', 'run-tests'):
        testSuite = unittest.TestSuite()
        items = conf.items("tests")
        for name, value in items:
            if not name.startswith("test"): continue
            name = name[4:]

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
                levels = map(int, conf.get(testconfname, "levels").split(","))
                processors = map(int, conf.get(testconfname, "processors").split(","))
                executables = conf.get(testconfname, "executables").split(",")

                testtuples = generateTuples(solvers, methods, levels, processors, executables)

                for testtuple in testtuples:
                    dougControlFile = doug.execution.ControlFile(filename=ctrlfname)
                    dougConfig = DOUGConfigParser(name='DOUG execution parameters')
                    # set/copy doug configuration from tests configuration
                    dougConfig.add_section('doug')
                    for name,value in conf.items('doug', raw=True):
                        dougConfig.set('doug', name, value)
                    dougConfig.add_section('doug-controls')
                    for name,value in conf._sections['doug-controls'].items():
                        dougConfig.set('doug-controls', name, value)
                    execution = doug.execution.DOUGExecution(dougConfig, dougControlFile)
                    test = dougtest.TestCase(execution)
                    #resultConfig=execution.run()
                    
                    #test = dougtest.MPITestCase(name+"_"+testconf, datadir, ctrlfname, solutionfname, conf, *testtuple)
                    testSuite.addTest(test)

        # run tests
        testRunner = dougtest.TestRunner(testResults)
        testRunner.run(testSuite)

# ----------------------------------------

import time
LOG.info("Running testscript at %s" % time.asctime())

try:
    # It is not a correct behaviour to find SVN revision here, but currently compilation and testing
    # data are saved to the same MySQL table. So, MySQLTestResult object represents both processes,
    # whereas it needs to know SVN revision on it's creation as tests running process.
    # In the future, compilation should recieve object and table of its own!
    if not conf.get('dougtest', 'info-svn'):
        try:
            revision = svnscripts.getRevision(conf.get('autotools','srcdir'))
            conf.set('dougtest', 'info-svn', str(revision))
        except ValueError, e:
            LOG.warn("Cannot parse svn revision: %s" % e)

    # create test result objects
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

    try:
        try:
            # run
            main(testResults)
        except ScriptException, e:
            LOG.critical("Error while running script: %s" % e)
            # again, this is compilation or some setup error, not actual test error
            # this must code be changed when compilation recieves its own table
            if saveMysql:
                mysqlTestResult.addError(None, (e.__class__,e,None))
            raise
        except Exception, e:
            LOG.critical("Unknown exception occured: %s" % e)
            # again, this is compilation or some setup error, not actual test error
            # this must code be changed when compilation recieves its own table
            if saveMysql:
                mysqlTestResult.addError(None, (e.__class__,e,None))
            raise

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

finally:
    LOG.info("Ended testscript at %s" % time.asctime())
