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
logging.basicConfig(level=logging.INFO)
import getopt
import sys
from ConfigParser import SafeConfigParser
from StringIO import StringIO
import os
import unittest
import doug.test as dougtest
import doug.execution
from doug.config import DOUGConfigParser
import subprocess

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
#  fine_method=1,2
#  num_iters=3
#  num_subdomains=1,2
#  processors=1,4
#  overlaps=-1
#  smoothers=0,1
#  executables=doug_geom,doug_aggr
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
            tname = name[4:]

            # read test data
            ctrlfname, solutionfname, testconfs = tuple(value.strip().split(" ", 2))
            ctrlfname = os.path.abspath(ctrlfname)
            solutionfname = os.path.abspath(solutionfname)
            datadir = os.path.dirname(ctrlfname)
            LOG.debug("Constructing test '%s'" % (tname, ))
            LOG.debug("Control file: %s" % ctrlfname)
            LOG.debug("Correct solution file: %s" % (solutionfname,))

            # read test configurations
            testconfs = testconfs.split(",")
            for testconf in testconfs:
                testconfname = "testconf_%s" % testconf
                testname = "%s_%s" % (tname, testconf)
                solvers = map(int, conf.get(testconfname, "solver").split(","))
                methods = map(int, conf.get(testconfname, "method").split(","))
                levels = map(int, conf.get(testconfname, "levels").split(","))
                fine_methods = map(int, conf.get(testconfname, "fine_methods").split(","))
                nsubdomains = map(int, conf.get(testconfname, "num_subdomains").split(","))
                processors = map(int, conf.get(testconfname, "processors").split(","))
                executables = conf.get(testconfname, "executables").split(",")
                overlaps = map(int, conf.get(testconfname, "overlaps").split(","))
                smoothers = map(int, conf.get(testconfname, "smoothers").split(","))

                testtuples = generateTuples(solvers, methods, levels, fine_methods, nsubdomains, processors, executables, overlaps, smoothers)

                for solver,method,level,fine_method, nsubdomain,nproc,executable,overlap,smoother in testtuples:
                    dougControlFile = doug.execution.ControlFile(filename=ctrlfname, basedir=os.path.dirname(ctrlfname))
                    dougConfig = DOUGConfigParser(name='DOUG execution parameters')
                    # set/copy doug configuration from tests configuration
                    dougConfig.add_section('doug')
                    for name,value in conf.items('doug', raw=True):
                        dougConfig.set('doug', name, value)
                    dougConfig.add_section('doug-controls')
                    for name,value in conf._sections['doug-controls'].items():
                        dougConfig.set('doug-controls', name, value)
                    dougConfig.add_section('doug-tests')
                    dougConfig.set('doug-tests', 'csolution_file', solutionfname)
                    # set test configuration
                    dougConfig.set('doug-controls', 'solver', str(solver))
                    dougConfig.set('doug-controls', 'method', str(method))
                    dougConfig.set('doug-controls', 'levels', str(level))
                    dougConfig.set('doug-controls', 'fine_method', str(fine_method))
                    dougConfig.set('doug-controls', 'num_subdomains', str(nsubdomain))
                    dougConfig.set('doug', 'nproc', str(nproc))
                    dougConfig.set('doug', 'executable', executable)
                    dougConfig.set('doug-controls', 'overlap', str(overlap))
                    dougConfig.set('doug-controls', 'smoothers', str(smoother))

                    # create DOUG execution object
                    execution = doug.execution.DOUGExecution(dougConfig, dougControlFile)
                    test = dougtest.TestCase(execution, testname)
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

    # pick up git version if not specified
    if not conf.has_option('dougtest', 'info-git') or \
       not conf.get('dougtest', 'info-git'):
        subp = subprocess.Popen(["git","log","--pretty=%H","HEAD^..HEAD"], stdout=subprocess.PIPE)
        subp.wait()
        if subp.returncode==0:
            gitversion = subp.communicate()[0].strip()
            conf.set('dougtest', 'info-git', gitversion)
            LOG.info("Set git version to %s" % gitversion)

    # set hostname
    if not conf.has_option('dougtest', 'info-server') or \
       not conf.get('dougtest', 'info-server'):
        import socket
        conf.set('dougtest', 'info-server', socket.gethostname())
        LOG.info("Set hostname to %s" % socket.gethostname())

    # try to recognise fortran compile through mpif90
    if not conf.has_option('dougtest', 'info-fc') or \
       not conf.get('dougtest', 'info-fc'):
        subp = subprocess.Popen(["mpif90","-showme"], stdout=subprocess.PIPE)
        subp.wait()
        if subp.returncode==0:
            fc = subp.communicate()[0].split(" ")[0]
            conf.set('dougtest', 'info-fc', fc)
            LOG.info("Set fc to %s" % fc)


    # create test result objects
    testResults = [unittest._TextTestResult(unittest._WritelnDecorator(sys.stderr), False, 1)]

    saveTar = conf.getboolean("testscript", "save-tar")
    saveMysql = conf.getboolean("testscript", "save-mysql")

    if saveTar:
        import doug.testtar as dougtesttar
        tarFileName = os.path.abspath(conf.get("testscript", "tar-file"))
        tarTestResult = dougtesttar.DougTarTestResult(tarFileName, conf)
        testResults.append(tarTestResult)

    if saveMysql:
        import doug.testmysql as dougtestmysql
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
