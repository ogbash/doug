import popen2
import re
import copy
from StringIO import StringIO
import os

import doug
from doug.config import DOUGConfigParser

from scripts import ScriptException

import logging
LOG = logging.getLogger('doug')

def getDefaultConfigContents():
    return doug.execution_conf_tmpl

class ControlFile:
    re_assignment = re.compile("(\S+)\s+(.*)")

    def __init__(self, filename=None, contents=None):
        self.name = filename
        self.options = {}

        if not contents:
            f = open(self.name, 'r')
            try:
                self._parse(f)
            finally:
                f.close()
        else:
            self._parse(StringIO(contents))

    def _parse(self, f):
        for line in f:
            match = self.re_assignment.match(line)
            if match:
                self.options[match.group(1)] = match.group(2)
            
    def save(self, filename):
        f = open(filename, "w")
        try:
            for key in self.options:
                f.write("%s %s\n" % (key, self.options[key]))
        finally:
            f.close()

_defaultControlFile = ControlFile(contents=doug.DOUG_ctl_tmpl)

class DOUGExecution:

    def __init__(self, config):
        self.config = DOUGConfigParser(name='DOUG execution')
        self.config.addConfigContents(doug.execution_conf_tmpl)
        self.config.addConfig(config)
        
        # create control file object
        self.controlFile = copy.deepcopy(_defaultControlFile)
        for option, value in self.config.items('doug-controls', nodefaults=True):
            self.controlFile.options[option] = value
        
        # output or other files, exception grabs it on exit
        self.files = []

    def setUp(self):
        LOG.debug("Preparing testing environment")
        self.workdir = self.config.get("doug", "cwd")
        if os.path.isdir(self.workdir):
            self.workdirExisted = True
        else:
            self.workdirExisted = False
            os.mkdir(self.workdir)
            LOG.debug("Working directory %s created" % self.workdir)
        self.preserveOutput = self.config.getboolean("doug", "preserveOutput")
        try:
            # create control file
            self.testctrlfname = os.path.abspath(os.path.join(self.workdir, 'DOUG-exec.ctl'))
            self.controlFile.save(self.testctrlfname)
            self.files.append((self.testctrlfname, "Control file"))

            # run mpiboot
            mpibootname = self.config.get("doug", "mpiboot")
            outfilename = self.config.get("doug", "mpiboot-outfilename")
            errfilename = self.config.get("doug", "mpiboot-errfilename")

            if mpibootname:
                LOG.debug("Setting up mpi")
                mpiboot = popen2.Popen3("%s > %s 2> %s" % (mpibootname, outfilename, errfilename))
                res = mpiboot.wait()
                if res:
                    raise ScriptException("Error running %s (%d)"
                                          "inspect output files (%s, %s) for error description."
                                          % (mpibootname, res, outfilename, errfilename))
        except:
            self._clean()
            raise
        
            
    def tearDown(self):
        try:
            mpihaltname = self.config.get("doug", "mpihalt")
            if mpihaltname:
                outfilename = self.config.get("doug", "mpihalt-outfilename")
                errfilename = self.config.get("doug", "mpihalt-errfilename")          
                LOG.debug("Shutting down mpi")
                mpihalt = popen2.Popen3("%s > %s 2> %s" % (mpihaltname, outfilename, errfilename))
                import time
                time.sleep(4) # lamhalt <=7.1.1 does not wait until whole universe is shut down
                res = mpihalt.wait()
                if res:
                    LOG.warn("Error running %s (%d)"
                         "inspect output files (%s, %s) for error description."
                         % (mpihaltname, res, outfilename, errfilename))
        except Exception, e:
            LOG.warn("Exception running mpihalt: %s" % e)
        
        LOG.debug("Cleaning after test")
        self._clean()


    def _clean(self):
        if not self.preserveOutput and hasattr(self, 'tmpdir') and self.tmpdir:
            # dangerous while developing -- os.system('rm -rf %s' % self.tmpdir)
            LOG.debug("Temporary directory %s deleted" % self.tmpdir)

    def run(self):
        self.setUp()
        try:
            self.runDOUG()
        finally:
            self.tearDown()

    def runDOUG(self):
        LOG.debug("Running DOUG")
        solver = self.config.getint('doug-controls', 'solver')
        nproc = self.config.getint('doug', 'nproc')
        levels = self.config.getint('doug-controls', 'levels')
        method = self.config.getint('doug-controls', 'method')
        LOG.info("solver=%d, method=%d, levels=%d, nproc=%d" % (solver, method, levels, nproc))
        mpirun = self.config.get("doug", "mpirun")
        main = self.config.getpath("doug", "executable")
        errfname = self.config.getpath("doug", "errfilename")
        outfname = self.config.getpath("doug", "outfilename")

        curdir = os.getcwd()
        
        try:
            LOG.debug("Changing directory to %s" % self.workdir)
            os.chdir(self.workdir)
            try:
                LOG.debug("Running %s -np %d %s -f %s -p" % (mpirun, nproc, main, self.testctrlfname))
                doug = popen2.Popen3('%s -np %d %s -f %s -p > %s 2> %s'%
                                     (mpirun, nproc, main, self.testctrlfname,
                                      outfname, errfname)
                                     )
                import time
                #time.sleep(1) # without this mpirun somehow gets HUP signal

                value = doug.wait()
                LOG.debug("Finished %s with code %d" % (mpirun, value))
                self.files.append((outfname, "%s standard output" % mpirun))
                self.files.append((errfname, "%s standard error" % mpirun))

                if value != 0:
                    se = ScriptException("Error occured while running doug (value=%d), "
                                             "inspect output files (%s, %s) for error description." %
                                             (value, outfname, errfname))
                    raise se

                # compare answers
            finally:
                LOG.debug("Changing directory to %s" % curdir)
                os.chdir(curdir)
        except ScriptException, e:
            for fn in self.files:
                e.addFile(*fn)
            raise

    def __str__(self):
        return "%s: solver=%d, method=%d, processors=%d" % \
               (self.testname, self.solver, self.method, self.nproc)
