import subprocess
import re
import copy
from StringIO import StringIO
import os

import doug
from doug.config import DOUGConfigParser, ControlFile

from scripts import ScriptException

import logging
LOG = logging.getLogger('doug')

_defaultConfig = None
def getDefaultConfig():
    global _defaultConfig
    if _defaultConfig is None:
        _defaultConfig = DOUGConfigParser(name="DOUG default config")
        _defaultConfig.addConfigContents(doug.execution_conf_tmpl)
    return _defaultConfig

def getDefaultControlFile(basedir):
    cf = ControlFile(contents=doug.DOUG_ctl_tmpl, basedir=basedir)
    cf.name = '<doug.DOUG_ctl_tmpl>'
    return cf

class DOUGExecution:

    def __init__(self, config, dougControls=None):
        self.workdir = os.path.abspath(config.get('doug', 'workdir'))
        self.config = DOUGConfigParser(name='DOUG execution', basedir=self.workdir)
        # default config
        self.config.addConfig(getDefaultConfig())
        self.config.addControlFile(getDefaultControlFile(self.workdir))
        
        # copy controls from control file
        if dougControls is not None:
            self.config.addControlFile(dougControls)
        # copy config
        self.config.addConfig(config)

        # output or other files, exception grabs it on exit
        self.files = []

        # how many test results are using this test, files are deleted only after last free() call
        self._inUse = 0

    def setUp(self):
        LOG.debug("Preparing testing environment")
        self.workdir = self.workdir
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
            controlFile = self.config.getControlFile(self.testctrlfname)
            controlFile.save(self.testctrlfname)
            self.files.append((self.testctrlfname, "Control file"))

            # run mpiboot
            mpibootname = self.config.get("doug", "mpiboot")
            outfilename = self.config.get("doug", "mpiboot-outfilename")
            errfilename = self.config.get("doug", "mpiboot-errfilename")

            if mpibootname:
                LOG.debug("Setting up mpi")
                mpiboot = subprocess.Popen("%s > %s 2> %s" % (mpibootname, outfilename, errfilename), shell=True)
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
                mpihalt = subprocess.Popen("%s > %s 2> %s" % (mpihaltname, outfilename, errfilename), shell=True)
                import time
                time.sleep(4) # lamhalt <=7.1.1 does not wait until whole universe is shut down
                res = mpihalt.wait()
                if res:
                    LOG.warn("Error running %s (%d)"
                         "inspect output files (%s, %s) for error description."
                         % (mpihaltname, res, outfilename, errfilename))
        except Exception, e:
            LOG.warn("Exception running mpihalt: %s" % e)

    def _clean(self):
        if not self.preserveOutput and not self.workdirExisted:
            os.system('rm -rf %s' % self.workdir)
            LOG.debug("Temporary directory %s deleted" % self.workdir)

    def acquire(self):
        self._inUse += 1

    def free(self):
        self._inUse -= 1
        if self._inUse == 0:
            self._clean()

    def run(self):
        return self.runDOUG()

    def runDOUG(self):
        LOG.debug("Running DOUG")
        nproc = self.config.getint('doug', 'nproc')
        solver = self.config.getint('doug-controls', 'solver')
        levels = self.config.getint('doug-controls', 'levels')
        method = self.config.getint('doug-controls', 'method')
        LOG.info("solver=%d, method=%d, levels=%d, nproc=%d" % (solver, method, levels, nproc))
        mpirun = self.config.get("doug", "mpirun")
        main = self.config.getpath("doug", "executable")
        errfname = self.config.getpath("doug", "errfilename")
        outfname = self.config.getpath("doug", "outfilename")
        solutionfname = self.config.getpath('doug-controls', 'solution_file')

        curdir = os.getcwd()

        result = DOUGConfigParser(self.config.defaults(), basedir=self.workdir)
        result.add_section('doug-result')
        try:
            LOG.debug("Changing directory to %s" % self.workdir)
            os.chdir(self.workdir)
            outf = open(outfname, "w")
            errf = open(errfname, "w")
            try:
                args = [mpirun, "-np", "%d"%nproc, main, "-f", self.testctrlfname, "-p"]
                LOG.debug("Running %s" % " ".join(args))
                doug = subprocess.Popen(args, stdout=outf, stderr=errf)
                    
                import time
                maxtime = self.config.getint('doug', 'max-time')
                for i in xrange(maxtime): # ~1 minute
                    time.sleep(1)
                    doug.poll()
                    if doug.returncode != None:
                        break
                else:
                    LOG.info("Terminating DOUG")
                    doug.terminate()
                    doug.wait()
                    
                value = doug.returncode
                LOG.debug("Finished %s with code %d" % (mpirun, value))
                self.files.append((outfname, "%s standard output" % mpirun))
                self.files.append((errfname, "%s standard error" % mpirun))
                result.setpath('doug-result', 'returnvalue', str(value))
                result.setpath('doug-result', 'outputfile', outfname)
                result.setpath('doug-result', 'errorfile', errfname)
                
                if value != 0:
                    se = ScriptException("Error occured while running doug (value=%d), "
                                             "inspect output files (%s, %s) for error description." %
                                             (value, outfname, errfname))
                    raise se

                if solutionfname and os.path.isfile(solutionfname):
                    result.setpath('doug-result', 'solutionfile', solutionfname)
                if solutionfname and os.path.isfile('aggr1.txt'):
                    result.setpath('doug-result', 'fineaggrsfile', 'aggr1.txt')
                if solutionfname and os.path.isfile('aggr2.txt'):
                    result.setpath('doug-result', 'coarseaggrsfile', 'aggr2.txt')
                    
                files = os.listdir(self.workdir)
                files = filter(lambda name: name.startswith('prof.'), files)
                if files:
                    result.setpath('doug-result', 'profilefile', files[0])
                
                # compare answers
                
            finally:
                outf.close()
                errf.close()
                LOG.debug("Changing directory to %s" % curdir)
                os.chdir(curdir)
        except ScriptException, e:
            for fn in self.files:
                e.addFile(*fn)
            raise

        return result

    def __str__(self):
        return "%s: solver=%d, method=%d, processors=%d" % \
               (self.testname, self.solver, self.method, self.nproc)
