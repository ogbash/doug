"""Module to help running autotools scripts.

Use run() procedure and pass list of configuration file names to it.
See defaultConfig member for default properties and reference.
"""

import popen2
import os
import time
import logging

defaultConfig="""
# This are default (ie reference) values for autotools helper package.
# See python ConfigParser package to parse this file.

[autotools]

# top source directory, ie where configure is located
srcdir: .
# top build directory, ie where configure is run (this may be different from srcdir)
builddir: .

# whether to run autogen.sh in srcdir directory
run-autogen: yes
# time to sleep between autogen completion polls (in seconds)
autogen-waittime: 3

autogen-outfilename: autogen.out
autogen-errfilename: autogen.err

# command arguments to configure script
configure-arguments: 

# time to sleep between autogen completion polls (in seconds)
configure-waittime: 6

configure-outfilename: configure.out
configure-errfilename: configure.err
"""
from StringIO import StringIO
defaultConfigFile = StringIO(defaultConfig)

if __name__=='__main__':
    logging.basicConfig()

LOG = logging.getLogger('autotools')

def run(configFileNames):
    "Compile project and build executables."
    LOG.info('Running autotools scripts')

    from ConfigParser import SafeConfigParser
    conf = SafeConfigParser()
    conf.readfp(defaultConfigFile)
    conf.read(configFileNames)
    
    ag = Autogen(conf)
    val = ag.run()
    if val: return val # error occured
    
    cn = Configure(conf)
    val = cn.run()
    if val: return val # error occured
    
    #mk = Make(conf)
    #val = mk.run()
    #if val: return val # error occured
    LOG.info('Exiting autotools scripts')

class Autogen:
    def __init__(self, conf):
        self.conf = conf

    def run(self):
        if not self.conf.getboolean('autotools', 'run-autogen'):
            return

        LOG.info('Running autogen')
        
        curdir = os.getcwd()
        waittime = self.conf.getint('autotools', 'autogen-waittime')
        outfname = self.conf.get('autotools', 'autogen-outfilename')
        errfname = self.conf.get('autotools', 'autogen-errfilename')
        srcdir = self.conf.get('autotools', 'srcdir')

        os.chdir(srcdir)
        LOG.debug('Changed directory to %s', srcdir)
        
        try:
            ag = popen2.Popen3('./autogen.sh > %s 2> %s' % (outfname, errfname))
            while ag.poll() == -1:
                LOG.debug("Waiting other %d seconds for ./autogen.sh to complete" % waittime)
                time.sleep(waittime)

            ag_value = ag.poll()
            if ag_value != 0:
                LOG.error("Error occured while running autogen.sh (value=%d), "
                          "inspect output files (%s, %s) for error description." %
                          (ag_value, outfname, errfname))
        
            return ag_value
        finally:
            os.chdir(curdir) 
            LOG.debug('Changed directory to %s', curdir)
               
class Configure:
    "Run and process configure output."

    def __init__(self, conf):
        self.conf = conf

    def run(self):
        LOG.info("Running configure")

        curdir = os.getcwd()
        waittime = self.conf.getint('autotools', 'configure-waittime')
        outfname = self.conf.get('autotools', 'configure-outfilename')
        errfname = self.conf.get('autotools', 'configure-errfilename')
        srcdir = os.path.abspath(self.conf.get('autotools', 'srcdir'))
        builddir = self.conf.get('autotools', 'builddir')
        argstr = self.conf.get('autotools', 'configure-arguments')

        os.chdir(builddir)
        LOG.debug('Changed directory to %s', srcdir)

        cmd = os.path.join(srcdir, "configure")
        
        try:
            cf = popen2.Popen3("%s %s > %s 2> %s" % (cmd, argstr, outfname, errfname))
            while cf.poll() == -1:
                LOG.debug("Waiting other %d seconds for ./configure to complete" % waittime)
                time.sleep(waittime)

            cf_value = cf.poll()
            if cf_value != 0:
                LOG.error("Error occured while running configure (value=%d), "
                          "inspect output files (%s, %s) for error description." %
                          (cf_value, outfname, errfname))
        
            return cf_value
        finally:
            os.chdir(curdir) 
            LOG.debug('Changed directory to %s', curdir)


class Make:
    "Run and process make output."

    def __init__(self, conf):
        self.conf = conf

    def run(self):
        LOG.info("Running make")        

if __name__=='__main__':
    logging.getLogger().setLevel(logging.DEBUG)
    print run([])
