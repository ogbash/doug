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
# or contact the author (University of Tartu, Faculty of Computer Science, Chair
# of Distributed Systems, Liivi 2, 50409 Tartu, Estonia, http://dougdevel.org,
# mailto:info(at)dougdevel.org)

"""Module to help running autotools scripts.

Use run() procedure and pass list of configuration file names to it.
See defaultConfig member for default properties and reference.
"""

import popen2
import os
import time
import logging
from scripts import ScriptException

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
# whether to run configure in builddir directory
run-configure: yes
# whether to run make in builddir directory
run-make: yes

# time to sleep between autogen completion polls (in seconds)
autogen-waittime: 3

autogen-outfilename: autogen.out
autogen-errfilename: autogen.err

# command arguments to configure script
configure-arguments: 

# time to sleep between configure completion polls (in seconds)
configure-waittime: 6

configure-outfilename: configure.out
configure-errfilename: configure.err

# time to sleep between make completion polls (in seconds)
make-waittime: 6

make-outfilename: make.out
make-errfilename: make.err
"""
from StringIO import StringIO
defaultConfigFile = StringIO(defaultConfig)

if __name__=='__main__':
    logging.basicConfig()

LOG = logging.getLogger('autotools')

def run(configFileNames):
    "Configure, compile project and build executables."
    LOG.info('Running autotools scripts')

    from ConfigParser import SafeConfigParser
    conf = SafeConfigParser()
    conf.readfp(defaultConfigFile)
    conf.read(configFileNames)
    
    if conf.getboolean('autotools', 'run-autogen'):
        ag = Autogen(conf)
        ag.run()
    
    if conf.getboolean('autotools', 'run-configure'):
        cn = Configure(conf)
        cn.run()
    
    if conf.getboolean('autotools', 'run-make'):
        mk = Make(conf)
        mk.run()

    LOG.info('Exiting autotools scripts')

class Autogen:
    def __init__(self, conf):
        self.conf = conf

    def run(self):
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
                raise ScriptException("Error occured while running autogen.sh (value=%d), "
                          "inspect output files (%s, %s) for error description." %
                          (ag_value, outfname, errfname))
        
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
        try:
            LOG.debug('Changed directory to %s', builddir)

            cmd = os.path.join(srcdir, "configure")
            cf = popen2.Popen3("%s %s > %s 2> %s" % (cmd, argstr, outfname, errfname))
            while cf.poll() == -1:
                LOG.debug("Waiting other %d seconds for ./configure to complete" % waittime)
                time.sleep(waittime)

            cf_value = cf.poll()
            if cf_value != 0:
                raise ScriptException("Error occured while running configure (value=%d), "
                          "inspect output files (%s, %s) for error description." %
                          (cf_value, outfname, errfname))
        
        finally:
            os.chdir(curdir) 
            LOG.debug('Changed directory to %s', curdir)


class Make:
    "Run and process make output."

    def __init__(self, conf):
        self.conf = conf

    def run(self):
        if not self.conf.getboolean('autotools', 'run-make'):
            return
        
        LOG.info("Running make")

        curdir = os.getcwd()
        waittime = self.conf.getint('autotools', 'make-waittime')
        outfname = self.conf.get('autotools', 'make-outfilename')
        errfname = self.conf.get('autotools', 'make-errfilename')
        builddir = self.conf.get('autotools', 'builddir')

        os.chdir(builddir)
        try:
            LOG.debug('Changed directory to %s', builddir)

            cmd = "make"
            mk = popen2.Popen3("%s > %s 2> %s" % (cmd, outfname, errfname))
            while mk.poll() == -1:
                LOG.debug("Waiting other %d seconds for make to complete" % waittime)
                time.sleep(waittime)

            mk_value = mk.poll()
            if mk_value != 0:
                raise ScriptException("Error occured while running make (value=%d), "
                          "inspect output files (%s, %s) for error description." %
                          (mk_value, outfname, errfname))

        finally:
            os.chdir(curdir) 
            LOG.debug('Changed directory to %s', curdir)

if __name__=='__main__':
    logging.getLogger().setLevel(logging.DEBUG)
    print run([])
