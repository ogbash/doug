#! /usr/bin/env python

from scripts import ScriptException
import logging
import os
import re

defaultConfig="""
[dougtest]
# directory where doug_main and doug_aggr are located
dougbindir:

# MPI
mpiboot: lamboot
mpirun: mpirun
"""

LOG = logging.getLogger("dougtest")



class ControlFile:	
	re_assignment = re.compile("(\S+)\s+(.*)")

	def __init__(self, filename):
		self.name = filename
		self.options = {}
		self._parse()

	def _parse(self):
		f = open(self.name, 'r')
		try:
			for line in f:
				match = self.re_assignment.match(line)
				if match:
					self.options[match.group(1)] = match.group(2)
		finally:
			f.close()
			
	def save(self, filename):
		f = open(filename, "w")
		try:
			for key in self.options:
				f.write("%s %s\n" % (key, self.options[key]))
		finally:
			f.close()

class TestCase:
	def __init__(self, datadir, ctrlfname, solutionfname, conf, method):
		self.datadir = datadir
		self.ctrlfname = ctrlfname
		self.solutionfname = solutionfname
		self.conf = conf
		self.method = method
	
	def setUp(self):
		LOG.debug("Preparing before test")
		tmpdir = os.tempnam(None, 'doug-')
		os.mkdir(tmpdir)
		LOG.debug("Temporary directory %s created" % tmpdir)
		self.tmpdir = tmpdir
		try:
			# copy control file
			self.controlFile = ControlFile(self.ctrlfname)
			self.controlFile.options['plotting'] = '0'
			self.controlFile.options['method'] = str(self.method)
			self.testctrlfname = os.path.join(self.tmpdir,
							  os.path.basename(self.ctrlfname))
			self.controlFile.save(self.testctrlfname)

			# make hard link to all other files
			files = os.listdir(self.datadir)
			for fname in files:
				fullfname = os.path.join(self.datadir, fname)

				if fname == os.path.basename(self.testctrlfname):
					continue
				if os.path.isdir(fullfname):
					continue
				
				LOG.debug("Creating hard link for %s" % fullfname)
				os.link(fullfname,
					os.path.join(self.tmpdir, fname))

			# run lamboot (temporary solution for developing)
			mpiboot = self.conf.get("dougtest", "mpiboot")
			res = os.spawnlp(os.P_WAIT, mpiboot, mpiboot)
			if res:
				raise ScriptException("Error running %s (%d)" % (mpiboot, res))
		except:
			self._clean()
			raise
		
			
	def tearDown(self):
		LOG.debug("Cleaning after test")
		self._clean()

	def _clean(self):
		if hasattr(self, 'tmpdir') and self.tmpdir:
			os.system('rm -rf %s' % self.tmpdir)
			LOG.debug("Temporary directory %s deleted" % self.tmpdir)

	def run(self, result=None):
		LOG.debug("Running test")
		mpirun = self.conf.get("dougtest", "mpirun")
		dougbindir = os.path.abspath(self.conf.get("dougtest", "dougbindir"))
		main = os.path.join(dougbindir, "doug_main")
		
		curdir = os.getcwd()
		
		LOG.debug("Changing directory to %s" % self.tmpdir)
		os.chdir(self.tmpdir)
		try:
			LOG.debug("Running %s -np 1 %s -f %s" % (mpirun, main, self.testctrlfname))
			res = os.spawnlp(os.P_WAIT, mpirun, mpirun, "-np", "1",
					 main, "-f", self.testctrlfname)
			LOG.debug("Finished %s with code %d" % (mpirun, res))
			if res:
				raise ScriptException("Error running %s (%d)" % (mpirun, res))
		finally:
			LOG.debug("Changing directory to %s" % curdir)
			os.chdir(curdir)
