#! /usr/bin/env python

import os, sys, popen2, time

class Repointerface:
	def __init__(self, debug):
		self.debug=debug
		#self.dougloc=dougloc

	def getcode(self, source):
		self.source=source
		#download the code
		print "getcode";
		print "source=",self.source
		if self.source.startswith("svn"):
			print "dougi asukoht svn-is"
			comm="svn co "+self.source
			print "reposource comm=",comm
			p=popen2.Popen3(comm, capturestderr=1)
			
			i=0
			while p.poll()==-1:
				if(i>10):
					print "Viga: ei suutnud 30 sekundiga koodi alla laadida. Valjun"
					sys.exit(3)
				i=i+1
				print "."
				time.sleep(3)
			#excode=pp.wait()
			excode=p.poll()
			errlines=p.childerr.readlines()
			if excode==0:
				print "Koodi allalaadimine onnestus, excode=",excode
				if len(errlines)>0:
					print "svn kask tagastas hoiatusi:"
					for line in errlines:
						print "svn hoiatus:"+line
				
			else:
				print "Koodi allalaadimine EI onnestunud, excode=",excode
				
				#code=p.wait()
				#print "svn tagastuskood:", code
				#erro=""
				#errlines=p.childerr.readlines()
				for line in errlines:
					line = line.strip()
					if(line.find("command not found")!=-1):
						#os viga
						print "Ei leia svn kasku. Kas svn on installeeritud?"
						print "Koodi allalaadimine ebaonnestus! Valjun."
						print "Veateade:", errlines
						#print "rida:"+line
						sys.exit(1)
					
					elif(line.startswith("svn:") or line.startswith("Unknown command:")):
						#svn viga
						print "Vigane svn syntaks voi vigased parameetrid."
						print "Koodi allalaadimine ebaonnestus. Valjun"
						print "Veateade:", errlines
						#print "rida:"+line
						sys.exit(1)
			
					else:
						print "nomatch errline while executing snv command:"+ line
		
			outlines=p.fromchild.readlines()
			for line in outlines:
				line = line.strip()
				print "out line:"+ line
			p.childerr.close()
			p.fromchild.close()
			del p
		else:
			#peaks kopima lihtsalt, eeldatavasti kohalikul kettal
			print "mujalt kui svnist kopeerimine pole veel implementeeritud. Kood tombamata. Valjun"
			sys.exit(1)
		
		#p.childerr.close()
		#p.fromchild.close()