#! /usr/bin/env python

import os, sys, popen2, time, shutil, getopt
#import interfaceAB, repointerface, errorhandler, comparator, watchdog
import repointerface, errorhandler, interfaceAB, watchdog
import subprocess   #alates v2.4

__DEBUG=True

def log(s):
	if(DEBUG==true):
		print s

def usage():
	print """usage: dougtest.py [options]
Start DOUG (regression) testing.
Options:
  -h or --help			Print this usage information
  -p<num>			Start each testcase using this 
  				list of process numbers.
				(default 1,4,7). Separate numbers 
				with commas(,) or (-).
  -d				Debug
  --testplace=<local|grid>	Specify where these tests should 
  				run(current default local)
  --nogridexec			Not implemented yet
  --nogridsave			Not implemented yet
  --nocompilation		Do not compile(make) DOUG. 
  				Assume compilation is already done 
				in ./doug dir.
  --makedeffile=<abslocalpath>	For specifying main Make.def file 
  				location(should be abs), which normally
				resides in doug/trunk/src. Script copies it 
				from <abslocalpath> to that dir. If not 
				specified, script will look for 
				./doug/trunk/src/Make.def or 
				./doug/trunk/src/Make.def.example(renames it)
				file. 
				If none is found, it will exit.
  --notest			Perform no testing activity.(no call 
  				to watchdog)
  --nocodegetting		No doug code download from repos.
	
	
     Default: dougtest.py -p1,4,7 --testplace=local 
	"""

class Dougtest:
	"DOUG testimise peamoodul"
	#self.__DEBUG=True
	#self.ABInterface
	#Repointerface=null
	#repoi=null
	#self.Errorhandler
	#self.Comparator
	#self.Watchdog
	#set repoi
	def log(self, s):
		if(self.DEBUG==true):
			print s
			
	

	def __init__(self):
		print "dougtest init";
		#default parameters
		self.id=time.strftime("%Y-%m-%d-%H-%M-%S")
		print "aja formaat ehk testi id:"+self.id
		#curtime=time.localtime(None)
		#self.id=str(curtime[0])+"-"+str(curtime[1])+"-"+str(curtime[2])+" "+str(curtime[3])+" "+str(curtime[4])
		self.nogrid=True #testplace= hoopis?
		self.testplace="local"
		self.DEBUG=False
		self.nocompilation=False
		self.dougaddr="svn://kheiron.at.mt.ut.ee/doug"
		self.examplesaddr="svn://kheiron.at.mt.ut.ee/doug_examples"
		self.makedefrelpath="doug/trunk/src" #relative to/in dougdir
		self.makedeffile=None
		#self.makedefrelpath="doug/trunk/src"
		self.makefilerelpath="doug/trunk"
		self.dougdir=self.getlastcomp(self.dougaddr)
		self.examplesdir="doug_examples"
		self.nocodegetting=False
		self.notest=False
		self.dougexerelpath="doug/trunk/bin"
		self.dougexefilename="DOUG_main"
		self.numberofcopies=[1, 4, 7]
		self.downloadtimelimit=10
		#self.maketimelimit=120
		self.maketimelimit=120
		self.testcasetimelimit=120
		self.downloadcheckinterval=3
		self.makecheckinterval=3
		self.testcasecheckinterval=5
		self.dougexe=None
		
		try:
        		opts, args = getopt.getopt(sys.argv[1:], "hdp:", ["help", "testplace=", "nogridexec", "nogridsave", "nocompilation", "makedeffile=", "notest", "nocodegetting", "make-timelimit=", "download-timelimit=", "testcase-timelimit=", "make-checkinterval=", "download-checkinterval=", "testcase-checkinterval=", "dougexe="])
		except getopt.GetoptError:
			# print help information and exit:
			usage()
			sys.exit(2)
		print "opts:",opts
		print "args:",args
		#kohalikud muutujad
		for o, a in opts:
			if o in ("-h", "--help"):
				usage()
				sys.exit()
			if o in ("-d"):
				self.DEBUG=True
			if o in ("--testplace"):
				if(a=="local" or a=="grid"):
					self.testplace=a #local voi grid, default local
			if o in ("--nocompilation"):
				self.nocompilation=True
			#if o in ("--dougloc"):
			#	print "siia ei peaks tulema"
			#	self.dougloc=a
			#if o in ("--examplesloc"):
			#	self.examplesloc=a
			if o in ("--makedeffile"):
				self.makedeffile=a #abs local path
			if o in ("--notest"):
				self.notest=True
			if o in ("--nocodegetting"):
				self.nocodegetting=True
			if o in ("-p"):
				self.numberofcopies=self.editoption(a)
			if o in ("--make-timelimit"):
				if(len(a)>0 and isdigit(a)==True and int(a)>=self.makecheckinterval):
					self.maketimelimit=int(a)
			if o in ("--testcase-timelimit"):
				if(len(a)>0 and isdigit(a)==True and int(a)>=self.testcasecheckinterval):
					self.testcasetimelimit=int(a)
			if o in ("--download-timelimit"):
				if(len(a)>0 and isdigit(a)==True and int(a)>=self.downloadcheckinterval):
					self.downloadtimelimit=int(a)
			if o in ("--make-checkinterval"):
				if(len(a)>0 and isdigit(a)==True and int(a)<=self.maketimelimit):
					self.makecheckinterval=int(a)
			if o in ("--download-checkinterval"):
				if(len(a)>0 and isdigit(a)==True and int(a)<=self.downloadtimelimit):
					self.downloadcheckinterval=int(a)
			if o in ("--testcase-checkinterval"):
				if(len(a)>0 and isdigit(a)==True and int(a)<=self.testcasetimelimit):
					self.testcasecheckinterval=int(a)
			if o in ("--dougexe"):
				if len(a)>0:
					self.dougexe=a
		
		
		#if self.makedefloc==None:
		#	self.makedefloc=os.path.join(self.dougdir, "trunk/src")
			
		#ifaceAB=interfaceAB.InterfaceAB()
		#errorh=errorhandler.Errorhandler(ifaceAB)
		self.repoi=repointerface.Repointerface(self.DEBUG)
		self.watchd=watchdog.Watchdog(self.repoi, self.testplace, self.DEBUG)
		
		#self.main()

	def main(self):
		print "numberofcopies=",self.numberofcopies
		print "dougtest main";
		self.curdirabspath=os.getcwd()
		self.dougtesthome=self.curdirabspath
		print "curdirabspath="+self.curdirabspath
		self.testdirabspath=self.dougtesthome+"/"+self.id
		os.mkdir(self.testdirabspath)
		self.makedefdirabspath=self.curdirabspath+"/"+self.makedefrelpath
		self.makefileabspath=self.curdirabspath+"/"+self.makefilerelpath
		self.makeoutfile=self.testdirabspath+"/dougtestmakeout.txt"
		self.makeerrfile=self.testdirabspath+"/dougtestmakeerr.txt"
		if self.dougexe is not None:
			#kontroll, kas path exists
			if os.path.exists(self.dougexe):
				self.dougexeabspath=self.dougexe
			else: self.dougexeabspath=self.curdirabspath+"/"+self.dougexerelpath+"/"+self.dougexefilename
		else: self.dougexeabspath=self.curdirabspath+"/"+self.dougexerelpath+"/"+self.dougexefilename
		
		
		#siin voiks options faili paigutada
		
		if self.nocodegetting==False:
			#download the code
			print "dougaddr=",self.dougaddr
			self.repoi.getcode(self.dougaddr)
		else: print "nocodegetting: koodi allalaadimist ei toimu"
		#kood saadud
		#eeldab, et kood kataloogis ./doug/trunk
		
		#ajut kopimine
		#command=["cp","--reply=yes", "/home/anti/proged/doug/trunk/src/Make.def", os.getcwd()+"/doug/trunk/src"]
		
		
		if self.nocompilation==False:
			#self.makedefdirabspath=self.curdirabspath+"/"+self.makedefrelpath
			if self.makedeffile is None: #rename Make.def.example->Make.def, ku juba ei eksisteeri
				#rename
				#if os.path.isfile("doug/trunk/src/Make.def.example"):
				
				if os.path.exists(self.makedefdirabspath+"/Make.def")==False:
					if  os.path.isfile(self.makedefdirabspath+"/Make.def.example"):
						p=popen2.Popen3("mv --reply=yes "+self.makedefdirabspath+"/Make.def.example "+self.makedefdirabspath+"/Make.def", capturestderr=1)
						#p.wait()
						i=0
						while p.poll()==-1:
							if(i>4):
								print "Viga: mv ei suutnud 10 sekundiga lopetada.Valjun"
								sys.exit(3)
							i=i+1
							time.sleep(i)
						mvexcode=p.poll()
						errlines=p.childerr.readlines()
						if mvexcode == 0:
							print "Make.def.example nime muutmise op edukas, excode=",mvexcode
							if(len(errlines)>0):
								print "mv kask andis hoiatusi:"
								for line in errlines:
									print "mv hoiatus:"+line
						else:
							print "mv Viga: Make.def.example nime muutmise op ebaonnestus, excode=",mvexcode
							for line in errlines:
								print "mv viga:"+line
							sys.exit(1)
						#errlines=p.childerr.readlines()
						
						p.childerr.close() #vb peaks kasutama so.close()
						p.fromchild.close()
						p.tochild.close()
						del p
					else:
						print "Struktuuri viga: Make.def faili pole tapsustatud"
						print "ja ei leidu ka asukohas "+self.makedefdirabspath
						print "Valjun."
						sys.exit(1)
				else: 
					print "info: Make.def fail eksisteerib."
					
					
			else: #--self.makedeffile is not none
				#p=popen2.Popen3("cp --reply=yes /home/anti/proged/doug/trunk/src/Make.def ./doug/trunk/src/", capturestderr=1)
				p=popen2.Popen3("cp --reply=yes "+self.makedeffile+" "+self.makedefdirabspath, capturestderr=1)
				#p=popen2.Popen3(command, capturestderr=1)
				i=0
				while p.poll()==-1:
					if(i>4):
						print "Viga: cp ei suutnud 10 sekundiga lopetada.Valjun"
						sys.exit(3)
					i=i+1
					time.sleep(i)
				cpexcode=p.poll()
				errlines=p.childerr.readlines()
				if cpexcode == 0:
					print "Make.def kopeerimise op edukas, excode=",cpexcode
					if(len(errlines)>0):
						print "cp kask andis hoiatusi:"
						for line in errlines:
							print "cp hoiatus:"+line
				else:
					print "cp Viga: Make.def kopeerimise op ebaonnestus, excode=",cpexcode
					for line in errlines:
						print "cp viga:"+line
					sys.exit(1)
				#errlines=p.childerr.readlines()
				
				p.childerr.close()
				p.fromchild.close()
				p.tochild.close()
				del p
			
			#shutil.copy2("/home/anti/proged/doug/trunk/src/Make.def", "/doug/trunk/src/")
			#print "kopeerimine onnestus"
			
			print "Peaks tulema kompileerimine"
			#kompileerimine
			#p1=popen2.Popen3("cd doug/trunk/", capturestderr=1)
			#ajutiselt no compilation
			
			#pp=popen2.Popen3("strace -o makeout.txt -ff -F -tt -v -s 128 make -C doug/trunk/ -f Makefile", capturestderr=1)
			comm="make -C "+self.makefileabspath+" -f Makefile 1>> "+self.makeoutfile+" 2>> "+self.makeerrfile
			print "kompileerimise(make) kask:",comm
			#pp=popen2.Popen3("make -C "+self.makefileloc+" -f Makefile 1>> dougtestmakeout.txt 2>> dougtestmakeerr.txt")
			pp=popen2.Popen3(comm)
			makeprid=pp.pid
			#command=["make", "-C", "doug/trunk/", "-f", "Makefile", "main"]
			#command=["doug/trunk/src/updatedeps.py", "-d"]
			
			#pp=popen2.Popen3(command, capturestderr=1)
			#print "enne waiti"
			#time.sleep(30)
			#vaike kontroll
			#i=0
			i=self.maketimelimit
			while pp.poll()==-1:
				if(i<0):
					print "Viga: make ei suutnud minutiga lopetada.Kill make ja valjun."
					os.kill(makeprid, signal.SIG_DFL)
					sys.exit(3)
				i=i-self.makecheckinterval
				print "."
				time.sleep(self.makecheckinterval)
			#excode=pp.wait()
			excode=pp.poll()
			fail=file(self.makeerrfile)
			errlines=fail.readlines()
			fail.close()
			
			del pp
			if excode!=0:
				print "Kompileerimine valjus veaga!"
				if len(errlines)>0:
					print "Veateade:"
					for line in errlines:
						print "komp errline:"+line
				sys.exit(1)
			else: 
				print "DOUG kompileerimise protsess valjus edukalt, excode=",excode
				if len(errlines)>0:
					print "Kompileerimine andmis hoiatusi:"
					for line in errlines:
						print "komp hoiatusrida:"+line
						
			
			#fail=file("dougtestmakeerr.txt")
			#errlines=fail.readlines()
			#fail.close()
			#if(len(errlines)>0):
			#	print "make Veateade:"
			#	for line in errlines:
			#		print "errline:"+line
			#	print "Valjun"
			#	sys.exit(2)
		
			#print "Kompileerimine edukalt lopetatud!"
			
		#DOUG_main KONTROLL
		
		if os.path.exists(self.dougexeabspath):
			print "DOUG exe kontroll: eksisteerib "+self.dougexeabspath
			if os.path.isfile(self.dougexeabspath):
				print "DOUG exe kontroll: "+self.dougexeabspath+" ON file."
			else:
				print "Ei leidnud DOUG exe faili: "+self.dougexeabspath+" ei ole fail! Valjun."
				sys.exit(1)
		else:
			print "DOUG exe kontroll: EI eksisteeri "+self.dougexeabspath+" teed! Valjun."
			sys.exit(1)
		
		#errlines=pp.childerr.readlines()
		#print "-------------Error lines:---------------"
		#for line in errlines:
		#	print "errline:"+line
		
		#outlines=pp.fromchild.readlines()
		#print "-------------Output lines:---------------"
		#for line in outlines:
		#	print "outline:"+line
		
		
		#excode=pp.wait()
		#print "LOPETATUD...excode:"
		#errlines=p.childerr.readline()
		#if(len(errlines)>0):
		#	print "Viga kompileerimisel. Exitcode=",excode
		#	print "Veateade", errlines
		#	print "Valjun"
		#	sys.exit(1)
		
		#makeargs=["make", "-C", "doug/trunk", "-f", "Makefire"]
		#errf=os.tmpfile()
		#serr=newfile("tmpfail")
		#print suka
		#p=subprocess.Popen(makeargs, stderr=serr)
		
		#returncode=p.wait()
		#print "returncode_", returncode
		#lines=serr.readlines()
		#if(len(lines)>0):
		#	print "Viga kompileerimisel."
		#	print lines
		
		#pp=popen2.Popen3("make -C doug/trunk/ -f Makefile")
		#print "childprocess pid:", pp.pid
		#for i in range(10):
		#	cod=pp.poll()
		#	if(cod==-1):
		#		print i,". kord: childprocess NOT exited."
		#	else:
		#		print i,". kord: childprocess DID exit, code=",cod
		#		sys.exit(0)
		#	time.sleep(10) #ootame 30s
		#print "make ei suutnud 300 sekundiga lopetada"
		
		
		#errlines=p.childerr.readlines()
		#if(len(errlines)>0):
		#	print "Viga kompileerimisel. Veateade:"
		#	print errlines
		#	print "Valjun."
		#	sys.exit(1)
		
		#print "Kompileerimine edukalt lopetatud."
		
		
		#kaimalykkamine, mpirun kasuga
		
		
		
		#peaks kasutama os.popen3(..) -e voi subprocess
		#pp.wait()
		#errlines=pp.childerr.readlines()
		#print "Errorlines:", errlines
		
		#if(len(errlines)>0):
		#	print "Viga kompileerimisel! Veateade:"
		#	print errlines
		#	print "Valjun."
		#	sys.exit(1)
		#print "Vigu ei esinenud. DOUG kompileeritud."
		
		
		#sakuri=popen2.Popen3("mpirun -np 1 n0 /home/anti/proged/doug/trunk/bin/DOUG_main -f /home/anti/proged/doug_examples/elemental/generated/UW0.5_4x4Q/omactrl.file")
		#for i in range(10):
		#	cod=sakuri.poll()
		#	if(cod==-1):
		#		print i,". kord: childprocess NOT exited."
		#	else:
		#		print i,". kord: childprocess DID exit, code=",cod
		#		sys.exit(0)
		#	time.sleep(10) #ootame 30s
		#print "make ei suutnud 300 sekundiga lopetada"
		
		
		#print "hakkan sulgema pipe-sid"
		#pp.tochild.close()
		#pp.fromchild.close()
		#pp.childerr.close()
		#print "pipe'd suletud, exit"
		#sys.exit(1)
		
		#kood tommatud, kompileeritud
		#nyyd testima!
		if self.notest==False:
			self.watchd.testit(examplesaddr=self.examplesaddr, examplesdir=self.examplesdir, dougexeabspath=self.dougexeabspath, testdirabspath=self.testdirabspath, dougtesthome=self.dougtesthome, numberofcopies=self.numberofcopies)
		else:
			print "notest: testimist ei toimu!"
		
	def getlastcomp(self, path):
		"Returns last component of path"
		if path is None: return None
		if path=="": return ""
		lastc=""
		if(path[-1]=="/"):
			path=path[:-1]
			
	def editoption(self, a):
		s=a.strip()
		if s==None or s=="": return []
		list=s.split(",")
		listofnumbers=[]
		for i in list:
			if len(i)!=0:
				if i.isdigit(): listofnumbers.append(int(i))
				else:
					abilist=i.split("-")
					eelmine=-1
					for j in abilist:
						if len(j)!=0:
							if(eelmine==-1):	
								eelmine=j
							else:
								eelmine=int(eelmine)
								j=int(j)
								if(eelmine<=j):
									listofnumbers=listofnumbers + range(eelmine, j+1)
								eelmine=-1
						
		return listofnumbers
		

print "ALGUS"
print "proov-__name__=", __name__
if __name__=="__main__":
	x=Dougtest()
	x.main()
print "LOPP"
#sys.exit(0)
