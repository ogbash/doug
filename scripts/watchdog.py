#! /usr/bin/env python
import sys, os, popen2, time, signal


class Watchdog:
	
	def log(self, s):
		if(self.DEBUG==True):
			print s
	
	def __init__(self, repoi, testplace, debug=False, logdest="stdout"):
		self.DEBUG=debug
		self.log("Watchdog init()")
		self.logdest=logdest
		self.testplace=testplace
		self.repoi=repoi
		self.log("Watchdog init lopp")
		
	
	def testit(self, examplesaddr, examplesdir, dougexeabspath, testdirabspath, dougtesthome, numberofcopies):
		self.log("Watchdog.testit() algus")
		self.examplesaddr=examplesaddr
		self.examplesdir=examplesdir
		self.dougexeabspath=dougexeabspath
		self.testdirabspath=testdirabspath
		self.dougtesthome=dougtesthome
		self.numberofcopies=numberofcopies
		
		#print "proov.__name__=", __name__
		if self.testplace=="local":
			self.localtest()
		elif selftestplace=="grid":
			self.gridtest()
		
		self.log("Watchdog.testit() lopp")
		
	def gridtest(self):
		pass
			
	def localtest(self):
		self.log("Watchdog.localtest() algus")
		
		#tomba ex
		print "examplesaddr=",self.examplesaddr
		self.repoi.getcode(self.examplesaddr)
		
		#iga testi kohta: 1, 3, 6 protsessoriga(-np6)
		#kus testid on
		#kontrollfaili asukoht teada
		self.ctrlfiledirrelpath="elemental/generated" #in examples dir
		self.ctrlfilename="ctrl.file"
		self.ctrlfileabspath=self.dougtesthome+"/"+self.examplesdir+"/"+self.ctrlfiledirrelpath+"/"+self.ctrlfilename
		self.testdatadirsrelpath="elemental/generated"
		ctr="doug_examples/elemental/generated/ctrl.file"
		logoutfile="log.out"
		logerrfile="log.err"
		
		#pyyab testijuhtumite kaustad saada
		self.testdatadirsabspath=self.dougtesthome+"/"+self.examplesdir+"/"+self.testdatadirsrelpath
		print "testdatadirsabspath="+self.testdatadirsabspath
		testdirslistrel=os.listdir(self.testdatadirsabspath)
		testdirslistabs=[]
		for dir in testdirslistrel:
			print "dir:",dir
			if dir.startswith(".")==False:
				ajut=self.testdatadirsabspath+"/"+dir
				if os.path.isdir(ajut):
					testdirslistabs.append(ajut)
				
		print "testdata dirlist with abs paths:",testdirslistabs
		#ctrlfile1="/home/anti/proged/doug_examples/assembled/Lap.ctl"; #assembled
		#testdatapath1="/home/anti/proged/doug_examples/assembled/"; #andmefailid koik koos
		#ctrlfile2="/home/anti/proged/doug_examples/elemental/generated/ctrl.file"
		#testdatapath2="/home/anti/proged/doug_examples/elemental/degenerate/"
		#ctrlfile3="/home/anti/proged/doug_examples/elemental/generated/ctrl.file"
		#testdatapath3="/home/anti/proged/doug_examples/elemental/generated/e4x4/"
		#ctrlfile4="/home/anti/proged/doug_examples/elemental/generated/ctrl.file"
		#testdatapath4="/home/anti/proged/doug_examples/elemental/generated/h4x4/"
		#dougasub="/home/anti/ope/baka/dougtest/doug/trunk/bin/DOUG_main"
		
		#paneb jadamisi toole
		
		#ks lam-mpi voiks olla installeeritud...kontrollib lamboot kasku, ega viga ei teki
		
		
		#kas lam-mpi tootab
		#lamboot log files:+abs paths to doughome
		lambooterr=self.testdirabspath+"/lamboot.err"
		lambootout=self.testdirabspath+"/lamboot.out"
		
		p=popen2.Popen3("lamboot 1>"+lambootout+" 2>"+lambooterr)
		lambootprid=p.pid
		i=0
		while p.poll()==-1:
			if(i>9):
				print "Viga: lamboot kask ei suutnud 30 sekundiga lopetada."
				print "Ei onnestunud kaivitada lam-mpi.Valjun"
				os.kill(lambootprid, signal.SIG_DFL)
				sys.exit(3)
			i=i+1
			print "."
			time.sleep(3)
		lambootexcode=p.poll()
		if lambootexcode!=0:
			print "Ei suutnud lami kaivitada kasuga lamboot, tagastas vea, excode=",lambootexcode
			print "Valjun."
			sys.exit(3)
		else: print "kask lamboot valjus normaalselt, lam keskkond yleval, excode",lambootexcode
		del p
		
		#laminfo testikausta 
		lamp=popen2.Popen3("laminfo 1>"+self.testdirabspath+"/laminfo.out 2>"+self.testdirabspath+"/laminfo.err")
		laminfoprid=lamp.pid
		i=0
		while lamp.poll()==-1:
			if(i>9):
				print "Viga: laminfo ei suutnud 30 sekundiga lopetada."
				os.kill(laminfoprid, signal.SIG_DFL)
				sys.exit(3)
			i=i+1
			time.sleep(3)
		
		laminfoexcode=lamp.poll()
		if laminfoexcode!=0:
			print "Hoiatus: laminfo abnormaalne valjumine, laminfoexcode=",laminfoexcode
		else: print "laminfo tagastus normaalne, excode=",laminfoexcode
		
		
		
		#TESTIMINE
		#ctrl faili kontroll
		#testijuhtumid: iga ctrl faili kohta testida 1, 4, 7(default) protsessoriga, aegu salv pole vaja
		if os.path.isfile(self.ctrlfileabspath):
			print "kontroll: ctrl file olemas."
		else:
			print "kontroll: ctrl faili ei leitud. Pole testida millegi pohjal.valjun."
			sys.exit(3)
		
		if len(testdirslistabs)==0:
			print "kontroll: Yhtegi kataloogi testandmefailidega ei leitud. Pole millegi pohjal testida.valjun."
			sys.exit(3)
		else:
			print "kontroll: testandmete kaustasid vahemalt 1."
			
		testiloendur=1
		for dir in testdirslistabs:
			print "------- Testcase "+str(testiloendur)+" -------"
			thistestcasedirabspath=self.testdirabspath+"/testcase"+str(testiloendur)
			os.mkdir(thistestcasedirabspath)
			infofileabspath=thistestcasedirabspath+"/testcase"+str(testiloendur)+".txt"
			infofile=open(infofileabspath, "a")
			infofile.writelines(["testdatafilesdirabspath="+dir+"\n", "ctrlfileabspath="+self.ctrlfileabspath+"\n"])
			infofile.close()
			lugeja=1
			for iter in self.numberofcopies:
				#kataloog sama testi erinevate protsesside arvuga jooksutamisel
				protsessidearv=lugeja
				print "testcase "+str(testiloendur)+" protsesside arvuga "+str(protsessidearv)
				absprocdir=thistestcasedirabspath+"/procs"+str(protsessidearv)
				os.mkdir(absprocdir)
				os.chdir(dir)
				command="mpirun -np "+str(protsessidearv)+" "+self.dougexeabspath+" -f "+self.ctrlfileabspath+" -q 1> "+logoutfile+" 2> "+logerrfile
				print "mpirun command: "+command
				mpi=popen2.Popen3(command)
				processid=mpi.pid
				i=0
				while mpi.poll()==-1:
					time.sleep(10)
					if(i>9):
						print "Viga: mpirun kask programmiga DOUG ei suutnud 90 sekundiga lopetada."
						print "Testijuhtum "+str(testiloendur)+" protsesside arvuga "+str(protsessidearv)+" ebaonnestus.Kill process id ",processid
						os.kill(processid, signal.SIG_DFL)
						#sys.exit(3)
					i=i+1
					print "."
					
				excode=mpi.poll()
				if excode!=0:
					print "Testijuhtum "+str(testiloendur)+" protsesside arvuga "+str(protsessidearv)+" on ebaonnestunud.excode=",excode
				else: print "Testijuhtum "+str(testiloendur)+" protsesside arvuga "+str(protsessidearv)+" valjus edukalt.excode=",excode
				
				#kopeeri failid
				#algul yldised logifailid
				#p.popen2.Popen3("mv -f "+logoutfile+" "+absprocdir)
				try:
					os.rename(logoutfile, absprocdir+"/"+logoutfile) #voiks lisada isfile checkid
					os.rename(logerrfile, absprocdir+"/"+logerrfile)
					for kord in range(protsessidearv):
						os.rename("log."+str(kord), absprocdir+"/log."+str(kord))
				except OSError:
					print "OS viga voi kask os.rename(2) 2. arg on dir."
					
				lugeja=lugeja+3
			
			print "Testcase "+str(testiloendur)+" lopetatud."
			testiloendur=testiloendur+1
				
				
		#os.chdir(dir)
		#command="mpirun -np 3 "+self.curdir+"/doug/trunk/bin/DOUG_main -f "+self.curdir+"/"+ctr+" -q 1> log.out 2> log.err"
		#os.chdir(testdatapath1)
		#os.chdir(testdatapath3)
		#comm="mpirun -np 1 n0 "+dougasub+" -f "+ctrlfile2+" -q 1> mlog.0 2> testcase1.err"
		#comm="mpirun -np 3 "+self.curdir+"/doug/trunk/bin/DOUG_main -f "+self.curdir+"/"+ctr+" -q 1> log.out 2> log.err"
		#print "KASK:"+comm
		#mpi=popen2.Popen3(comm, capturestderr=1)
		#i=0
		#while mpi.poll()==-1:
		#	if(i>9):
		#		print "Viga: mpirun kask programmiga DOUG ei suutnud 30 sekundiga lopetada."
		#		print "Testijuhtum ebaonnestus."
		#		sys.exit(3)
		#	i=i+1
		#	time.sleep(10)
		#errlines=mpi.childerr.readlines()
		#if len(errlines)>0:
		#	print "mpirun veateade:"
		#	for line in errlines:
		#		print line
		#	print "ei suutnud lopetada runmpi kasku.valjun."
		#	sys.exit(3)
		print "mpirun lopetanud"
		os.chdir(self.dougtesthome)
		#dr=self.testdir+"/"+"test"+str(i)
		#print "dr:",dr
		#os.mkdir(dr)
		#ppp=popen2.Popen3("mv -f log* "+dr)
		#cood=ppp.wait()
		#if cood==-1:
		#	print "logide liigutamine ebaonnestus"
		#else: print "logid liigutatud"
	
		#TESTIMINE
		
		self.log("Watchdog.localtest() lopp")
		
		
	def gridtest(self):
		pass
		
	def laminfokontroll(self):
		pass
		
	def mpirunkontroll(self):
		pass