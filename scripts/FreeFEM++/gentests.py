import os
import shutil
import subprocess

def executeFreeFEMTest(test, N):
	FILE = open("test.conf", "w")
	FILE.write(str(N))
	FILE.write(" ")
	FILE.close()
	subprocess.call(["FreeFem++ " + test + ".edp"], shell=True)
	shutil.move("assembled_mtx.txt", test + "_" + str(N) + "_mtx.txt")
	shutil.move("assembled_sol.txt", test + "_" + str(N) + "_sol.txt")
	shutil.move("assembled_rhs.txt", test + "_" + str(N) + "_rhs.txt")
	os.unlink("test.conf")

executeFreeFEMTest("poisson1", 20)
executeFreeFEMTest("poisson1", 40)
executeFreeFEMTest("poisson1", 60)
executeFreeFEMTest("poisson1", 80)

executeFreeFEMTest("poisson2", 20)
executeFreeFEMTest("poisson2", 40)
executeFreeFEMTest("poisson2", 60)
executeFreeFEMTest("poisson2", 80)

executeFreeFEMTest("poisson3", 20)
executeFreeFEMTest("poisson3", 40)
executeFreeFEMTest("poisson3", 60)
executeFreeFEMTest("poisson3", 80)
