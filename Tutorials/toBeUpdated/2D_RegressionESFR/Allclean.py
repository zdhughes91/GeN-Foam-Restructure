import shutil
import sys
import os

whiteList = ["0", "0.orig", "system", "constant", "dependencies", "referenceResults", "Allclean.py", "Allrun.py", "runDict", "readme", "Allrun.sh", "Allclean.sh", "theory.pdf"]

for i in os.listdir(os.getcwd()) :
	if i not in whiteList :
		try :
			os.remove(i)
		except :
			try :
				shutil.rmtree(i) 
			except :
				pass

if "0.orig" in os.listdir(os.getcwd()) :
	try :
		shutil.rmtree("0")
	except :
		pass
	shutil.copytree("0.orig", "0")
	shutil.rmtree("0.orig")