# Python 2
# point this at any one CEF file in a folder full
# it will compile and return [name, n, RT, ref, area, height, CAS] from each one
# these are collected in a pandas df, CEF_df
def parse_cef(Folder = None ):
	import glob, os
	from tkinter import Tk
	from tkinter.filedialog import askopenfilename
	from tkinter.filedialog import askdirectory
	import numpy
	import platform
	import pandas as pd
	debug = True
	
	# point this at a CEF file
	# function returns (name, n, RT, ref, area, height, CAS)
	def parse_CEF(cef):
		n = 0
		ref = []
		RT = []
		area = []
		height = []
		name = []
		CAS = []
		f = open(cef, 'r')
		lines = f.readlines()
		for i in lines:
			if "      <Location " in i:
				x = i.split(' ')
				mm = x[7].split('=')
				mmm = mm[1].strip('\"')
				ref.insert(n,float(mmm))
				rt = x[8].split('=')
				rtt = rt[1].strip('\"')
				RT.insert(n, float(rtt))
				a = x[9].split('=')
				aa = a[1].strip('\"')
				area.insert(n, float(aa))
				y = x[9].split('=')
				yy = y[1].strip('\"')
				height.insert(n, float(yy))
			elif '        <Molecule name=' in i:
				nom = i.split('=')
				nomm = nom[1].split('\"')
				name.insert(n, nomm[1])
			elif '            <Accession db=' in i:
				c = i.split('id=')
				cc = c[1].split('\"')
				CAS.insert(n, cc[1])
			elif '</Compound>' in i:
				n = n +1

		return(name, n, RT, ref, area, height, CAS)

	#cef = '/media/chris/data/Chris_stuff/Google_Drive/KODE/my_Python_code/mass_spec/Agilent unknowns proc/drive-download-20170717T093726Z-001/20170621_MCFvXX3_CJP_250uM_AA_STD-AllHits.cef'
	#print parse_CEF(cef)

	# select any CEF file in your folder of results
	Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
	
	if Folder is None: Folder = askdirectory(title='Please select a directory contains all cef files')  # strip out the folder path
	os.chdir(Folder)
	files = []
	for file in glob.glob("*.cef"):  # compile a list of all cef files in the folder
		if debug: print(file)
		files.append(file)


	CEF_df = pd.DataFrame()

	for i in files:  # parse each file and add to the list of lists
		if debug: print('processing file ' + str(i))
		n = i.split('.')
		dat = parse_CEF(Folder + '/' + i)
		# function returns [name, n, RT, ref, area, height, CAS]
		compounds = {'name':dat[0]}
		compounds['n'] = dat[1]
		compounds['RT'] = dat[2]
		compounds['ref'] = dat[3]
		compounds['area'] = dat[4]
		compounds['height'] = dat[5]
		compounds['CAS'] = dat[6]
		CEF_df = CEF_df.append(pd.DataFrame.from_dict(compounds))


	if debug:
		print('processing complete')
		return((CEF_df))
		#names = list(CEF_df['name'])
		#for i in names:
		#	print(i)