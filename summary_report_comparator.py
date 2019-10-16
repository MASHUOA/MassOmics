from tkFileDialog import askopenfilename
from Tkinter import Tk
import pandas as pd
import tkSimpleDialog as simpledialog
from Tkinter import *
import glob, os


# get the first file
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
#inFile = askopenfilename(filetypes=[('1st summary report xlsx', '.xlsx')],
    #title='select a summary report to compare') # show an "Open" dialog box and return the path to the selected file
inFile = 'E:/DataProcessing/2019/Sarah Mitchell/CLIMB_SPME_GCMS_2/CLIMB_unknowns_fullset - WORKING__filterv08_summaryReport.xlsx'
sourceFolder = os.path.dirname(inFile)  # strip out the folder path
print('inFile1: ', inFile)
print('sourceFolder: ' + sourceFolder)
outFile = sourceFolder + '/summary_reporter_comparator.xlsx'
print('outFile: ', outFile)
print(' ')

# list the sheets in the chosen file
dd = pd.ExcelFile(inFile)
shits = dd.sheet_names
#shitNumber = 0
#for i in shits:
    #print('sheet number ' + str(shitNumber) + ': ' + i)
    #shitNumber +=1
## choose a sheet to load
#print(' ')

#shit = simpledialog.askinteger("select number of sheet to load", "select number from list")
shit = 2
df1 = pd.read_excel(inFile, sheet_name=shits[shit])
#df1.set_index('Name')
print(list(df1))


# get the second file
#inFile = askopenfilename(filetypes=[('1st summary report xlsx', '.xlsx')],
    #title='select a summary report to compare') # show an "Open" dialog box and return the path to the selected file
inFile = 'E:/DataProcessing/2019/Sarah Mitchell/CLIMB_SPME_GCMS_2/massOmics/Summary report_NIST v03 - Sarah good_Ref_Ions.xlsx'

print('inFile2: ', inFile)
print(' ')

# list the sheets in the chosen file
dd = pd.ExcelFile(inFile)
shits = dd.sheet_names
#shitNumber = 0
#for i in shits:
    #print('sheet number ' + str(shitNumber) + ': ' + i)
    #shitNumber +=1
#shit = simpledialog.askinteger("select number of sheet to load", "select number from list")
shit = 0
#print(' ')

df2 = pd.read_excel(inFile, sheet_name=shits[shit])
#df2.set_index('Name')
print(list(df2))
print(' ')


## combine the compound names and create a df
#compounds = list(df1['Name'])
##print('compounds: ', compounds)
##print(' ')
#print('compunds length: ', len(compounds))
##print(' ')
##print(' ')

#compounds = compounds + list(df2['Name'])
##print('compounds: ', compounds)
##print(' ')
#print('compunds length: ', len(compounds))
##print(' ')
##print(' ')

#compounds = list(set(compounds))
##print('compounds: ', compounds)
##print(' ')
#print('compunds length: ', len(compounds))
##print(' ')
##print(' ')

#for i in compounds:
	#compound1 = df1[df1['Name'] == i]
	#compound2 = df2[df2['Name'] == i]
	#compound_row = pd.concat([compound1, compound2], axis=1)
	#print(compound_row)
	#quit()


# pandas is awesome
combi = df1.merge(df2, on='Name', how='outer')


# Create a Pandas Excel writer using XlsxWriter as the engine
writer = pd.ExcelWriter(outFile, engine='xlsxwriter')
combi.to_excel(writer, sheet_name='awesome')
# Close the Pandas Excel writer and output the Excel file.
writer.save()

print(' ')
print('awesomeness accomplished')