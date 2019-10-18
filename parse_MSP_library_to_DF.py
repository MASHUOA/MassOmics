# pass this a list of MSP files and it will return the data parsed to a DF

import pandas as pd
import glob, os
import platform

# test data
#f = 'E:/DataProcessing/2019/Sarah Mitchell/CLIMB_SPME_GCMS_2/massOmics/CLIMB_unknowns_fullset - WORKING_badhitschecked__filtered_v05-correctedCAS.csv'
# f = '/media/chris/data/Chris_stuff/Google_Drive/KODE/my_Python_code/mass_spec/MSL_parse/massOmics/test/test.csv'
#f = 'E:/Chris Pook/code/MSL_library_writer/massOmics/CLIMB_unknowns_fullset - WORKING_badhitschecked__filtered_v05-correctedCAS.csv'


def MSP_parse(fileList, debug):

    if platform.system() == 'Windows':
        nuRow = '\n'
    else:
        nuRow = '\r\n'

    MSL_df = pd.DataFrame(columns = [])
    df_row = 0

    for i in fileList:
        stringDict = {}
        spectrumDict = {}
        names_rows = []
        print(' ')
        print('~ *~ *~ *~ *~ *~ * ~* ~* ~* ~* ~* ~* ~* ~* ~*')
        mspLib = open(i, 'r')
        lines = mspLib.readlines()
        mzList = []
        intensityList = []
        finale = 0
        basepeak_mz = 0
        for e in lines:
             if ":" in e:
                 finale = 0
                 splitString = e.split(':')
                 if len(splitString) >2:
                     for a in range(2,len(splitString)):
                         splitString[1] = splitString[1] + splitString[a]
                 stringDict[splitString[0]] = splitString[1].strip()
             elif e ==nuRow:
                 finale +=1
                 if finale >1:
                    break
                 if debug: print(stringDict['Name'])
                 stringDict['mzList'] = mzList
                 stringDict['intensityList'] = intensityList
                 stringDict['Source'] = i
                 stringDict['basepeak'] = basepeak_mz
                 if debug: print( stringDict)
                 MSL_df = MSL_df.append(stringDict, ignore_index=True)
                 mzList = []
                 intensityList = []
                 stringDict = {}
             else:
                xx = e.strip().split(' ')
                mzList.append(xx[0].strip())
                intensityList.append(xx[1].strip())
                if xx[1] == '999.00':
                    basepeak_mz = xx[0]

    return MSL_df


