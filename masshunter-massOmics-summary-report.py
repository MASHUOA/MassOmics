# feed this your WORKING xlsx, your original output xlsx and one of your MSP sub-libraries

# anything with an incidence less than the threshold will be ditched
# secondary filter uses SNR
# tertiary QC information regarding RT RSD is compiled and presented to help identify
# good hits that you can accept and bad hits that need manual scrutiny and curation
# MSP libraries will be compiled into a DF

import pandas as pd
from collections import Counter
from tkinter import Tk
from tkinter.filedialog import askopenfilename
from parse_MSP_library_to_DF import MSP_parse
import glob, os
from fuzzywuzzy import process


# SET THE THRESHOLD FOR INCIDENCES OF EACH COMPOUND, BELOW WHICH TO EXCLUDE THAT TARGET
incidenceThreshold = 2

# SNR threshold
SNRthreshold = 10

# do you want to export a summary for downstream processing with massOmics?
massOmics_summaryReport = True


width = 15
RT_shift = 0.5
debug = False

# get your excel file
Tk().withdraw() # we don't want full GUI, so keep the root window from appearing
f = askopenfilename(filetypes=[('excel file', '.xlsx')],
    title='select your WORKING xlsx results file') # show an "Open" dialog box and return the path to the selected file

original_results_xlsx = askopenfilename(filetypes=[('excel file', '.xlsx')],
    title='select your original xlsx results file') # show an "Open" dialog box and return the path to the selected file
orig_df = pd.read_excel(original_results_xlsx)
orig_names_set = set(orig_df['Compound Name'])

MSP_file = askopenfilename(filetypes=[('MSP library file', '.MSP')],
    title='select the MSP library file') # show an "Open" dialog box and return the path to the selected file
Folder = os.path.dirname(MSP_file)  # strip out the folder path
#  print('Folder: ', Folder)
os.chdir(Folder)
MSP_file_list = []
#  print('MSP files to parse: ')
for file in glob.glob("*.MSP"):  # compile a list of all msp files in the folder CASE SENSITIVE
    MSP_file_list.append(file)
    #  print(str(file))

# parse the MSP to a df
MSP = MSP_parse(MSP_file_list, False)
MSP_names = set(list(MSP['Name']))
# if debug:  print('length of MSP df: ' + str(len(MSP)))


aa = f.split('.')
outFile  = aa[0] + "__filterv08_summaryReport.xlsx"

d = pd.read_excel(f)
columnNames = list(d)

# filter out any features with SNR less than the threshold and those from blanks
dFiltered = d[(d['Base Peak SNR'] >= SNRthreshold) & (d['Type'] != 'Blank')]

# RT threshold for RSD
threshold = 0.02    # 2%

# put all the names of compounds you're starting with into a list
allNames = list(set(dFiltered['Compound Name']))

# list to hold the filtered compound data
filtered_compound_data = []

# use a Counter to count incidences of each name
# https://stackoverflow.com/questions/2600191/how-to-count-the-occurrences-of-a-list-item/23909767
nn = Counter(dFiltered['Compound Name'])
n, m = nn.keys(), nn.values()

# iterate through the compounds and compile the ones that exceed the criteria to a list
for i in range(len(n)):
    if m[i] >= incidenceThreshold:
         #  if debug:print('n[i]: ' + n[i])
        subset = dFiltered[dFiltered['Compound Name'] == n[i]]
        # rationalise the compound name against the original xlsx output & extract a sub with the basepeak
        real_name_tuple= process.extract(n[i], orig_names_set, limit=1)[0]
        #  if debug: print('real_name_tuple: ' + str(real_name_tuple[0]) + ', ' + str(real_name_tuple[1]))
        basepeak_sub = MSP[MSP['Name'] == real_name_tuple[0]]
        CAS_sub = orig_df[orig_df['Compound Name'] == real_name_tuple[0]]
        fuzzy_basepeak_mz = 'name_not_matched'
        fuzzy_CAS = 'name_not_matched'
        if real_name_tuple[1] > 92:
            try: fuzzy_basepeak_mz = max(list(set(basepeak_sub['basepeak'])))
            except ValueError:  fuzzy_basepeak_mz = str(int(max(list(set(CAS_sub['Base Peak MZ'])))))
            try: fuzzy_CAS = set(CAS_sub['CAS#']).pop()
            except ValueError:  fuzzy_CAS = 'ValueError'
        #  if debug:
            #  print('fuzzy_basepeak_mz: ' + str(fuzzy_basepeak_mz))
            #  print('length of basepeak_sub: ', len(basepeak_sub))
            #  print('fuzzy_CAS: ' + str(fuzzy_CAS))
            #  print('length of CAS_sub: ', len(CAS_sub))
        # if columnNames contains 'Best Hit' then summarise that data
        if columnNames.count('Best Hit') >0:    # calculate nHits & nBestHits
            hitList = list(subset['Best Hit'])
            nBestHits = hitList.count(True)
            bestHitRatio = float(nBestHits) / float(m[i])    # ratio of best hits to incidence
        else:    # use filler data
            nBestHits = 'no_best_hit'
            bestHitRatio = 'no_best_hit'

        data = [n[i],    # compile the compound data to a list
            real_name_tuple[0],
            fuzzy_CAS,
            m[i],
            subset['Component RT'].mean(),
            subset['Component RT'].median(),
            subset['Component RT'].std(),
            subset['Component RT'].std() / subset['Component RT'].mean(),
            subset['Base Peak MZ'].mean(),
            subset['Base Peak MZ'].std()/subset['Base Peak MZ'].mean(),
            subset['Match Factor'].mean(),
            subset['Match Factor'].std()/subset['Match Factor'].mean(),
            nBestHits,
            bestHitRatio,
            fuzzy_basepeak_mz]
        filtered_compound_data.append(data)


# output some summary stats
# #  print 'length of allNames: ' + str(len(allNames))
# #  print 'length of filteredNames: ' + str(len(filtered_compound_data))

dfColumns = ['WORKING_name',
            'fuzzy_name',
            'fuzzy_CAS',
            'sample_incidences',
            'RT_mean',
            'RT_median',
            'RT_SD',
            'RT_RSD',
            'mz_mean',
            'mz_RSD',
            'MF_mean',
            'MF_RSD',
            'nBestHits',
            'bestHitsRatio',
            'fuzzy_basepeak_mz']

filtered_compound_data_df = pd.DataFrame.from_records(filtered_compound_data, columns = dfColumns)
good_df = filtered_compound_data_df[filtered_compound_data_df['RT_RSD'] <= threshold]
bad_df = filtered_compound_data_df[filtered_compound_data_df['RT_RSD'] > threshold]
#  #  print 'number of goodHits: ' + str(len(good_df['WORKING_name']))
#  print 'number of badHits: ' + str(len(bad_df['WORKING_name']))


# write the output df to excel
writer = pd.ExcelWriter(outFile, engine='xlsxwriter') # Create a Pandas Excel writer using XlsxWriter as the engine
good_df.to_excel(writer, sheet_name='GC-MS hit filter v08 GOOD HITS')
bad_df.to_excel(writer, sheet_name='BAD HITS')


# ~ ~ ~   ~ ~ ~   ~ ~ ~   ~ ~ ~   ~ ~ ~   ~ ~ ~
if massOmics_summaryReport:
    CAS_list = []    # list for the CAS numbers
    name_list = list(filtered_compound_data_df['WORKING_name'])
    for i in name_list:
        subset = orig_df[orig_df['Compound Name'] == i]
        # if debug:  print(len(subset))
        CAS = list(subset['CAS#'])
        try:
            #  print('CAS[0]: ' + str(CAS[0]))
            CAS_list.append(str(CAS[0]))
        except IndexError:
            #  print('CAS for ' + i + ' missing')
            CAS_list.append('missing')
    filtered_compound_data_df['CAS'] = CAS_list

    peak_width_list = []
    for ie in range(0,len(CAS_list)):
        peak_width_list.append(width)
    # if debug:  print('peak_width_list length: ' + str(len(peak_width_list)))
    filtered_compound_data_df['width'] = peak_width_list

    massOmics_summary_report_columns = ['Name',
                                        'Ref.ion',
                                        'Library.match',
                                        'Total.ID',
                                        'RT.median',
                                        'RT.shfit.Lower',
                                        'RT.shfit.upper',
                                        'Peak.width',
                                        'CAS']

    MO_data_dict = {'Name': filtered_compound_data_df['WORKING_name'],
                    'Ref.ion':filtered_compound_data_df['fuzzy_basepeak_mz'],
                    'Library.match':filtered_compound_data_df['MF_mean'],
                    'Total.ID':filtered_compound_data_df['sample_incidences'],
                    'RT.median':filtered_compound_data_df['RT_median'],
                    'RT.shfit.Lower':filtered_compound_data_df['RT_median'] -RT_shift,
                    'RT.shfit.upper':filtered_compound_data_df['RT_median'] +RT_shift,
                    'Peak.width':filtered_compound_data_df['width'] +RT_shift,
                    'CAS':filtered_compound_data_df['fuzzy_CAS']}
    #massOmics_df = pd.DataFrame(MO_data_dict, columns= [massOmics_summary_report_columns])
    massOmics_df = pd.DataFrame(MO_data_dict)
    # if debug:
    	#  print('MO_data_dict[Name]: ', MO_data_dict['Name'][:3])
    	#  print('massOmics_df header: ' ,massOmics_df.head)
    	#  print('massOmics_df length: ' , massOmics_df.head)

    massOmics_df.to_excel(writer, sheet_name='Summary_report_NIST')


# Close the Pandas Excel writer and output the Excel file.
writer.save()

#  print('script completed!')
