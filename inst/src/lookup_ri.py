
from bs4 import BeautifulSoup as bs
import numpy as np
import sys
from urllib import request as req
import re 

base_url = 'https://webbook.nist.gov/cgi/cbook.cgi?ID='
gas_base_url = 'https://webbook.nist.gov'

def contain_gas_chrom(tag):
  '''
  custom search gas chromatography link for beautiful soup
  '''
  if tag.name == 'a' and '#Gas-Chrom' in tag.attrs['href']:
    return True
  else:
    return False

def get_gasChrom_url(cas_number):
  url = base_url + cas_number
  # print (url)
  page = req.urlopen(url)
  soup = bs(page, 'html.parser')
  try:
    gas_chrom_href = soup.find(contain_gas_chrom).attrs['href']
  except:
    gas_chrom_href = None
  # print (gas_chrom_href)

  return gas_chrom_href

def get_ri(cas_number, column_type):
  '''
    Get average RI from nistwebbook
    cas_number:  CAS number of the compound
    comlumn_type: '5ms' or 'wax'    
  '''
  try:
    gas_url = gas_base_url + get_gasChrom_url(cas_number)
    # print (gas_url)
    page = req.urlopen(gas_url)
  except:
    return -1

  gas_soup = bs(page, 'html.parser')
  exp_rows = gas_soup.find_all('tr', class_='exp')
  # return exp_rows
  RIs = []
  if column_type == '5ms':
    reg = '(-5)[^0-9]'
  elif column_type == 'wax':
    reg = '(wax)'
  for row in exp_rows:
    tds = row.find_all('td')
    column =tds[1].text
    match = re.search(reg, column, re.IGNORECASE)
    if not match:
      continue
    # print (column)
    # print (tds)
    ri = float(tds[2].text)
    if ri > 600:
      RIs.append(ri)

  RI = -1
  if len(RIs) == 0:
    return -1
  else:
    return np.mean(RIs)



def main():
    '''
    Returns average retention index from nist webbook for 5ms or wax column
    usage:  lookup_ri casnumber 5ms or lookup_ri casnumber wax
    '''
    if len(sys.argv) != 3:
        print ('Usage: lookup_ri casNumber 5ms or lookup_ri casNumber wax')
        return 0

    cas = sys.argv[1]
    column = sys.argv[2]
    # print (cas, column)
    ri = get_ri(cas, column)
    print (ri)
    return 


if __name__ == '__main__': 
    main()
