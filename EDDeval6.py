GCTLs = {'benzene':1,
         'acetone':6300,
         'acenaphthene':20,
         'acenaphthylene':210,
         'anthracene':2100,
         'benzo(a)anthracene':0.05,
         'benzo(a)pyrene':0.2,
         'benzo(b)fluoranthene':0.05,
         'benzo(g,h,i)perylene':210,
         'benzo(k)fluoranthene': 0.5,
         'chrysene':4.8,
         'dibenz(a,h)anthracene':0.005,
         'dibromochloromethane':0.4,
         'fluoranthene':280,
         'fluorene':280,
         'indeno(1,2,3-cd)pyrene':0.05,
         'isopropylbenzene':0.8,
         'cumene':0.8,
         '1-methylnaphthalene':28,
         '1-methylnaphthalene':28,
         '2-methylnaphthalene':28,
         'naphthalene':14,
         'phenanthrene':210,
         'pyrene':210,
         'ethylbenzene':30,
         'toluene':40,
         'xylenes, total':20,
         'xylenes- total':20,
         'dibromoethane, 1,2-':0.02,
         'edb':0.02,
         'dichloroethane, 1,2-':3,
         'mtbe':20,
         'methyl-t-butyl ether':20,
         'trphs':5000,
         'total recoverable pet. hydrocarbons':5000,
         'arsenic':10,
         'cadmium':5,
         'chromium (total)':100,
         'lead':15}

NADCs = {'benzene':100,
         'acenaphthene':200,
         'acenaphthylene':2100,
         'anthracene':21000,
         'benzo(a)anthracene':5, 
         'benzo(a)pyrene':20,
         'benzo(b)fluoranthene':5,
         'benzo(g,h,i)perylene':2100,
         'benzo(k)fluoranthene': 50, 
         'chrysene':480,
         'dibenz(a,h)anthracene':0.5,
         'fluoranthene':2800,
         'fluorene':2800,
         'indeno(1,2,3-cd)pyrene':5,
         'isopropylbenzene':8,
         'cumene':8,
         '1-methylnaphthalene':280,
         '1-methylnaphthalene':280,
         '2-methylnaphthalene':280,
         'naphthalene':140,
         'phenanthrene':2100,
         'pyrene':2100, 
         'ethylbenzene':300,
         'toluene':400,
         'xylenes, total':200,
         'xylenes- total':200,
         'dibromoethane, 1,2-':2,
         'edb':2,
         'dichloroethane, 1,2-':300, 
         'mtbe':200,
         'methyl-t-butyl ether':200,
         'trphs':50000,
         'total recoverable pet. hydrocarbons':50000,
         'arsenic':100,
         'cadmium':50,
         'chromium (total)':1000,
         'lead':150}

SCTLs = {'acenapthene':2400000,
         'acenaphthylene':1800000,
         'anthracene':21000000,
         'benzo(b)anthracene': None,
         'benzo(a)pyrene':100,
         'benzo(b)fluoranthene':None,
         'benzo(g,h,i)perylene':2500000,
         'benzo(k)fluoranthene':None,
         'chrysene':None,
         'dibenzo(a,h)anthracene':None,
         'fluoranthene':3200000,
         'fluorene':2600000,
         'indeno(1,2,3-cd)pyrene':None,
         '1-methylnaphthalene':200000,
         '2-methylnaphthalene':210000,
         'naphthalene':55000,
         'phenanthrene':2200000,
         'pyrene':2400000,
         'benzene':1200,
         'ethylbenzene':1500000,
         'toluene':7500000,
         'total xylenes':130000,
         '1,2-dichloroethane':500,
         'mtbe':4400000,
         'methyl-t-butyl ether':4400000,
         'trph':460000,
         'arsenic':2100,
         'cadmium':82000,
         'chromium':210,
         'lead':400000}

Leachabilities = {'acenapthene':2.1,
         'acenaphthylene':27,
         'anthracene':2500,
         'benzo(a)anthracene': .8,
         'benzo(a)pyrene':8,
         'benzo(b)fluoranthene':2.4,
         'benzo(g,h,i)perylene':32000,
         'benzo(k)fluoranthene':24,
         'chrysene':77,
         'dibenzo(a,h)anthracene':.7,
         'fluoranthene':1200,
         'fluorene':160,
         'indeno(1,2,3-cd)pyrene':6.6,
         '1-methylnaphthalene':3.1,
         '2-methylnaphthalene':8.5,
         'naphthalene':1.2,
         'phenanthrene':250,
         'pyrene':880,
         'benzene':.007,
         'ethylbenzene':.6,
         'toluene':.5,
         'total xylenes':.2,
         '1,2-dichloroethane':0.01,
         'mtbe':.090,
         'methyl-t-butyl ether':.090,
         'trph':340,
         'arsenic':2.1,
         'cadmium':7.5,
         'chromium':38,
         'lead':400}

import csv
import sys
from operator import itemgetter
import os
os.system("mode con: cols=220")
import fileinput
try:
        from GTCLdictcleaned import GCTLs
        print('import of expanded gctls successful')
except:
        print('import of expanded gctl dictionary did not work')

try:
        from Leachability import Leachability as Leachabilities
        print('import of expanded leachability successful')        
except:
        print('import of expanded leachability dictionary did not work')

try:
        from DERCTLs import DERCTLs as SCTLs
        print('import of expanded sctls successful')
except:
        print('import of expanded sctl dictionary did not work')
        
HitSummary = []
GCTLSummary = []
NADCSummary = []
LeachabilitySummary = []
SCTLSummary = []
DetectionSummary = []
EDD = []


def display_table(list):
        
        sortedlist = sorted(list,key=itemgetter('Parameter'))
        sortedlist = sorted(sortedlist,key=itemgetter('TestSite_Name'))
        for count, row in enumerate(sortedlist):
                GCTL = ''
                NADC = ''
                GCTL_Exceedance = ''
                NADC_Exceedance = ''
                Leachability_exceedance = ''
                SCTL_Exceedance = ''
                Detection = ''
                SCTL = ''
                Leachability = ''
                                
                if row['Parameter'] in ['SAMPLE DEPTH','Dissolved Oxygen', 'Specific Conductance', 'Temperature, Water', 'Turbidity', 'pH']:
                        continue
                if row['Qualifier']!= 'U':
                        Detection = 'Detection'
                else:
                        Detection = ''
                   
                
                if row['Matrix'] == 'W':        # Everything within this if is for water samples
                        if count == 0:
                            print("{:10}{:<20}{:<25}{:<39}{:>7}{:^5}{:<10}{:^10}{:^10}{:^25}{:^25}{:^25}".format \
                      ('Fac_ID','Loc_Name', 'Samp_Date','Param','R','Q',\
                       'U', 'GCTL ug/L', 'NADC ug/L', 'Dtd', 'GCTL Ex', 'NADC Ex'))
                        if GCTLs.get(row['Parameter'].lower()) != None:         # if a GCTL is in the GCTL list set its value to the var GCTL
                                GCTL = GCTLs.get(row['Parameter'].lower())
                        else:                                                   # if a GCTL is not in the list then the GCTL variable is made blank
                                GCTL = ''
                        if NADCs.get(row['Parameter'].lower()) != None:         # if a NADC is not in the list then set ths value to the var NADC
                                NADC = NADCs.get(row['Parameter'].lower())      
                        else:
                                NADC = ''                                       # if a NADC is not in the list then the NADC variable is made blank
                        if   GCTLs.get(row['Parameter'].lower()) != None:        
                                if   row['Units'].lower() == 'ug/l' and float(row['Result']) > float(GCTLs.get(row['Parameter'].lower())):
                                        GCTL_Exceedance = 'GCTL_Exceedance'
                                      
                                        
                                elif row['Units'].lower() == 'mg/l' and float(row['Result'])*1000 > float(GCTLs.get(row['Parameter'].lower())):
                                        GCTL_Exceedance = 'GCTL_Exceedance'
                                       
                                else:
                                        GCTL_Exceedance = ''
                        

                                
                        if   NADCs.get(row['Parameter'].lower()) != None:
                                if    row['Units'].lower() == 'ug/l' and float(row['Result']) > float(NADCs.get(row['Parameter'].lower())):
                                        NADC_Exceedance = 'NADC_Exceedance'
                                       
                                        
                                elif  row['Units'].lower() == 'mg/l' and float(row['Result'])*1000 > float(NADCs.get(row['Parameter'].lower())):
                                        NADC_Exceedance = 'NADC_Exceedance'
                                        
                                else:
                                        NADC_Exceedance = ''
                        print("{:10}{:<20}{:<25}{:<39}{:>7}{:^5}{:<10}{:^10}{:^10}{:^25}{:^25}{:^25}".format \
                          (row['Facility_ID'],row['TestSite_Name'], row['Sample_Date'], row['Parameter'], row['Result'], row['Qualifier'],\
                           row['Units'], GCTL, NADC, Detection, GCTL_Exceedance, NADC_Exceedance))        
                
                elif row['Matrix'] == 'S':      # Everything within this if is for soil samples
                        if count == 0:
                            print("{:10}{:<20}{:<25}{:<39}{:>7}{:^5}{:<10}{:^10}{:^10}{:^25}{:^25}{:^25}".format \
                      ('Fac_ID','Loc_Name', 'Samp_Date','Param','R','Q',\
                       'U', 'Lech mg/Kg', 'SCTL mg/Kg', 'Dtd', 'Leachability Ex', 'SCTL Ex'))
                        if SCTLs.get(row['Parameter'].lower()) != None:
                                SCTL = SCTLs.get(row['Parameter'].lower())
                        else:
                                SCTL = ''        
                        if Leachabilities.get(row['Parameter'].lower()) != None:
                                Leachability = Leachabilities.get(row['Parameter'].lower())
                        else:
                                Leachability = ''
                        if   Leachabilities.get(row['Parameter'].lower()) != None:        
                                if   row['Units'].lower() == 'ug/kg' and float(row['Result'])/1000 > float(Leachabilities.get(row['Parameter'].lower())):
                                        Leachability_exceedance = 'Leachability'
                                      
                                        
                                elif row['Units'].lower() == 'mg/kg' and float(row['Result']) > float(Leachabilities.get(row['Parameter'].lower())):
                                        Leachability_exceedance = 'Leachability'
                                       
                                else:
                                        Leachability_exceedance = ''

                                
                        if   SCTLs.get(row['Parameter'].lower()) != None:
                                if    row['Units'].lower() == 'ug/kg' and float(row['Result'])/1000 > float(SCTLs.get(row['Parameter'].lower())):
                                        SCTL_Exceedance = 'SCTL_Exceedance'
                                       
                                        
                                elif  row['Units'].lower() == 'mg/kg' and float(row['Result']) > float(SCTLs.get(row['Parameter'].lower())):
                                        SCTL_Exceedance = 'SCTL_Exceedance'
                                        
                                else:
                                        SCTL_Exceedance = ''
                        print("{:10}{:<20}{:<25}{:<39}{:>7}{:^5}{:<10}{:^10}{:^10}{:^25}{:^25}{:^25}".format \
                                  (row['Facility_ID'],row['TestSite_Name'], row['Sample_Date'], row['Parameter'], row['Result'], row['Qualifier'],\
                                   row['Units'], Leachability, SCTL, Detection, Leachability_exceedance, SCTL_Exceedance))             
                else:                           # Anything that wasnt S or W is skipped over entirely
                        continue
                        
with fileinput.input(sys.argv[1:]) as fileset:
	reader = csv.reader(fileset)
	for count, row in enumerate(reader):
		if count < 1:
			continue
		elif count == 1:
			fieldnames = row
			print(fieldnames)
			break
	reader = csv.DictReader(fileset,fieldnames=fieldnames)
	for count, row in enumerate(reader):
                if fileinput.filelineno() == 1 or fileinput.filelineno() == 2:
                        continue
                EDD.append(row)
                
                if row['Qualifier']!= 'U':
                        Detection = 'Detection'
                        DetectionSummary.append(row)
                else:
                        Detection = ''

                        
                if row['Matrix'] == 'W':        # Everything within this if is for water samples
                        if   GCTLs.get(row['Parameter'].lower()) != None:        
                                if   row['Units'].lower() == 'ug/l' and float(row['Result']) > float(GCTLs.get(row['Parameter'].lower())):
                                        if row not in HitSummary:
                                                HitSummary.append(row)
                                        GCTLSummary.append(row)
                                        
                                elif row['Units'].lower() == 'mg/l' and float(row['Result'])*1000 > float(GCTLs.get(row['Parameter'].lower())):
                                        if row not in HitSummary:
                                                HitSummary.append(row)
                                        GCTLSummary.append(row)
                        if   NADCs.get(row['Parameter'].lower()) != None:
                                if    row['Units'].lower() == 'ug/l' and float(row['Result']) > float(NADCs.get(row['Parameter'].lower())):
                                        upper_Exceedance = 'NADC_Exceedance'
                                        if row not in HitSummary:
                                                HitSummary.append(row)
                                        NADCSummary.append(row)
                                        if row in GCTLSummary:
                                                GCTLSummary.pop()
                                                
                                elif  row['Units'].lower() == 'mg/l' and float(row['Result'])*1000 > float(NADCs.get(row['Parameter'].lower())):
                                       
                                        if row not in HitSummary:
                                                HitSummary.append(row)
                                        NADCSummary.append(row)
                                        if row in GCTLSummary:
                                                GCTLSummary.pop()
                                         
                elif row['Matrix'] == 'S':      # Everyting within this if is for soil samples
                        if   Leachabilities.get(row['Parameter'].lower()) != None:        
                                if   row['Units'].lower() == 'ug/kg' and float(row['Result'])/1000 > float(Leachabilities.get(row['Parameter'].lower())):
                                        if row not in HitSummary:
                                                HitSummary.append(row)
                                        LeachabilitySummary.append(row)
                                        
                                elif row['Units'].lower() == 'mg/kg' and float(row['Result']) > float(Leachabilities.get(row['Parameter'].lower())):
                                        if row not in HitSummary:
                                                HitSummary.append(row)
                                        LeachabilitySummary.append(row)
                        if   SCTLs.get(row['Parameter'].lower()) != None:
                                if    row['Units'].lower() == 'ug/kg' and float(row['Result'])/1000 > float(SCTLs.get(row['Parameter'].lower())):
                                        if row not in HitSummary:
                                                HitSummary.append(row)
                                        SCTLSummary.append(row)
                                        if row in LeachabilitySummary:
                                                LeachabilitySummary.pop()
                                                
                                elif  row['Units'].lower() == 'mg/kg' and float(row['Result']) > float(SCTLs.get(row['Parameter'].lower())):
                                        if row not in HitSummary:
                                                HitSummary.append(row)
                                        SCTLSummary.append(row)
                                        if row in LeachabilitySummary:
                                                LeachabilitySummary.pop()

                else:                           # This is for anything that isnt S or W, it gets skipped over silently
                        continue

print(" HitSummary ".center(220,"*"))
print('')
display_table(HitSummary)
print('')
print('')
print(' GCTL_Exceedances '.center(220,'*'))
print('')
display_table(GCTLSummary)
print('')
print('')
print(' NADC_Exceedance '.center(220,'*'))
print('')
display_table(NADCSummary)
print('')
print('')
print(' SCTLSummary '.center(220,'*'))
print('')
display_table(SCTLSummary)
print('')
print('')
print(' LeachabilitySummary '.center(220,'*'))
print('')
display_table(LeachabilitySummary)
print('')
print('')
print(' DetectionSummary '.center(220,'*'))
print('')
display_table(DetectionSummary)
print('')
print('')

input('Press Enter to close')

	
