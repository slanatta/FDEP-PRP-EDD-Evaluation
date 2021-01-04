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
Leachability = {'acenapthene':2100,
         'acenaphthylene':27000,
         'anthracene':2500000,
         'benzo(a)anthracene': 800,
         'benzo(a)pyrene':8000,
         'benzo(b)fluoranthene':2400,
         'benzo(g,h,i)perylene':32000000,
         'benzo(k)fluoranthene':24000,
         'chrysene':77000,
         'dibenzo(a,h)anthracene':700,
         'fluoranthene':1200000,
         'fluorene':160000,
         'indeno(1,2,3-cd)pyrene':6600,
         '1-methylnaphthalene':3100,
         '2-methylnaphthalene':8500,
         'naphthalene':1200,
         'phenanthrene':250000,
         'pyrene':880000,
         'benzene':7,
         'ethylbenzene':600,
         'toluene':500,
         'total xylenes':200,
         '1,2-dichloroethane':10,
         'mtbe':90,
         'methyl-t-butyl ether':90,
         'trph':340000,
         'arsenic':2100,
         'cadmium':7500,
         'chromium':38000,
         'lead':400000}
import csv
import sys
from operator import itemgetter
import os
import fileinput
try:
        from GTCLdictcleaned import GCTLs
        input('import of expanded gctls successful,  PRESS ENTER TO CONTINUE')
except:
        input('import of expanded gctl dictionary did not work, PRESS ENTER TO CONTINUE')
try:
        from DERCTLs import DERCTLs as SCTLs
        input('import of expanded sctls successful,  PRESS ENTER TO CONTINUE')
except:
        input('import of expanded sctl dictionary did not work, PRESS ENTER TO CONTINUE')
        
HitSummary = []
GCTLSummary = []
NADCSummary = []
LeachabilitySummary = []
SCTLSummary = []
DetectionSummary = []
EDD = []


def display_table(list):
        lower_Exceedance = ''
        upper_Exceedance = ''
        lowerlevel = ''
        upperlevel = ''
        Detection = ''

        print("{:10}{:<20}{:<25}{:<39}{:>7}{:^5}{:<10}{:>10}{:>10}{:^25}{:^25}{:^25}".format \
                      ('Fac_ID','Loc_Name', 'Samp_Date','Param','R','Q',\
                       'U', 'G/L', 'NA/S', 'Dtd', 'G/L Ex', 'NA/S Ex'))
        
        sortedlist = sorted(list,key=itemgetter('Parameter'))
        sortedlist = sorted(sortedlist,key=itemgetter('TestSite_Name'))
        for count, row in enumerate(sortedlist):
                lower_Exceedance = ''
                upper_Exceedance = ''
                lowerlevel = ''
                upperlevel = ''
                if row['Parameter'] in ['SAMPLE DEPTH','Dissolved Oxygen', 'Specific Conductance', 'Temperature, Water', 'Turbidity', 'pH']:
                        continue
                if row['Qualifier']!= 'U':
                        Detection = 'Detection'
                else:
                        Detection = ''

                        


                if row['Matrix'] == 'W':
                        if GCTLs.get(row['Parameter'].lower()) != None:
                                lowerlevel = GCTLs.get(row['Parameter'].lower())
                        else:
                                lowerlevel = ''
                        if NADCs.get(row['Parameter'].lower()) != None:
                                upperlevel = NADCs.get(row['Parameter'].lower())
                        else:
                                upperlevel = ''
                        if   GCTLs.get(row['Parameter'].lower()) != None:        
                                if   row['Units'] == 'ug/L' and float(row['Result']) > float(GCTLs.get(row['Parameter'].lower())):
                                        lower_Exceedance = 'GCTL_Exceedance'
                                      
                                        
                                elif row['Units'] == 'mg/L' and float(row['Result'])*1000 > float(GCTLs.get(row['Parameter'].lower())):
                                        lower_Exceedance = 'GCTL_Exceedance'
                                       
                                else:
                                        lower_Exceedance = ''

                                
                        if   NADCs.get(row['Parameter'].lower()) != None:
                                if    row['Units'] == 'ug/L' and float(row['Result']) > float(NADCs.get(row['Parameter'].lower())):
                                        upper_Exceedance = 'NADC_Exceedance'
                                       
                                        
                                elif  row['Units'] == 'mg/L' and float(row['Result'])*1000 > float(NADCs.get(row['Parameter'].lower())):
                                        upper_Exceedance = 'NADC_Exceedance'
                                        
                                else:
                                        upper_Exceedance = ''
                                
                elif row['Matrix'] == 'S':
                        if SCTLs.get(row['Parameter'].lower()) != None:
                                upperlevel = SCTLs.get(row['Parameter'].lower())
                        else:
                                upperlevel = ''        
                        if Leachability.get(row['Parameter'].lower()) != None:
                                lowerlevel = Leachability.get(row['Parameter'].lower())
                        else:
                                lowerlevel = ''
                        if   Leachability.get(row['Parameter'].lower()) != None:        
                                if   row['Units'].lower() == 'ug/kg' and float(row['Result']) > float(Leachability.get(row['Parameter'].lower())):
                                        lower_Exceedance = 'Leachability'
                                      
                                        
                                elif row['Units'].lower() == 'mg/kg' and float(row['Result'])*1000 > float(Leachability.get(row['Parameter'].lower())):
                                        lower_Exceedance = 'Leachability'
                                       
                                else:
                                        lower_Exceedance = ''

                                
                        if   SCTLs.get(row['Parameter'].lower()) != None:
                                if    row['Units'].lower() == 'ug/kg' and float(row['Result']) > float(SCTLs.get(row['Parameter'].lower())):
                                        upper_Exceedance = 'SCTL_Exceedance'
                                       
                                        
                                elif  row['Units'].lower() == 'mg/kg' and float(row['Result'])*1000 > float(SCTLs.get(row['Parameter'].lower())):
                                        upper_Exceedance = 'SCTL_Exceedance'
                                        
                                else:
                                        upper_Exceedance = ''
                              
                else:
                        continue
                        
                
                
                print("{:10}{:<20}{:<25}{:<39}{:>7}{:^5}{:<10}{:>10}{:>10}{:^25}{:^25}{:^25}".format \
                      (row['Facility_ID'],row['TestSite_Name'], row['Sample_Date'], row['Parameter'], row['Result'], row['Qualifier'],\
                       row['Units'], lowerlevel, upperlevel, Detection, lower_Exceedance, upper_Exceedance))

                

os.system("mode con: cols=220")
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

                        
                if row['Matrix'] == 'W':
                        if GCTLs.get(row['Parameter'].lower()) != None:
                                GCTLno = GCTLs.get(row['Parameter'].lower())
                        else:
                                GCTLno = ''
                        if NADCs.get(row['Parameter'].lower()) != None:
                                NADCno = NADCs.get(row['Parameter'].lower())
                        else:
                                NADCno = ''

                        if   GCTLs.get(row['Parameter'].lower()) != None:        
                                if   row['Units'] == 'ug/L' and float(row['Result']) > float(GCTLs.get(row['Parameter'].lower())):
                                        lower_Exceedance = 'GCTL_Exceedance'
                                        if row not in HitSummary:
                                                HitSummary.append(row)
                                        GCTLSummary.append(row)
                                        
                                elif row['Units'] == 'mg/L' and float(row['Result'])*1000 > float(GCTLs.get(row['Parameter'].lower())):
                                        lower_Exceedance = 'GCTL_Exceedance'
                                        if row not in HitSummary:
                                                HitSummary.append(row)
                                        GCTLSummary.append(row)
                                else:
                                        lower_Exceedance = ''
##                        
##
                        if   NADCs.get(row['Parameter'].lower()) != None:
                                if    row['Units'] == 'ug/L' and float(row['Result']) > float(NADCs.get(row['Parameter'].lower())):
                                        upper_Exceedance = 'NADC_Exceedance'
                                        if row not in HitSummary:
                                                HitSummary.append(row)
                                        NADCSummary.append(row)
                                        if row in GCTLSummary:
                                                GCTLSummary.pop()
                                                
                                elif  row['Units'] == 'mg/L' and float(row['Result'])*1000 > float(NADCs.get(row['Parameter'].lower())):
                                        upper_Exceedance = 'NADC_Exceedance'
                                        if row not in HitSummary:
                                                HitSummary.append(row)
                                        NADCSummary.append(row)
                                        if row in GCTLSummary:
                                                GCTLSummary.pop()
                                else:
                                        upper_Exceedance = ''           
##
##                        
                elif row['Matrix'] == 'S':
                        if SCTLs.get(row['Parameter'].lower()) != None:
                                SCTLno = SCTLs.get(row['Parameter'].lower())
                        else:
                                SCTLno = ''        
                        if Leachability.get(row['Parameter'].lower()) != None:
                                Leachabilityno = Leachability.get(row['Parameter'].lower())
                        else:
                                Leachabilityno = ''
                        if   Leachability.get(row['Parameter'].lower()) != None:        
                                if   row['Units'].lower() == 'ug/kg' and float(row['Result']) > float(Leachability.get(row['Parameter'].lower())):
                                        lower_Exceedance = 'Leachability'
                                        if row not in HitSummary:
                                                HitSummary.append(row)
                                        LeachabilitySummary.append(row)
                                        
                                elif row['Units'].lower() == 'mg/kg' and float(row['Result'])*1000 > float(Leachability.get(row['Parameter'].lower())):
                                        lower_Exceedance = 'Leachability'
                                        if row not in HitSummary:
                                                HitSummary.append(row)
                                        LeachabilitySummary.append(row)
                                else:
                                        lower_Exceedance = ''
##                        
##
                        if   SCTLs.get(row['Parameter'].lower()) != None:
                                if    row['Units'].lower() == 'ug/kg' and float(row['Result']) > float(SCTLs.get(row['Parameter'].lower())):
                                        upper_Exceedance = 'SCTL_Exceedance'
                                        if row not in HitSummary:
                                                HitSummary.append(row)
                                        SCTLSummary.append(row)
                                        if row in LeachabilitySummary:
                                                LeachabilitySummary.pop()
                                                
                                elif  row['Units'].lower() == 'mg/kg' and float(row['Result'])*1000 > float(SCTLs.get(row['Parameter'].lower())):
                                        upper_Exceedance = 'SCTL_Exceedance'
                                        if row not in HitSummary:
                                                HitSummary.append(row)
                                        SCTLSummary.append(row)
                                        if row in LeachabilitySummary:
                                                LeachabilitySummary.pop()
                                else:
                                        upper_Exceedance = ''           
                                
                else:
                        continue

                        
                
#                print(fileinput.filename(),fileinput.fileno(),fileinput.filelineno(), fileinput.lineno(), fileinput.isfirstline(), row['TestSite_Name'])




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

	
