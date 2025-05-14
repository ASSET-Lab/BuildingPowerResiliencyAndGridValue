import pandas as pd,os,numpy as np

#Add cooling types from EIA Form 860, 6_2_EnviroEquip.xlsx for WECC or NY.
#In WECC, all units at a given power plant have same cooling type except Scattergood (404), Comanche (470), Tracy (2336), and Dave Johnston (4158).  
#We use the dominant or more recent cooling type for each of these plants (404: DC, 470: RI, 2336: RI, 4158: ON).
#In NY, all units at given power plant have same cooling type except Ravenswood (61130). 
#We use the dominant cooling type for this plant (ON).
def addCoolingTypes(genFleet, region, deratedPTs=['Combined Cycle','Coal Steam']):
    if region != 'WECC' and region != 'NY': 
        print('Cooling types not set up!')
        genFleet['CoolingType'] = 'NA' #need this for compression in SetupGeneratorFleet
    else:
        #Load 860 data
        cool860 = pd.read_excel(os.path.join('Data','EIA860','6_2_EnviroEquip_Y2021.xlsx'),sheet_name='Cooling',header=1)
        #Trim last row
        if 'NOTE:' in cool860.iloc[-1]['Utility ID']: cool860 = cool860.iloc[:-1]
        cool860['Plant Code'] = cool860['Plant Code'].astype(int)
        #Drop duplicates in cool860 so only use first row of each plant (see comment above)
        cool860.drop_duplicates('Plant Code',inplace=True)
        #Add cooling type to gen fleet
        genFleet = genFleet.merge(cool860[['Plant Code','Cooling Type 1']],left_on='ORIS Plant Code',right_on='Plant Code',how='left')
        genFleet.rename(columns={'Cooling Type 1':'CoolingType'},inplace=True),genFleet.drop(['Plant Code'],axis=1,inplace=True)
        #Simplify cooling types so left with nan, RC (recirculating), ON (once through) & DC (dry cool)
        genFleet['CoolingType'] = genFleet['CoolingType'].replace({'RI':'RC','HRI':'DC'})
        #For CC and coal plants, which we derate, fill in missing cooling types with RC (the most common cooling type in WECC (120 RC vs 18/5 DC/ON)
        genFleet.loc[(genFleet['PlantType'].isin(deratedPTs)) & (genFleet['CoolingType'].isnull()),'CoolingType'] = 'RC'
    return genFleet