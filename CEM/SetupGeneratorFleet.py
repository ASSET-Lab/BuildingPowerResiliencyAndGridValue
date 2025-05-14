#Michael Craig, 7 July 2020

import csv, os, copy, operator, random, pandas as pd, numpy as np

def setupGeneratorFleet(region,startYear,fuelPrices,stoEff,stoMinSOC,stoFTLabels,stoInCE):
    #Import NEEDS (base fleet) and strip down to 1 fuel
    # genFleet = pd.read_excel(os.path.join('Data','needs_v6_06-30-2020.xlsx'),sheet_name='NEEDS v6_active',header=0)
    genFleet = importNEEDS()
    genFleet['FuelType'] = genFleet['Modeled Fuels'].str.split(',', expand=True)[0]
    #Import and extract data from EIA 860
    genFleet = addEIA860Data(genFleet,region,stoFTLabels) 
    #Model solar thermal facilities as solar PV
    if 'Solar Thermal' in genFleet['PlantType'].unique():
        print('Modeling all solar thermal facilities as solar PV')
        genFleet.loc[genFleet['PlantType']=='Solar Thermal','PlantType'] = 'Solar PV'
    #Add parameters
    genFleet.loc[genFleet['FuelType'].isin(stoFTLabels),'Efficiency'] = stoEff
    genFleet.loc[genFleet['FuelType'].isin(stoFTLabels),'Minimum Energy Capacity (MWh)'] = stoMinSOC
    genFleet = addFuelPrices(genFleet,startYear,fuelPrices)
    genFleet = addEmissionsRates(genFleet) 
    #Remove storage if not including in CE
    if not stoInCE: genFleet.drop(genFleet.loc[genFleet['FuelType'].isin(stoFTLabels)].index,inplace=True)  
    return genFleet

def compressAndAddSizeDependentParams(genFleet,compressFleet,regElig,contFlexInelig,regCostFrac,stoFTLabels,stoPTLabels):
    #Combine generators to reduce size of generator fleet
    if compressFleet == True: 
        genFleet,compressedGens = compressNonStorageUnits(genFleet)
        if genFleet.loc[genFleet['FuelType'].isin(stoFTLabels)].shape[0]>0: 
            genFleet,compressedGens = compressStorageUnits(genFleet,stoFTLabels,compressedGens)
        #Save all generators that were compressed
        compressedGens = pd.concat(compressedGens)
    #Add other parameters
    genFleet = addUnitCommitmentParameters(genFleet,'PhorumUCParameters.csv') 
    genFleet = addUnitCommitmentParameters(genFleet,'StorageUCParameters.csv')
    genFleet = addRandomOpCostAdder(genFleet)
    genFleet = addVOM(genFleet,stoPTLabels)  # genFleet = addVOMandFOM(genFleet,stoPTLabels) 
    genFleet = calcOpCost(genFleet)
    genFleet = addRegResCostAndElig(genFleet,regElig,regCostFrac)
    genFleet = addReserveEligibility(genFleet,contFlexInelig)
    #Add retirement tracking columns
    genFleet['Retired'] = False
    for c in ['YearAddedCE','YearRetiredByCE','YearRetiredByAge']: genFleet[c] = 0
    #Add unique code used as GAMS symbol
    genFleet['GAMS Symbol'] = genFleet['ORIS Plant Code'].astype(str) + "+" + genFleet['Unit ID'].astype(str)
    #Drop duplicates if present (1 duplicated unit in EI NEEDS data)
    if genFleet.duplicated(subset='GAMS Symbol').sum()>0:
        print('Dropped ' + str(genFleet.duplicated(subset='GAMS Symbol').sum()) + ' duplicated units')
        genFleet.drop_duplicates(subset='GAMS Symbol',inplace=True,ignore_index=True)
    return genFleet,compressedGens

################################################################################
def importNEEDS():
    active = pd.read_excel(os.path.join('Data','NEEDS rev 02-14-2023.xlsx'),sheet_name='NEEDS v621_Active',header=0)
    retiring = pd.read_excel(os.path.join('Data','NEEDS rev 02-14-2023.xlsx'),sheet_name='NEEDS v621_Retired_Through2028',header=0)
    new = pd.read_excel(os.path.join('Data','NEEDS rev 02-14-2023.xlsx'),sheet_name='NEEDSv621_NewCapacity_Hardwired',header=0)
    genFleet = pd.concat([active,retiring,new])
    return genFleet

def addEIA860Data(genFleet,region,stoFTLabels,missingStoDuration=4,missingPSStoDuration=8): #statesForAnalysis
    gens,plants,storage = importEIA860()
    genFleet = genFleet.merge(plants[['Plant Code','Latitude','Longitude']],left_on='ORIS Plant Code',right_on='Plant Code',how='left')
    genFleet = fillMissingCoords(genFleet)
    genFleet = genFleet.merge(storage[['Plant Code','Generator ID','Nameplate Energy Capacity (MWh)','Maximum Charge Rate (MW)','Maximum Discharge Rate (MW)','Technology']],left_on=['ORIS Plant Code','Unit ID'],right_on=['Plant Code','Generator ID'],how='left')
    #Isolate area of interest
    needsRegions = mapInterconnToNEEDSRegions()[region]
    genFleet = genFleet.loc[genFleet['Region Name'].str.contains('|'.join(needsRegions))]
    genFleet.reset_index(inplace=True,drop=True)
    #Replace storage plant type from NEEDS w/ technology from 860
    stoRowsMatched = genFleet.loc[(genFleet['FuelType'].isin(stoFTLabels)) & (genFleet['Nameplate Energy Capacity (MWh)'].isnull() == False)]
    genFleet.loc[stoRowsMatched.index,'PlantType'] = stoRowsMatched['Technology']
    #Fill in missing storage parameters
    #Fill in rest of storage units (and PS charge & discharge rate)
    stoRowsMissingMatch = genFleet.loc[(genFleet['FuelType'].isin(stoFTLabels)) & (genFleet['Nameplate Energy Capacity (MWh)'].isnull())]
    genFleet.loc[stoRowsMissingMatch.index,'Maximum Charge Rate (MW)'] = genFleet['Capacity (MW)']
    genFleet.loc[stoRowsMissingMatch.index,'Maximum Discharge Rate (MW)'] = genFleet['Capacity (MW)']
    genFleet.loc[stoRowsMissingMatch.index,'Nameplate Energy Capacity (MWh)'] = genFleet['Capacity (MW)'] * missingStoDuration  
    #Pumped hydropower (as of Y2018 release) is not included in storage dataset. Use 8 hour duration for them.
    genFleet.loc[genFleet['FuelType'] == 'Pumped Storage','Nameplate Energy Capacity (MWh)'] = genFleet['Capacity (MW)'] * missingPSStoDuration    
    genFleet['Nameplate Energy Capacity (MWh)'] = genFleet['Nameplate Energy Capacity (MWh)'].astype(float)
    return genFleet

#Fill generators with missing lat/lon (not in EIA 860) w/ coords from other plants in same county or state
def fillMissingCoords(genFleet):
    missingCoordRows = genFleet.loc[genFleet['Latitude'].isna()]
    for idx,row in missingCoordRows.iterrows():
        county,state = row['County'],row['State Name']
        otherRows = genFleet.loc[(genFleet['County']==county) & (genFleet['State Name']==state)]
        otherRowsWithCoords = otherRows.loc[otherRows['Latitude'].isna()==False]
        if otherRowsWithCoords.shape[0]==0: otherRowsWithCoords = genFleet.loc[genFleet['State Name']==state] 
        lat,lon = otherRowsWithCoords['Latitude'].median(),otherRowsWithCoords['Longitude'].median()
        genFleet.loc[idx,'Latitude'],genFleet.loc[idx,'Longitude'] = lat,lon
    return genFleet

def mapInterconnToNEEDSRegions(): #make sure list is returned in dict
    return {'ERCOT':['ERC'],'WECC':['WEC'],'EI':['FRCC','MIS','NENG','NY','PJM','SPP','S_C','S_D','S_SOU','S_VACA'],
            'NY':['NY_Z']}

# def importEIA860():
#     dir860 = os.path.join('Data','EIA860')
#     gens860 = pd.read_excel(os.path.join(dir860,'3_1_Generator_Y2018.xlsx'),sheet_name='Operable',header=1)
#     sto860 = pd.read_excel(os.path.join(dir860,'3_4_Energy_Storage_Y2018.xlsx'),sheet_name='Operable',header=1)
#     plants860 = pd.read_excel(os.path.join(dir860,'2___Plant_Y2018.xlsx'),sheet_name='Plant',header=1)
#     #Remove solar thermal from storage (no storage values given)
#     sto860 = sto860.loc[sto860['Technology'] != 'Solar Thermal with Energy Storage']
#     return gens860,plants860,sto860

def importEIA860():
    dir860 = os.path.join('Data','EIA860')
    gens860 = pd.read_excel(os.path.join(dir860,'3_1_Generator_Y2021.xlsx'),sheet_name='Operable',header=1)
    sto860 = pd.read_excel(os.path.join(dir860,'3_4_Energy_Storage_Y2021.xlsx'),sheet_name='Operable',header=1)
    plants860 = pd.read_excel(os.path.join(dir860,'2___Plant_Y2021.xlsx'),sheet_name='Plant',header=1)
    #Remove solar thermal from storage (no storage values given)
    sto860 = sto860.loc[sto860['Technology'] != 'Solar Thermal with Energy Storage']
    #Rename batteries to Battery Storage
    sto860['Technology'] = sto860['Technology'].replace({'Batteries':'Battery Storage'})
    return gens860,plants860,sto860
################################################################################

################################################################################
#COMPRESS FLEET BY COMBINING SMALL UNITS BY REGION
def compressNonStorageUnits(genFleetAll,maxSizeToCombine=200,maxCombinedSize=10000,firstYr=1975,lastYr=2026,stepYr=10):
    genFleetAll['FuelType2'] = genFleetAll['FuelType']
    genFleetAll.loc[genFleetAll["PlantType"] == "Combined Cycle", "FuelType2"] = "Combined Cycle"
    genFleetAll['hrGroup'] = genFleetAll.groupby(['region','FuelType2'])['Heat Rate (Btu/kWh)'].transform(lambda x: pd.qcut(x, 4, duplicates='drop'))
    # genFleetAll['hrGroup'] = genFleetAll.groupby(['region', 'FuelType'])['Heat Rate (Btu/kWh)'].transform(lambda x: pd.qcut(x.rank(method='first'), 4, duplicates='drop'))
    startRegionCap,startFuelCap = genFleetAll.groupby(['region']).sum()['Capacity (MW)'],\
                                              genFleetAll.groupby(['FuelType2']).sum()['Capacity (MW)']
    rowsToDrop,rowsToAdd,compressedGens = list(),list(),list()
    for region in genFleetAll['region'].unique():
        genFleet = genFleetAll.loc[genFleetAll['region']==region]
        for fuel in ['Distillate Fuel Oil','Natural Gas','Combined Cycle','Residual Fuel Oil','Bituminous','Subbituminous','Lignite']: #'Geothermal'
            genFleetFuel = genFleet.loc[genFleet['FuelType2'] == fuel]
            for hr in genFleetFuel['hrGroup'].unique():
                startHrCap = genFleetFuel.groupby(['Heat Rate (Btu/kWh)']).sum()['Capacity (MW)']
                fuelRows = genFleetFuel.loc[(genFleetFuel['hrGroup'] == hr) & (genFleetFuel['FuelType2'] == fuel) & (genFleetFuel['Capacity (MW)'] < maxSizeToCombine)]
                yearIntervals = [yr for yr in range(firstYr,lastYr,stepYr)]
                for endingYear in yearIntervals:
                    beginningYear = 0 if endingYear == firstYr else endingYear-stepYr
                    fuelRowsYears = fuelRows.loc[(fuelRows['On Line Year']>beginningYear) & (fuelRows['On Line Year']<=endingYear)]
                    if fuelRowsYears.shape[0]>1:
                        runningCombinedSize,rowsToCombine = 0,list()
                        for index, row in fuelRowsYears.iterrows():
                            if (runningCombinedSize + row['Capacity (MW)'] > maxCombinedSize):
                                newRow,idxsToDrop,compressedGens = aggregateRows(rowsToCombine,compressedGens)
                                rowsToAdd.append(newRow),rowsToDrop.extend(idxsToDrop)
                                runningCombinedSize,rowsToCombine = row['Capacity (MW)'],[row]
                            else:
                                runningCombinedSize += row['Capacity (MW)']
                                rowsToCombine.append(row)
                        if len(rowsToCombine)>1:
                            newRow,idxsToDrop,compressedGens = aggregateRows(rowsToCombine,compressedGens)
                            rowsToAdd.append(newRow),rowsToDrop.extend(idxsToDrop)
                endHrCap = genFleetFuel.groupby(['Heat Rate (Btu/kWh)']).sum()['Capacity (MW)']
                assert (startHrCap.astype(int).equals(endHrCap.astype(int)))
    assert(len(set(rowsToDrop))==len(rowsToDrop))
    genFleetAll.drop(index=rowsToDrop,inplace=True)
    genFleetAll = genFleetAll.append(pd.DataFrame(rowsToAdd))
    genFleetAll.reset_index(drop=True,inplace=True)
    endRegionCap,endFuelCap = genFleetAll.groupby(['region']).sum()['Capacity (MW)'], \
                                       genFleetAll.groupby(['FuelType2']).sum()['Capacity (MW)']
    assert(startRegionCap.astype(int).equals(endRegionCap.astype(int)))
    assert(startFuelCap.astype(int).equals(endFuelCap.astype(int)))

    # For the units with cost of 0:
    rowsToDrop2, rowsToAdd2 = list(), list()
    for region in genFleetAll['region'].unique():
        genFleet = genFleetAll.loc[genFleetAll['region'] == region]
        for fuel in ['Landfill Gas','MSW','Biomass','Non-Fossil Waste','Fossil Waste','Geothermal']:
            fuelRows = genFleet.loc[(genFleet['FuelType'] == fuel) & (genFleet['Capacity (MW)'] < maxSizeToCombine) & (genFleet['PlantType'] != 'Combined Cycle')]
            yearIntervals = [yr for yr in range(firstYr, lastYr, stepYr)]
            for endingYear in yearIntervals:
                beginningYear = 0 if endingYear == firstYr else endingYear - stepYr
                fuelRowsYears = fuelRows.loc[(fuelRows['On Line Year'] > beginningYear) & (fuelRows['On Line Year'] <= endingYear)]
                if fuelRowsYears.shape[0] > 1:
                    runningCombinedSize, rowsToCombine = 0, list()
                    for index, row in fuelRowsYears.iterrows():
                        if (runningCombinedSize + row['Capacity (MW)'] > maxCombinedSize):
                            newRow, idxsToDrop,compressedGens = aggregateRows(rowsToCombine,compressedGens)
                            rowsToAdd2.append(newRow), rowsToDrop2.extend(idxsToDrop)
                            runningCombinedSize, rowsToCombine = row['Capacity (MW)'], [row]
                        else:
                            runningCombinedSize += row['Capacity (MW)']
                            rowsToCombine.append(row)
                    if len(rowsToCombine) > 1:
                        newRow, idxsToDrop, compressedGens = aggregateRows(rowsToCombine,compressedGens)
                        rowsToAdd2.append(newRow), rowsToDrop2.extend(idxsToDrop)
    assert (len(set(rowsToDrop2)) == len(rowsToDrop2))
    genFleetAll.drop(index=rowsToDrop2, inplace=True)
    genFleetAll = genFleetAll.append(pd.DataFrame(rowsToAdd2))
    genFleetAll.reset_index(drop=True, inplace=True)
    endRegionCap, endFuelCap = genFleetAll.groupby(['region']).sum()['Capacity (MW)'], genFleetAll.groupby(['FuelType2']).sum()['Capacity (MW)']
    assert (startRegionCap.astype(int).equals(endRegionCap.astype(int)))
    assert (startFuelCap.astype(int).equals(endFuelCap.astype(int)))

    return genFleetAll,compressedGens
                
def aggregateRows(rowsToCombine,compressedGens):
    #Create 1 aggregate row representing rowsToCombine
    rowsToCombine = pd.DataFrame(rowsToCombine)
    capacWts = rowsToCombine['Capacity (MW)']/rowsToCombine['Capacity (MW)'].sum()
    newRow = rowsToCombine.iloc[0].copy()
    newRow['Capacity (MW)'] = rowsToCombine['Capacity (MW)'].sum()
    for p in ['CO2EmRate(lb/MMBtu)','Heat Rate (Btu/kWh)','Latitude','Longitude','On Line Year']: #'NOxEmRate(lb/MMBtu)','SO2EmRate(lb/MMBtu)'
        newRow[p] = (rowsToCombine[p]*capacWts).sum() if p != 'On Line Year' else int((rowsToCombine[p]*capacWts).sum())
    newRow['Unit ID'] = str(newRow['Unit ID'])+'COMBINED'

    #Save units that are being compressed
    rowsToCombine['UnitCompressedInto'] = newRow['ORIS Plant Code'].astype(str) + "+" + newRow['Unit ID']
    compressedGens.append(rowsToCombine)
    return newRow,rowsToCombine.index,compressedGens

#Combine all storage types into a single storage unit, ignoring diff storage techs
def compressStorageUnits(genFleetAll,stoFTLabels,compressedGens):
    #Compile total values to check fleet at end
    startRegionCap = genFleetAll.groupby(['region']).sum()['Capacity (MW)']
    stoFleet = genFleetAll.loc[genFleetAll['FuelType'].isin(stoFTLabels)]
    startEStoCap = stoFleet.groupby(['region']).sum()['Nameplate Energy Capacity (MWh)']
    startStoCap = stoFleet.groupby(['region']).sum()['Capacity (MW)']
    #Compress storage rows
    rowsToDrop,rowsToAdd = list(), list()
    for region in genFleetAll['region'].unique():
        stoRegionRows = stoFleet.loc[stoFleet['region'] == region]
        if stoRegionRows.shape[0] > 0:
            newRow, idxsToDrop, compressedGens = aggregateStoRows(stoRegionRows,compressedGens)
            rowsToDrop.extend(idxsToDrop),rowsToAdd.append(newRow)
    genFleetAll.drop(index=rowsToDrop, inplace=True)
    genFleetAll = genFleetAll.append(pd.DataFrame(rowsToAdd))
    genFleetAll.reset_index(drop=True, inplace=True)
    #Check total values haven't changed
    endRegionCap = genFleetAll.groupby(['region']).sum()['Capacity (MW)']
    stoFleet = genFleetAll.loc[genFleetAll['FuelType'].isin(stoFTLabels)]
    endEStoCap = stoFleet.groupby(['region']).sum()['Nameplate Energy Capacity (MWh)']
    endStoCap = stoFleet.groupby(['region']).sum()['Capacity (MW)']
    assert(startRegionCap.astype(int).equals(endRegionCap.astype(int)))
    assert(startEStoCap.astype(int).equals(endEStoCap.astype(int)))
    assert(startStoCap.astype(int).equals(endStoCap.astype(int)))

    return genFleetAll,compressedGens

def aggregateStoRows(rowsToCombine,compressedGens):
    rowsToCombine = pd.DataFrame(rowsToCombine)
    capacWts = rowsToCombine['Capacity (MW)']/rowsToCombine['Capacity (MW)'].sum()
    newRow = rowsToCombine.iloc[0].copy() #inherit properties from first row
    newRow['Capacity (MW)'] = rowsToCombine['Capacity (MW)'].sum()
    newRow['Maximum Charge Rate (MW)'] = rowsToCombine['Maximum Charge Rate (MW)'].sum()
    newRow['Maximum Discharge Rate (MW)'] = rowsToCombine['Maximum Discharge Rate (MW)'].sum()
    newRow['Nameplate Energy Capacity (MWh)'] = rowsToCombine['Nameplate Energy Capacity (MWh)'].sum()
    newRow['On Line Year'] = rowsToCombine['On Line Year'].median()
    for p in ['CO2EmRate(lb/MMBtu)','Heat Rate (Btu/kWh)']: newRow[p] = 0
    newRow['Unit ID'] = str(newRow['Unit ID'])+'COMBINED'

    #Save units that are being compressed
    rowsToCombine['UnitCompressedInto'] = newRow['ORIS Plant Code'].astype(str) + "+" + newRow['Unit ID']
    compressedGens.append(rowsToCombine)
    return newRow,rowsToCombine.index,compressedGens
################################################################################

################################################################################
#ADD VARIABLE O&M COSTS (based on plant type)
#Set which types of wind, solar, & battery techs are used from ATB 
def getATBTechDetailsForWindSolarBattery():
    #Technology types
    solarTech,windTech,offwindTech,batteryTech = 'Utility PV','Land-Based Wind','Offshore Wind','Utility-Scale Battery Storage' #'Utility PV' for regular runs, 'Utility-Scale PV-Plus-Battery' for CESM2-LE runs!
    #Tech details 
    windTechDetail,offwindTechDetail,solarTechDetail,batteryTechDetail = 'Class5', 'Class5', 'Class5', '4Hr Battery Storage'
    return solarTech,windTech,batteryTech,windTechDetail,solarTechDetail,batteryTechDetail,offwindTech,offwindTechDetail

#Get data from NREL's ATB for current year of analysis and input into new techs df
def addVOM(genFleet,stoPTLabels,scenario='Moderate',dataYear=2020):
    print('Using VOM of 0.1$/MWh for storage to eliminate simultaneous charging & discharging')

    #Import tech detail labels for ATB
    solarTech,windTech,batteryTech,windTechDetail,solarTechDetail,batteryTechDetail,offwindTech,offwindTechDetail = getATBTechDetailsForWindSolarBattery()

    #Align plant type labels from existing fleet to ATB rows using 'technology' and 'techdetail' columns in ATB CSV
    ptToATBTechAlias = {'Solar PV':solarTech,'Onshore Wind':windTech,'Battery Storage':batteryTech,'Nuclear':'Nuclear',
                    'Coal Steam':'Coal','Coal Steam CCS':'Coal','Combined Cycle':'Natural Gas','Combined Cycle CCS':'Natural Gas','Combustion Turbine':'Natural Gas',
                    'IGCC':'Coal','Geothermal':'Geothermal','Hydro':'Hydropower','Pumped Storage':'Pumped Storage Hydropower','Biomass':'Biopower',
                    'Municipal Solid Waste':'Biopower','Non-Fossil Waste':'Coal','O/G Steam':'Coal','Landfill Gas':'Biopower','Pet. Coke':'Coal',
                    'Fossil Waste':'Coal','Tires':'Coal','Fuel Cell':batteryTech,'Energy Storage':batteryTech}
    feToATBTechDetail = {'Solar PV':solarTechDetail,'Onshore Wind':windTechDetail,'Battery Storage':batteryTechDetail,'Nuclear':'Nuclear',
                    'Coal Steam':'newAvgCF','Coal Steam CCS':'CCS90AvgCF','Combined Cycle':'CCAvgCF','Combined Cycle CCS':'CCCCSAvgCF','Combustion Turbine':'CTAvgCF',
                    'IGCC':'IGCCAvgCF','Geothermal':'HydroBinary','Hydro':'NPD2','Pumped Storage':'NatlClass5','Biomass':'Dedicated',
                    'Municipal Solid Waste':'Dedicated','Non-Fossil Waste':'newAvgCF','O/G Steam':'newAvgCF','Landfill Gas':'Dedicated','Pet. Coke':'newAvgCF',
                    'Fossil Waste':'newAvgCF','Tires':'newAvgCF','Fuel Cell':batteryTechDetail,'Energy Storage':batteryTechDetail}

    #Load ATB CSV
    atb = pd.read_csv(os.path.join('Data','NewPlantData','ATBe.csv'),index_col=0,header=0)

    #Filter ATB rows
    atb = atb.loc[atb['core_metric_variable']==dataYear] 
    atb = atb.loc[atb['core_metric_case']=='Market'] 
    atb = atb.loc[atb['scenario']==scenario]

    #Add values for each plant type
    for pt in genFleet['PlantType'].unique():
        #Get rows for tech type
        techRows = atb.loc[atb['technology_alias']==ptToATBTechAlias[pt]]
        techRows = techRows.loc[techRows['techdetail']==feToATBTechDetail[pt]]

        #Extract parameters (use .iloc[0] because of redundancies in data that do not change values of our parameters of interest)
        if 'Variable O&M' in techRows['core_metric_parameter'].unique():
            vom = techRows.loc[techRows['core_metric_parameter']=='Variable O&M']['value'].iloc[0]
        elif pt in stoPTLabels:
            vom = 0.1 #0 #$/MWh 
        else:
            vom = 0 #$/MWh    
        # vom = techRows.loc[techRows['core_metric_parameter']=='Variable O&M']['value'].iloc[0] if 'Variable O&M' in techRows['core_metric_parameter'].unique() else 0 #$/MWh

        #Add parameters
        genFleet.loc[genFleet['PlantType']==pt,'VOM($/MWh)'] = vom #already in $/MWh

    #Align costs w/ target year costs
    genFleet['VOM($/MWh)'] = convertCostToTgtYr('vom',genFleet['VOM($/MWh)'])
    return genFleet
################################################################################

################################################################################
#ADD UNIT COMMITMENT PARAMETERS
#Based on fuel and plant type; data from PHORUM
def addUnitCommitmentParameters(genFleet,fname):
    ucData = readCSVto2dList(os.path.join('Data',fname))
    for ucHeader in ['MinDownTime(hrs)','RampRate(MW/hr)','StartCost($)','MinLoad(MWh)']:
        if ucHeader not in genFleet.columns: #only initialize once
            genFleet[ucHeader] = genFleet['Capacity (MW)'] if (ucHeader in ['RampRate(MW/hr)','MinLoad(MWh)']) else 0
        phorumParamName = mapHeadersToPhorumParamNames()[ucHeader]
        for index,row in genFleet.iterrows():
            (fuel,plantType,size) = (row['FuelType'],row['PlantType'],float(row['Capacity (MW)']))
            phorumValue = getMatchingPhorumValue(ucData,fuel,plantType,size,phorumParamName)
            if phorumValue is not None: #input files don't have all plant types
                valToAdd = phorumValue if ucHeader == 'MinDownTime(hrs)' else phorumValue*size
                if ucHeader == 'StartCost($)': valToAdd = convertCostToTgtYr('startup',valToAdd)
                genFleet.loc[index,ucHeader]=valToAdd
    return genFleet

#Read CSV to 2d list. Input: full file name including dir (str). Output: 2d list.
def readCSVto2dList(fileNameWithDir):
    with open(fileNameWithDir,'r') as f:
        f = csv.reader(f)
        f = list(f)
    return f

def getMatchingPhorumValue(ucData,fuel,plantType,size,paramName):
    if plantType == 'Fuel Cell': plantType = 'Combustion Turbine'
    fuel = mapFuels()[fuel]
    phorumPropertyNameCol = ucData[0].index('PropertyName')
    phorumFuelCol = ucData[0].index('Fuel')
    phorumPlantTypeCol = ucData[0].index('PlantType')
    phorumLowerSizeCol = ucData[0].index('LowerPlantSizeLimit')
    phorumUpperSizeCol = ucData[0].index('UpperPlantSizeLimit')
    phorumValueCol = ucData[0].index('PropertyValue')
    phorumProperties = [row[phorumPropertyNameCol] for row in ucData[1:]]
    phorumFuels = [row[phorumFuelCol] for row in ucData[1:]]
    phorumPlantTypes = [row[phorumPlantTypeCol] for row in ucData[1:]]
    phorumLowerSizes = [int(row[phorumLowerSizeCol]) for row in ucData[1:]]
    phorumUpperSizes = [int(row[phorumUpperSizeCol]) for row in ucData[1:]]
    phorumValues = [float(row[phorumValueCol]) for row in ucData[1:]]
    for idx in range(len(phorumProperties)):
        if (phorumProperties[idx] == paramName and phorumFuels[idx] == fuel and 
            (phorumPlantTypes[idx] in plantType or phorumPlantTypes[idx] == 'All') and 
            (phorumLowerSizes[idx] <= size and phorumUpperSizes[idx] > size)):
            return float(phorumValues[idx])

#Return dictionary of fleet fuel : UC fuel
def mapFuels():
    return {'Bituminous': 'Coal', 'Petroleum Coke': 'Pet. Coke',
        'Subbituminous': 'Coal', 'Lignite': 'Coal', 'Natural Gas': 'NaturalGas',
        'Distillate Fuel Oil': 'Oil', 'Hydro': 'Hydro', 'Landfill Gas': 'LF Gas',
        'Biomass': 'Biomass', 'Solar': 'Solar', 'Non-Fossil Waste': 'Non-Fossil',
        'MSW': 'MSW', 'Pumped Storage': 'Hydro', 'Residual Fuel Oil': 'Oil',
        'Wind': 'Wind', 'Nuclear Fuel': 'Nuclear', 'Coal': 'Coal','Energy Storage':'Storage',
        'Hydrogen':'Storage','Storage':'Storage','Fossil Waste':'Oil','Tires':'Non-Fossil',
        'Waste Coal':'Coal','Geothermal':'Geothermal'}
    #EIA860 fuels
    # fleetFuelToPhorumFuelMap = {'BIT':'Coal','PC':'Pet. Coke','SGC':'Coal',
    #         'SUB':'Coal','LIG':'Coal','RC':'Coal','NG':'NaturalGas','OG':'NaturalGas',
    #         'BFG':'NaturalGas','DFO':'Oil','RFO':'Oil','KER':'Oil','BLQ':'Oil',
    #         'WAT':'Hydro','LFG':'LF Gas','OBG':'LF Gas','WH':'NaturalGas','PUR':'NaturalGas',
    #         'AB':'Biomass','WDS':'Biomass','SUN':'Solar','Non-Fossil Waste':'Non-Fossil',
    #         'MSW':'MSW','WND':'Wind','NUC':'Nuclear','MWH':'MWH'}

def mapHeadersToPhorumParamNames():
    return {'MinDownTime(hrs)':'Min Down Time','RampRate(MW/hr)':'Ramp Rate',
            'StartCost($)':'Start Cost','MinLoad(MWh)':'Min Stable Level'}
################################################################################

################################################################################
#ADD FUEL PRICES
def addFuelPrices(genFleet,currYear,fuelPrices):
    if currYear > 2050: currYear = 2050
    fuelPrices = fuelPrices.loc[currYear] if currYear in fuelPrices.index else fuelPrices.iloc[-1]
    fuelPrices = convertCostToTgtYr('fuel',fuelPrices)
    prices = fuelPrices.to_dict()
    fuelMap = mapFuelsToAEOPrices()
    genFleet['FuelPrice($/MMBtu)'] = [prices[fuelMap[f]] if (f in fuelMap and fuelMap[f] in prices) else (prices[f] if f in prices else 0) for f in genFleet['FuelType']]
    return genFleet

def mapFuelsToAEOPrices():
    return {'Bituminous': 'Steam Coal', 'Petroleum Coke': 'Steam Coal','Coal':'Steam Coal',
        'Subbituminous': 'Steam Coal', 'Lignite': 'Steam Coal','Nuclear Fuel': 'Uranium'}

#Convert dollar years to 2012 dollars
#CPI from Minneapolis Fed, https://www.minneapolisfed.org/about-us/monetary-policy/inflation-calculator/consumer-price-index-1913-
#Inputs: name of parameter of dollar year to convert, parameter value (cost)
#Outputs: cost in target year dollars
def convertCostToTgtYr(paramName,cost,targetDollarYear=2022):    
    paramDollarYears = {'startup':2011,'vom':2020,'fom':2020,'occ':2020,'fuel':2022}
    cpiValues = {2011:224.9,2020:258.8,2021:271,2022:294.4}
    return doConversion(paramName,cost,paramDollarYears,targetDollarYear,cpiValues)

# #Convert dollar year
def doConversion(paramName,cost,paramDollarYears,targetDollarYear,cpiValues):
    paramDollarYear = paramDollarYears[paramName]
    (cpiTgtYear,cpiParamYear) = (cpiValues[targetDollarYear],cpiValues[paramDollarYear])
    return cost*cpiTgtYear/cpiParamYear
################################################################################

################################################################################
#ADD RANDOM OP COST ADDER TO FLEET IN NEW COLUMN
#This function deprecated as of 9/19/22 by Michael Craig - increases computation time. 
#If someone chooses to add at later date, note that new generators added to fleet by CE model 
#do not have a random op cost added. (See ProcessCEResults.py, addNewTechRowToFleet function).
#If use value of 0.05, max addition to op cost of gen in fleet is 0.19%.  
def addRandomOpCostAdder(genFleet,ocAdderMin=0,ocAdderMax=0):
    random.seed()
    genFleet['RandOpCostAdder($/MWh)'] = pd.Series(np.random.uniform(ocAdderMin,ocAdderMax,genFleet.shape[0]))
    return genFleet
################################################################################

################################################################################
def calcOpCost(genFleet):
    genFleet['OpCost($/MWh)'] = genFleet['FuelPrice($/MMBtu)']*genFleet['Heat Rate (Btu/kWh)']/1000+genFleet['VOM($/MWh)']+genFleet['RandOpCostAdder($/MWh)']
    return genFleet
################################################################################

################################################################################
#ADD REG OFFER COST AND ELIGIBILITY
def addRegResCostAndElig(genFleet,regElig,regCostFrac):
    genFleet['RegOfferElig'] = 0
    genFleet.loc[genFleet['PlantType'].str.contains('|'.join(regElig)),'RegOfferElig'] = 1
    genFleet['RegOfferCost($/MW)'] = regCostFrac*genFleet['OpCost($/MWh)']*genFleet['RegOfferElig']
    return genFleet
################################################################################

################################################################################
def addReserveEligibility(genFleet,contFlexInelig):
    genFleet['FlexOfferElig'],genFleet['ContOfferElig'] = 1,1
    genFleet.loc[genFleet['FuelType'].str.contains('|'.join(contFlexInelig)),'FlexOfferElig'] = 0
    genFleet.loc[genFleet['FuelType'].str.contains('|'.join(contFlexInelig)),'ContOfferElig'] = 0
    return genFleet
################################################################################
    
################################################################################
def addEmissionsRates(genFleet):
    emissionRates = pd.read_excel(os.path.join('Data','co2_vol_mass_updated.xls'),sheet_name='Sheet1',index_col=0,skiprows=2,usecols='A,F')
    emissionRates = emissionRates[emissionRates.columns[0]] #convert to Series
    fuelMap = fuelMapEmissions()
    genFleet['CO2EmRate(lb/MMBtu)'] = [emissionRates[fuelMap[f]] if (f in fuelMap and fuelMap[f] in emissionRates) else (emissionRates[f] if f in emissionRates else 0) for f in genFleet['FuelType']]
    return genFleet

def fuelMapEmissions():
    return {'MSW':'Municiple Solid Waste','Biomass':'Municiple Solid Waste','Landfill Gas':'Natural Gas',
            'Distillate Fuel Oil':'Other petroleum & miscellaneous','Residual Fuel Oil':'Other petroleum & miscellaneous',
            'Waste Coal':'Bituminous','Fossil Waste':'Other petroleum & miscellaneous','Non-Fossil Waste':'Other petroleum & miscellaneous',
            'Petroleum Coke':'Petroleum coke'}
################################################################################
