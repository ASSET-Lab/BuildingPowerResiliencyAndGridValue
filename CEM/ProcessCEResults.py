#Michael Craig
#October 4, 2016
#Process CE results by: 1) save new builds, 2) add new builds to gen fleet, 
#3) determine which units retire due to economics

import copy, os, random, pandas as pd, numpy as np
from CreateFleetForCELoop import *
from GAMSAuxFuncs import *

########### STORE BUILD DECISIONS FROM CAPACITY EXPANSION ######################
#Inputs: running list of CE builds (2d list), CE model output as GAMS object, 
#curr CE year
#Outputs: new gen builds by technology (list of tuples of (techtype, # builds))
def saveCEBuilds(capacExpModel,resultsDir,currYear):
    newGenerators = extract1dVarResultsFromGAMSModel(capacExpModel,'vN')
    # print('New generators:',newGenerators)
    newStoECap = extract1dVarResultsFromGAMSModel(capacExpModel,'vEneBuiltSto')
    newStoPCap = extract1dVarResultsFromGAMSModel(capacExpModel,'vPowBuiltSto')
    newLines = extract1dVarResultsFromGAMSModel(capacExpModel,'vLinecapacnew')   
    for n,d in zip(['vN','vEneBuiltSto','vPowBuiltSto','vLinecapacnew'],[newGenerators,newStoECap,newStoPCap,newLines]):
        pd.Series(d).to_csv(os.path.join(resultsDir,n+str(currYear)+'.csv'))
    return newGenerators,newStoECap,newStoPCap,newLines
                
########### ADD CAPACITY EXPANSION BUILD DECISIONS TO FLEET ####################
#Adds generators to fleet
def addNewGensToFleet(genFleet,newGenerators,newStoECap,newStoPCap,newTechs,currYear):
    if 'SiteLocUsed' not in genFleet.columns: genFleet['SiteLocUsed'] = False
    for tech,newBuilds in newGenerators.items():
        if newBuilds>0: 
            techRow = newTechs.loc[newTechs['GAMS Symbol']==tech].copy()
            #Add new info to tech row   
            techRow['Unit ID'],techRow['YearAddedCE'],techRow['Retirement Year'],techRow['SiteLocUsed'] = '1',currYear,currYear+techRow['Lifetime(years)'],False
            techRow['On Line Year'],techRow['Retired'],techRow['YearRetiredByCE'],techRow['YearRetiredByAge'] = currYear,False,0,0
            #Add rows to genFleet by building full units then the remaining partial unit
            if techRow['PlantType'].values[0] != 'Hydrogen':
                #Add rows for each full build
                while newBuilds > 1: 
                    genFleet = addNewTechRowToFleet(genFleet,techRow)    
                    newBuilds -= 1
                #Add row for partial build
                techRow['Capacity (MW)'] *= newBuilds
                techRow['Nameplate Energy Capacity (MWh)'] *= newBuilds
                genFleet = addNewTechRowToFleet(genFleet,techRow)
            else:
                #Add seasonal storage (hydrogen) by evenly dividing added E & P capacity among new units (E capacity is separate variable)
                numNewH2Facilities = int(np.ceil(newBuilds))
                for newH2Facility in range(numNewH2Facilities):
                    techRow['Nameplate Energy Capacity (MWh)'] = newStoECap[tech]/numNewH2Facilities*1000 #1000 to go from GWh to MWh
                    techRow['Capacity (MW)'] = newStoPCap[tech]/numNewH2Facilities*1000 #1000 to go from GW to MW
                    genFleet = addNewTechRowToFleet(genFleet,techRow)
    genFleet.reset_index(inplace=True,drop=True)
    return genFleet

#Concats a new generator onto genFleet using params from techRow
def addNewTechRowToFleet(genFleet,techRow):
    # if techRow['FuelType'].values[0] not in ['Wind','Solar']: genFleet,techRow = addLocToTechRow(genFleet,techRow) #adds lat/lon location to technology
    if techRow['Latitude'].values[0] == '': 
        print('**in processresults')
        print(techRow)
        genFleet,techRow = addLocToTechRow(genFleet,techRow) #adds lat/lon location to technology

    techRow['ORIS Plant Code'] = int(genFleet['ORIS Plant Code'].max())+1
    techRow['GAMS Symbol'] = techRow['ORIS Plant Code'].astype(str) + "+" + techRow['Unit ID'].astype(str)
    genFleet = pd.concat([genFleet,techRow])
    return genFleet

#Adds lat/lon coordinate to new non-RE generators (new RE generators are built for a specific lat/lon coord)
def addLocToTechRow(genFleet,techRowDf,thermalPTs=['Combined Cycle','Coal Steam','Combustion Turbine']):
    techRow = techRowDf.iloc[0] #turn into series
    regionGens = genFleet.loc[genFleet['region']==techRow['region']]
    retired = regionGens.loc[regionGens['Retired']==True]
    retiredOpenSites = retired.loc[retired['SiteLocUsed']==False]    
    matchingPTFT = retiredOpenSites.loc[(retiredOpenSites['PlantType']==techRow['PlantType']) 
                                        & (retiredOpenSites['FuelType']==techRow['FuelType'])]
                                        
    if matchingPTFT.shape[0] > 0: #if any gens match criteria
        matchingPTFTCap = matchingPTFT.loc[matchingPTFT['Capacity (MW)'] > techRow['Capacity (MW)']]
        if matchingPTFTCap.shape[0] > 0: 
            site = matchingPTFTCap.iloc[np.random.randint(low=0,high=matchingPTFTCap.shape[0])]
        else:
            site = matchingPTFT.iloc[np.random.randint(low=0,high=matchingPTFT.shape[0])]
    else:
        matchingFT = retiredOpenSites.loc[retiredOpenSites['FuelType']==techRow['FuelType']]
        if matchingFT.shape[0] > 0:
            matchingFTCap = matchingFT.loc[matchingFT['Capacity (MW)'] > techRow['Capacity (MW)']]
            if matchingFTCap.shape[0] > 0: 
                site = matchingFTCap.iloc[np.random.randint(low=0,high=matchingFTCap.shape[0])]
            else:
                site = matchingFT.iloc[np.random.randint(low=0,high=matchingFT.shape[0])]
        else:
            matchingThermal = retiredOpenSites.loc[retiredOpenSites['PlantType'].isin(thermalPTs)]
            if matchingThermal.shape[0] > 0:
                matchingThermalCap = matchingThermal.loc[matchingThermal['Capacity (MW)'] > techRow['Capacity (MW)']]
                if matchingThermalCap.shape[0] > 0: 
                    site = matchingThermalCap.iloc[np.random.randint(low=0,high=matchingThermalCap.shape[0])]
                else:
                    site = matchingThermal.iloc[np.random.randint(low=0,high=matchingThermal.shape[0])]
            else:
                allThermal = regionGens.loc[regionGens['PlantType'].isin(thermalPTs)]
                site,lat,lon = None,allThermal['Latitude'].mean(),allThermal['Longitude'].mean() 
    if site is not None: 
        lat,lon = site['Latitude'],site['Longitude']
        genFleet.loc[site.name,'SiteLocUsed'] = True
    techRow['Latitude'],techRow['Longitude'] = lat,lon
    techRow = techRow.to_frame().T
    return genFleet,techRow

########### ADD NEW LINE CAPACITIES TO LINE LIMITS #############################
def addNewLineCapToLimits(lineLimits, newLines, gwToMW = 1000):
    for line,newCapacity in newLines.items():
        lineLimits.loc[lineLimits['GAMS Symbol']==line,'TotalCapacity'] += newCapacity*gwToMW #CE solves for GW; scale to MW
    return lineLimits

########### RETIRE UNITS BY CF #################################################        
def retireUnitsByCF(genFleet,hoursForCE,capacExpModel,currYear,prm,prmHour,retirementCFCutoff=0.001, nonRetPTs = ['Geothermal','Solar PV','Onshore Wind','Hydro','Pumped Storage','Batteries','Nuclear']):
    #Filter out already retired units + units that just came online
    genFleetOnline = genFleet.loc[genFleet['Retired']==False]

    #Get total generation from CE run for existing units (not newly built this round)
    gen = getGenExistingUnits(genFleetOnline,hoursForCE,capacExpModel)
    genNewREAtPRM = getGenNewREAtPRMHour(capacExpModel,prmHour)
    
    #Retire units based on generation
    ptEligForRetireByCF = genFleetOnline['PlantType'].unique()
    ptEligForRetireByCF = [pt for pt in ptEligForRetireByCF if pt not in nonRetPTs]
    unitsRetireCF = selectRetiredUnitsByCFAndPRM(retirementCFCutoff,gen,genNewREAtPRM,genFleetOnline,prm,prmHour,ptEligForRetireByCF,currYear)
    print('Num units & units w/ CF-based retirements after CE in ' + str(currYear) + ':' + str(len(unitsRetireCF)) + '\n',unitsRetireCF)
    
    #Mark retired units
    genFleet.loc[genFleet['GAMS Symbol'].isin(unitsRetireCF),'YearRetiredByCE'] = currYear
    genFleet.loc[genFleet['GAMS Symbol'].isin(unitsRetireCF),'Retired'] = True
    return genFleet

#Get generation from prior CE run by plant
def getGenExistingUnits(genFleetOnline,hoursForCE,capacExpModel):
    gen = pd.DataFrame(columns=genFleetOnline['GAMS Symbol'],index=hoursForCE.index)
    for rec in capacExpModel.out_db['vGen']: gen.loc[rec.key(1),rec.key(0)] = rec.level
    return gen

#Get generation from new RE units
def getGenNewREAtPRMHour(capacExpModel,prmHour):
    genNewREAtPRM = 0
    for rec in capacExpModel.out_db['vGentech']: 
        if ('Wind' in rec.key(0) or 'Solar' in rec.key(0)) and (rec.key(1) == prmHour):
            genNewREAtPRM += rec.level 
    return genNewREAtPRM*1000 #scale GWh to MWh

#Determines which units retire based on CF < threshold & PRM maintained. 
def selectRetiredUnitsByCFAndPRM(retirementCFCutoff,gen,genNewREAtPRM,genFleet,prm,prmHour,ptEligForRetireByCF,currYear):
    #Get capacity that counts towards PRM, using PRM-hour-specific generation for wind & solar
    nonRECapac = genFleet.loc[~genFleet['FuelType'].isin(['Wind','Solar'])]
    nonRECapac = nonRECapac['Capacity (MW)'].sum()

    reSymbols = genFleet.loc[genFleet['FuelType'].isin(['Wind','Solar'])]['GAMS Symbol']
    reGen = gen[reSymbols]
    reGen = gen[reSymbols].loc[prmHour].sum()*1000 #scale from GWh to MWh

    #Get surplus capacity relative to PRM
    surplusCap = reGen + genNewREAtPRM + nonRECapac - prm

    #Get total CF by each unit
    totalGen = gen.sum()
    capacs = pd.Series(genFleet['Capacity (MW)'].values,index=genFleet['GAMS Symbol'])
    cfs = totalGen*1000/(capacs*gen.shape[0])

    #Sort plants below cutoff by operational cost; calculate cumulative sum; and drop plants until dropped capacity reaches surplus capacity
    gensEligRet = genFleet.loc[genFleet['PlantType'].isin(ptEligForRetireByCF)]
    gensEligRet = gensEligRet.loc[genFleet['YearAddedCE'] == 0] #only retire non-newly-built units #!= currYear] 
    gensEligRet.index = gensEligRet['GAMS Symbol']
    gensEligRet['CF'] = cfs
    gensEligRet.sort_values('CF',inplace=True,ascending=True)
    gensEligRet = gensEligRet.loc[gensEligRet['CF']<retirementCFCutoff]
    gensEligRet.sort_values('OpCost($/MWh)',inplace=True,ascending=False)
    gensEligRet['CumCap(MW)'] = gensEligRet['Capacity (MW)'].cumsum()
    gensToDrop = gensEligRet.loc[gensEligRet['CumCap(MW)']<surplusCap]

    return gensToDrop.index.values