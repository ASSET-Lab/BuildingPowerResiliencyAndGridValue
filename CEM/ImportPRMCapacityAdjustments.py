import os, copy, datetime, pandas as pd, geopandas as gpd, datetime as dt, numpy as np, xarray as xr
from os import path
from ImportNonREMetVars import importTGWMetVars
from CalculateDerates import getClosestCellCoordsWithMet,shiftToTz

def importPRMCapacityAdjustments(genFleetForCE, newTechsCE, demand, prmBasis, region, nonCCReanalysis, 
                                    weatherYears, compressedGens, climateScenario, tz, defaultFOR=0.05):
    # if climateScenario != None and len(climateScenario)>1: sys.exit('PRM capacity adjustments not built for multiple CESM members!')

    #Load FOR regressions
    forsRegression,forPTMatching = loadFORRegressions()

    #Wind & solar not eligible for PRM if based on net demand because already accounted for in calculating planning reserve!
    prmEligWindSolar = 1 if prmBasis == 'demand' else 0

    #Wind and solar FORs equal default FOR
    windFOR,solarFOR = defaultFOR,defaultFOR

    #Initialize FORs
    # fors = pd.DataFrame(defaultFOR,index=demand.index,columns=genFleetForCE['GAMS Symbol'])
    # forsTechs = pd.DataFrame(defaultFOR,index=demand.index,columns=newTechsCE['GAMS Symbol'])

    #Get temperatures
    temps = importTGWMetVars(weatherYears,climateScenario)
    tVar = 'TREFHT'

    #Calculate FORs for existing units
    fors = calculateFORs(temps,genFleetForCE,compressedGens,tVar,forsRegression,forPTMatching)

    #Calculate FORs for new units
    forsTechs = calculateNewTechFORs(temps,newTechsCE,tVar,forsRegression,forPTMatching)

    #Shift from UTC to local timezone
    fors = shiftToTz(fors,weatherYears,tz)
    forsTechs = shiftToTz(forsTechs,weatherYears,tz)

    #Fill in missing indices in derates (last few hours of year due to timeshift)
    missing = demand.index.difference(fors.index)
    fors = pd.concat([fors,pd.DataFrame(0,index=missing,columns=fors.columns)])
    forsTechs = pd.concat([forsTechs,pd.DataFrame(0,index=missing,columns=forsTechs.columns)])

    return fors,windFOR,solarFOR,forsTechs,prmEligWindSolar

def loadFORRegressions():
    #Import regression-based relationships for TDFORs. Taken from Murphy et al. (https://www.sciencedirect.com/science/article/pii/S0306261919311870)
    forsRegression = pd.read_excel(os.path.join('Data','TDFORRelationships.xlsx'),index_col=0,header=0)
    #Extend FORs index so 1 degree intervals, then interpolate
    forsRegression = forsRegression.reindex(range(forsRegression.index.min(),forsRegression.index.max()+1))
    forsRegression.interpolate(inplace=True)
    #Want FORs as fractions, but given as percents in spreadsheet
    forsRegression /= 100 

    #Create dict mapping from fleet plant types to plant type codes in FOR sheet
    forPTMatching = {'Coal Steam':'ST','Coal Steam CCS':'ST','Biomass':'ST','O/G Steam':'ST','Geothermal':'ST',
        'Nuclear':'NU',
        'Landfill Gas':'CT','Combustion Turbine':'CT',
        'Combined Cycle':'CC','Combined Cycle CCS':'CC',
        'Fossil Waste':'Other','Non-Fossil Waste':'Other','Fuel Cell':'Other','Municipal Solid Waste':'Other','Solar PV':'Other','Onshore Wind':'Other',
        'Offshore Wind':'Other','Hydro':'HD'}
    return forsRegression,forPTMatching

def calculateFORs(temps,fleet,compressedGens,tVar,forsRegression,forPTMatching):    
    #For each generator w/ a TDFOR, replace temperatures w/ TDFOR
    fors = list()
    for c in fleet['GAMS Symbol']:      
        genRow = fleet.loc[fleet['GAMS Symbol']==c].squeeze()

        if 'COMBINED' not in genRow['GAMS Symbol']: #If not a combined unit, get FOR that unit
            # fors.append(calculateFORsForGen(temps,genRow,forsRegression,forPTMatching,tVar))#.values #use values due to datetime format mismatch error
            tempFors = calculateFORsForGen(temps,genRow,forsRegression,forPTMatching,tVar)#.values #use values due to datetime format mismatch error
        else: #If a combined unit, get FORs for constituent units and average
            constituentUnits = compressedGens.loc[compressedGens['UnitCompressedInto']==c]
            constituentFORs = [calculateFORsForGen(temps,constituentUnits.iloc[i],forsRegression,forPTMatching,tVar) for i in range(constituentUnits.shape[0])]
            tempFors = pd.concat(constituentFORs,axis=1).mean(axis=1)#.values #use values due to datetime format mismatch error
        tempFors.name = genRow['GAMS Symbol']
        fors.append(tempFors)

    fors = pd.concat(fors,axis=1)

    return fors 

#Calculate FORs for given generator using temps
def calculateFORsForGen(temps,genRow,forsRegression,forPTMatching,tVar):
    pt,lat,lon = genRow['PlantType'],genRow['Latitude'],genRow['Longitude']
    # print(genRow)

    #Get closest coords w/ non-NAN temperatures
    lat,lon = getClosestCellCoordsWithMet(temps,lat,lon,tVar)

    #Get temperatures & dt index from coordinate, then store in Series
    genTemps = temps.sel({'lat':lat,'lon':lon},method='nearest')
    genTemps,times = np.array(genTemps.variables[tVar][:]),genTemps['time'] #in time,lat,lon; swap so lat,lon,time
    genTemps = pd.Series(genTemps,index=times)

    #If CESM data, T given in K; convert to C and round to 5 degree intervals
    if np.nanmax(genTemps.values)>150: #screens for K values
        genTemps -= 273.15
        genTemps = genTemps.round() #5*((genTemps/5).round())
        genTemps[genTemps>35],genTemps[genTemps<-15] = 35,-15

    #Get FOR type
    forType = forPTMatching[pt] if pt in forPTMatching else 'Other'

    #Replace temperatures w/ FORs        
    tToFOR = forsRegression[forType].to_dict()
    genFORs = genTemps.replace(tToFOR)
    # genFORs.name = genRow['GAMS Symbol']
    return genFORs 

#Calculate FORs for new techs using location of techs
def calculateNewTechFORs(temps,newTechsCE,tVar,forsRegression,forPTMatching):
    forsTechs = list()
    for c in newTechsCE['GAMS Symbol']:      
        techRow = newTechsCE.loc[newTechsCE['GAMS Symbol']==c].squeeze()
        tempFors = calculateFORsForGen(temps,techRow,forsRegression,forPTMatching,tVar)#.values #use values due to datetime format mismatch error
        tempFors.name = c
        forsTechs.append(tempFors)

    forsTechs = pd.concat(forsTechs,axis=1)

    return forsTechs

# #Calculate FORs for new techs by taking average FORs of all existing generators
# #of same plant type and region.
# def calculateNewTechFORs(forsTechs,fors,genFleetForCE,newTechsCE):
#     for tech in forsTechs.columns:
#         #Get tech's plant type and region
#         pt = newTechsCE.loc[newTechsCE['GAMS Symbol'] == tech,'PlantType'].values[0]
#         region = newTechsCE.loc[newTechsCE['GAMS Symbol'] == tech,'region'].values[0]

#         #Skip wind & solar techs, since don't have TDFORs
#         if 'wind' not in tech and 'solar' not in tech:

#             #Get existing generators of same pt and region
#             gensOfPT = genFleetForCE.loc[genFleetForCE['PlantType']==pt]
#             gensOfPTAndRegion = gensOfPT.loc[gensOfPT['region']==region]

#             #If existing generators of same pt and region, get FORs and average them
#             if gensOfPTAndRegion.shape[0]>0:
#                 gensFORs = fors[gensOfPTAndRegion['GAMS Symbol']]
#                 forsTechs[tech] = gensFORs.mean(axis=1)
#     return forsTechs

