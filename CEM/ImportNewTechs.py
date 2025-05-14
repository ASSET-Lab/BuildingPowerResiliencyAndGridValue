#Michael Craig
#October 4, 2016
#Function imports data for new technologies eligible for construction in capacity expansion model

import os, pandas as pd
from SetupGeneratorFleet import addUnitCommitmentParameters,addFuelPrices,addRandomOpCostAdder,calcOpCost,addRegResCostAndElig,addReserveEligibility,convertCostToTgtYr,getATBTechDetailsForWindSolarBattery
from CalculateDerates import getClosestCellCoordsWithMet
from ImportNonREMetVars import importCESMMetVars,importTGWMetVars

def getNewTechs(regElig,regCostFrac,currYear,stoInCE,seasStoInCE,fuelPrices,yearIncDACS,yearIncHydrogen,yearIncCCS,yearIncNuclear,
                transRegions,contFlexInelig,weatherYears,climateScenario,genFleet,ira,lbToShortTon,onlyNSPSUnits=False,
                allowCoalWithoutCCS=False,newPlantDataDir=os.path.join('Data','NewPlantData'),
                iraITCTechs = ['Battery Storage','Hydrogen','Solar PV'],iraPTCTechs=['Onshore Wind','Offshore Wind','Nuclear'],
                ira45QTechs = ['Combined Cycle CCS']):
    print('in ImportNewTechs, using non-NSPS-compliant units so CTs are investment option!')
    if currYear > 2050: currYear = 2050
    #Read in new techs and add parameters
    newTechsCE = pd.read_excel(os.path.join(newPlantDataDir,'NewTechFramework.xlsx'))
    newTechsCE = extractATBDataForCurrentYear(newTechsCE,newPlantDataDir,currYear)
    newTechsCE = addHydrogenCostAndEffParams(newTechsCE)
    newTechsCE = addUnitCommitmentParameters(newTechsCE,'PhorumUCParameters.csv') 
    newTechsCE = addUnitCommitmentParameters(newTechsCE,'StorageUCParameters.csv')
    newTechsCE = addFuelPrices(newTechsCE,currYear,fuelPrices)
    if currYear >= yearIncDACS: newTechsCE = addDACS(newTechsCE,fuelPrices,currYear)
    newTechsCE = addRandomOpCostAdder(newTechsCE)
    newTechsCE = calcOpCost(newTechsCE)
    newTechsCE = addRegResCostAndElig(newTechsCE,regElig,regCostFrac)
    newTechsCE = addReserveEligibility(newTechsCE,contFlexInelig)
    #Modify costs per IRA (https://www.whitehouse.gov/cleanenergy/clean-energy-tax-provisions/)
    print('Not including IRA!')
    # if ira: 
    ## IF UNCOMMENT THEN PUT WIND WITH CAPEX SAVINGS
    #     newTechsCE.loc[newTechsCE['PlantType'].isin(iraITCTechs),'CAPEX($/MW)'] *= 0.7 #apply 30% ITC - assumes meet labor requirements
    #     newTechsCE.loc[newTechsCE['PlantType'].isin(iraPTCTechs),'OpCost($/MWh)'] -= 15 #apply 15 $/MWh PTC - assumes meet labor requirements (1.5 cents/kWh * 5 for apprenticeship)
    #     newTechsCE.loc[newTechsCE['PlantType'].isin(ira45QTechs),'OpCost($/MWh)'] -= (newTechsCE.loc[newTechsCE['PlantType'].isin(ira45QTechs),'CO2EmRate(lb/MMBtu)']
    #                                                         *newTechsCE.loc[newTechsCE['PlantType'].isin(ira45QTechs),'Heat Rate (Btu/kWh)']/1e6*1000/lbToShortTon*85) #apply $85/ton 45Q - assumes meet labor requirements
    #Discount costs
    for c,l in zip(['CAPEX($/MW)','FOM($/MW/yr)'],['occ','fom']): newTechsCE[c] = convertCostToTgtYr(l,newTechsCE[c])
    #Filter plants
    if allowCoalWithoutCCS == False: newTechsCE = newTechsCE.loc[newTechsCE['PlantType'] != 'Coal Steam']
    if onlyNSPSUnits: newTechsCE = newTechsCE.loc[newTechsCE['NSPSCompliant'] == 'Yes']
    if not stoInCE: newTechsCE = newTechsCE.loc[newTechsCE['FuelType'] != 'Energy Storage']
    if not seasStoInCE: newTechsCE = newTechsCE.loc[newTechsCE['PlantType'] != 'Hydrogen']
    if currYear < yearIncHydrogen: newTechsCE = newTechsCE.loc[newTechsCE['PlantType'] != 'Hydrogen']
    if currYear < yearIncCCS: newTechsCE = newTechsCE.loc[~newTechsCE['PlantType'].str.contains('CCS')]
    if currYear < yearIncNuclear: newTechsCE = newTechsCE.loc[newTechsCE['PlantType'] != 'Nuclear']
    #If not a climate change analysis, copy non-RE units for each region; otherwise, copy for each met var cell.
    # if climateScenario == None: newTechsCE = repeatNonRETechOptionsForEachRegion(newTechsCE,transRegions)
    # else: newTechsCE = repeatNonRETechOptionsForEachCell(newTechsCE,transRegions,genFleet,climateScenario,weatherYears)
    newTechsCE = repeatNonREAndNonThermalTechOptionsForEachRegion(newTechsCE,transRegions)
    newTechsCE = repeatThermalTechOptionsForEachCell(newTechsCE,transRegions,genFleet,climateScenario,weatherYears)
    newTechsCE.reset_index(inplace=True,drop=True)
    #Scale up capital costs for more than 1 weather year
    if len(weatherYears)>1: newTechsCE = scaleUpCapitalCostsForMultiYears(newTechsCE,weatherYears)
    return newTechsCE

#Get data from NREL's ATB for current year of analysis and input into new techs df
def extractATBDataForCurrentYear(newTechsCE,newPlantDataDir,currYear,scenario='Moderate'):
    print('Using VOM of 0.1$/MWh for storage to eliminate simultaneous charging & discharging')

    #Get tech details for wind, solar, battery in ATB
    solarTech,windTech,batteryTech,windTechDetail,solarTechDetail,batteryTechDetail,offwindTech,offwindTechDetail = getATBTechDetailsForWindSolarBattery()

    #Create mapping from our technologies to names in ATB
    ptToATBTechAlias = {'Solar PV':solarTech,'Onshore Wind':windTech,'Offshore Wind':offwindTech,'Battery Storage':batteryTech,'Nuclear':'Nuclear',
                    'Coal Steam CCS':'Coal','Combined Cycle':'Natural Gas','Combined Cycle CCS':'Natural Gas','Combustion Turbine':'Natural Gas'}
    feToATBTechDetail = {'Solar PV':solarTechDetail,'Onshore Wind':windTechDetail,'Offshore Wind':offwindTechDetail,'Battery Storage':batteryTechDetail,'Nuclear':'Nuclear',
                    'Coal Steam CCS':'CCS90AvgCF','Combined Cycle':'CCAvgCF','Combined Cycle CCS':'CCCCSAvgCF','Combustion Turbine':'CTAvgCF'}
    ptToATBHRs = {'Coal Steam CCS':'Coal-90%-CCS','Combined Cycle':'NG F-Frame CC','Combined Cycle CCS':'NG F-Frame CC 90% CCS','Combustion Turbine':'NG F-Frame CT','Nuclear':'Nuclear - AP1000'}

    #Get plant types of storage
    stoPTs = newTechsCE.loc[newTechsCE['ThermalOrRenewableOrStorage']=='storage']

    #Import ATB
    atb = pd.read_csv(os.path.join(newPlantDataDir,'ATBe.csv'),index_col=0,header=0)
    hrs = pd.read_csv(os.path.join(newPlantDataDir,'ATBHeatRatesJuly2022Edition.csv'),index_col=0,skiprows=1,header=0) #MMBtu/MWh

    #Filter ATB rows
    atb = atb.loc[atb['core_metric_variable']==currYear]
    atb = atb.loc[atb['core_metric_case']=='Market'] 
    atb = atb.loc[atb['scenario']==scenario]
    hrs = hrs.loc[hrs['scenario']==scenario]
    hrs = hrs[str(currYear)]

    #Add values for each plant type
    for pt in newTechsCE['PlantType']:
        if pt in ptToATBTechAlias:
            #Get rows for tech type
            techRows = atb.loc[atb['technology_alias']==ptToATBTechAlias[pt]]
            techRows = techRows.loc[techRows['techdetail']==feToATBTechDetail[pt]]

            #Extract parameters (use .iloc[0] because of redundancies in data that do not change values of our parameters of interest)
            cap = techRows.loc[techRows['core_metric_parameter']=='CAPEX']['value'].iloc[0] #$/kW
            fom = techRows.loc[techRows['core_metric_parameter']=='Fixed O&M']['value'].iloc[0] #$/kW-yr
            if 'Variable O&M' in techRows['core_metric_parameter'].unique():
                vom = techRows.loc[techRows['core_metric_parameter']=='Variable O&M']['value'].iloc[0]
            elif pt in stoPTs:
                vom = .1 #$/MWh
            else:
                vom = 0 #$/MWh
            # vom = techRows.loc[techRows['core_metric_parameter']=='Variable O&M']['value'].iloc[0] if 'Variable O&M' in techRows['core_metric_parameter'].unique() else 0 #$/MWh

            #Get heat rate from separate file
            hr = hrs[ptToATBHRs[pt]] if pt in ptToATBHRs else 0

            #Add parameters
            newTechsCE.loc[newTechsCE['PlantType']==pt,'CAPEX($/MW)'] = cap*1000 #$/kW to $/MW
            newTechsCE.loc[newTechsCE['PlantType']==pt,'FOM($/MW/yr)'] = fom*1000 #$/kW-yr to $/MW-yr
            newTechsCE.loc[newTechsCE['PlantType']==pt,'VOM($/MWh)'] = vom #already in $/MWh
            newTechsCE.loc[newTechsCE['PlantType']==pt,'Heat Rate (Btu/kWh)'] = hr*1000
    return newTechsCE

#Hydrogen not in ATB, so add parameters for current year of analysis and input into new techs df.
#Use separate capex values for power (pcapex) and energy (ecapex).
#Units: CAPEX($/MW), FOM($/MW/yr), VOM($/MWh), Efficiency, ECAPEX($/MWH).
#Source: Dowling, Joule, 2020, Role of long-duration energy storage in variable renewable electricity systems.
def addHydrogenCostAndEffParams(newTechsCE,pcapex=1058000,fom=0,vom=0,eff=.49,ecapex=160):
    newTechsCE.loc[newTechsCE['PlantType'] == 'Hydrogen','CAPEX($/MW)'] = pcapex
    newTechsCE.loc[newTechsCE['PlantType'] == 'Hydrogen','FOM($/MW/yr)'] = fom
    newTechsCE.loc[newTechsCE['PlantType'] == 'Hydrogen','VOM($/MWh)'] = vom
    newTechsCE.loc[newTechsCE['PlantType'] == 'Hydrogen','Efficiency'] = eff
    newTechsCE.loc[newTechsCE['PlantType'] == 'Hydrogen','ECAPEX($/MWH)'] = ecapex
    return newTechsCE

def scaleUpCapitalCostsForMultiYears(newTechsCE,weatherYears):
    numYearsInCE = len(weatherYears)
    newTechsCE['CAPEX($/MW)'] *= numYearsInCE
    newTechsCE['FOM($/MW/yr)'] *= numYearsInCE
    newTechsCE['ECAPEX($/MWH)'] *= numYearsInCE
    return newTechsCE

#Add DACS parameters to new techs. DACS gen & capac values are negative in optim, so 
#all cost values (op, fom, cap) are negative (so that gen or cap * cost = positive) and emission rate is 
#positive (so that gen * ER = negative). Data based on Keith's Joule paper.
def addDACS(newTechsCE,fuelPrices,currYear):
    dacsCap = -500
    #CO2 ems removal rate NET of emissions from NG for heat. 
    #Keith: 1 t CO2 removal requires 366 kWh e [burned 5.25 GJ NG is also captured; not included in 1 t co2]
    #1 t co2/366 kwh * 2000 lb/t * 1000 kwh/mwh * mwh/3.412 mmbtu = 1601.5 lb CO2/MMBtu
    #If include heat, 1 t co2/(366 kwh + 5.25gj * .277 mwh/gj * 1000 kwh/mwh) = 1 t co2/1820.25 kwh
    #               1 t co2/1820.25 kwh * 2000lb/t * 1000kwh/mwh * mwh/3.412 mmbtu = 322 lb co2/mmbtu
    #Realmonte: 1 t CO2 removed per 500 kWh e (given as 1.8 GJ electricity, 1 GJ = 0.27778 MWh) (also 8.1 GJ heat)
    #1 t co2/1.8 gj * gj/.277 mwh * mwh/1000 kwh * 2000 lb/t * 1000 kwh/mwh * mwh/3.412 mmbtu = 1172.3 lb co2/mmbtu
    dacsNetEmsRate = 322 #lb CO2/MMBtu ; 1601.5 for keith, 1172.3 for realmonte, 322 for keith w/ NG heat
    dacsHR = 3412 #btu/kWh; just conversion factor
    #Set op costs as VOM because HR is not an accurate value
    vom = 26/366*1000 #$/MWh; Keith gives $26/t CO2; instead use $26/366 kWh
    fuelPrices = fuelPrices.loc[currYear] if currYear in fuelPrices.index else fuelPrices.iloc[-1]
    fuelPrices = convertCostToTgtYr('fuel',fuelPrices)
    ngPrice = fuelPrices.to_dict()['Natural Gas']   
    natGasCost = ngPrice * 5.25/366 * 1000 * 0.947 #$/MWh; 5.25 GJ NG/366 kwh given in Keith; * ng price ($/mmbtu) * conversions
    totalOpCost = natGasCost + vom
    capCost = 779.5*1e6/40.945 #2086000000/40.945 # 779.5*1e6/40.945 #$779.5M buys 0.98Mt co2/yr; @ 366kwh/1 t co2, that is 0.98*366*1e6/8760 = 40.945 MW
    #Add row to new techs df
    newRow = {'PlantType':['DAC'],'DataSource':['handCalc'],'FuelType':['DAC'],'Capacity (MW)':[dacsCap],
        'Heat Rate (Btu/kWh)':[dacsHR],'CAPEX($/MW)':[-capCost],'FOM($/MW/yr)':[0],'VOM($/MWh)':[-totalOpCost],
        'NSPSCompliant':['Yes'],'CO2EmRate(lb/MMBtu)':[dacsNetEmsRate],'Lifetime(years)':[30],
        'FuelPrice($/MMBtu)':[0],'RampRate(MW/hr)':[abs(dacsCap)],'MinLoad(MWh)':[0],'MinDownTime(hrs)':[0],'StartCost($)':[0]}
    newTechsCE = pd.concat([newTechsCE,pd.DataFrame(newRow)])
    newTechsCE.reset_index(drop=True,inplace=True)
    return newTechsCE

#For each non-wind & non-solar tech option, repeat per region. (New W&S use lat/long coords later.)
def repeatNonREAndNonThermalTechOptionsForEachRegion(newTechsCE,transRegions):
    # newTechsRE = newTechsCE.loc[newTechsCE['PlantType'].isin(['Onshore Wind','Offshore Wind','Solar PV'])].copy()
    # newTechsNotRE = newTechsCE.loc[~newTechsCE.index.isin(newTechsRE.index)].copy()
    newREAndThermalTechs = newTechsCE.loc[newTechsCE['ThermalOrRenewableOrStorage'].isin(['thermal','renewable'])].copy()
    newTechsPerRegion = newTechsCE.loc[~newTechsCE.index.isin(newREAndThermalTechs.index)].copy()
    l = [newREAndThermalTechs]
    for r in transRegions:
        regionTechs = newTechsPerRegion.copy()
        regionTechs['region'] = r
        l.append(regionTechs)
    return pd.concat(l)

#Add non-RE tech investment option for each cell that has an existing large thermal unit by p-region
def repeatThermalTechOptionsForEachCell(newTechsCE,transRegions,genFleet,climateScenario,weatherYears,existingPTs=['Coal Steam','Combined Cycle','Nuclear']):
    #Separate RE from non-RE
    newThermalTechs = newTechsCE.loc[newTechsCE['ThermalOrRenewableOrStorage']=='thermal'].copy()
    otherTechs = newTechsCE.loc[~newTechsCE.index.isin(newThermalTechs.index)].copy()
    #Import temperatures
    metVars = importTGWMetVars([weatherYears[0]],climateScenario) #only need 1 weather year b/c just need the set of lat/longs
    #Loop through regions
    l = [otherTechs]
    for r in transRegions:
        #Get existing gens (retirement status does not matter) of right PT in region
        gensRegionPT = genFleet.loc[genFleet['PlantType'].isin(existingPTs)]
        gensRegionPT = gensRegionPT.loc[gensRegionPT['region']==r]

        cellsWithGens = list()
        for i,g in gensRegionPT.iterrows():
            #Get closest cell coords w/ non-NAN temperatures
            lat,lon = g['Latitude'],g['Longitude']
            lat,lon = getClosestCellCoordsWithMet(metVars,lat,lon)
            cell = metVars.sel({'lat':lat,'lon':lon},method='nearest') 
            #Extract cell's lat & lon & add to list
            lat,lon = float(cell['lat'].values),float(cell['lon'].values)
            cellsWithGens.append((lat,lon))
        #Get unique cells with generators in region
        cellsWithGens = list(set(cellsWithGens))
        #Add a thermal gen copy at each cell w/in region
        for c in cellsWithGens:
            cellTechs = newThermalTechs.copy()
            cellTechs['region'],cellTechs['Latitude'],cellTechs['Longitude'] = r,c[0],c[1]
            l.append(cellTechs)
    return pd.concat(l)

########## OLD ######################
#Account for ITC in RE cap costs
#http://programs.dsireusa.org/system/program/detail/658
def modRECapCostForITC(newTechsCE,currYear):
    if currYear > 2050: currYear = 2050
    windItc,windItcYear = .21,2020 #wind ITC expires at 2020; .21 is average of 2016-2019 ITCs
    solarItcInit,solarItcIndef, solarItcYear = .3,.1,2020 #solar ITC doesn't expire, but goes from .3 to .1
    if currYear <= windItcYear: 
        modRECost(newTechsCE,windItc,'Onshore Wind')
        modRECost(newTechsCE,windItc,'Offshore Wind')
    if currYear <= solarItcYear: modRECost(newTechsCE,solarItcInit,'Solar PV')
    else: modRECost(newTechsCE,solarItcIndef,'Solar PV')
    
def modRECost(newTechsCE,itc,plantType):
    ptCol = newTechsCE[0].index('PlantType')
    capexCol = newTechsCE[0].index('CAPEX($/MW)')
    ptRow = [row[0] for row in newTechsCE].index(plantType)
    newTechsCE[ptRow][capexCol] *= (1-itc)
