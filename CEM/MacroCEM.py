import sys, os, csv, operator, copy, time, random, warnings, numpy as np, datetime as dt, pandas as pd
from os import path; from gams import *
from SetupGeneratorFleet import setupGeneratorFleet,compressAndAddSizeDependentParams
from AddCoolingTypes import addCoolingTypes
from ProcessHydro import processHydro
from UpdateFuelPriceFuncs import updateFuelPricesAndCosts
from ImportDemand import importDemand,createTimeBlocks
from DemandFuncsCE import getHoursForCE
from IsolateDataForCE import isolateDataInCEHours,isolateDataInCEBlocks
from ImportNewTechs import getNewTechs
from RetireUnitsCFPriorCE import retireUnitsCFPriorCE
from CreateFleetForCELoop import createFleetForCurrentCELoop
from GetRenewableCFs import getREGen
from GetNewRenewableCFs import getNewRenewableCFs
from AddWSSitesToNewTechs import addWSSitesToNewTechs
from ProcessCEResults import saveCEBuilds,addNewGensToFleet,addNewLineCapToLimits,retireUnitsByCF
from ScaleRegResForAddedWind import scaleRegResForAddedWind
from CombinePlants import combineWindSolarStoPlants
from GAMSAddSetToDatabaseFuncs import *
from GAMSAddParamToDatabaseFuncs import *
from InitializeOnOffExistingGensCE import initializeOnOffExistingGens
from ReservesWWSIS import calcWWSISReserves
from GetIncResForAddedRE import getIncResForAddedRE
from SaveCEOperationalResults import saveCapacExpOperationalData
from WriteTimeDependentConstraints import writeTimeDependentConstraints
from WriteBuildVariable import writeBuildVariable
from WriteBuildingUpgradeFix import writeBuildingUpgradeFix
from CreateEmptyReserveDfs import createEmptyReserveDfs
from SetupTransmissionAndZones import setupTransmissionAndZones, defineTransmissionRegions
from DefineReserveParameters import defineReserveParameters
from CalculateDerates import calculateLineThermalDerates,calculatePlantCapacityDerates
from ImportPRMCapacityAdjustments import importPRMCapacityAdjustments
from AddRegionalGasPricesFromJill import addRegionalGasPricesFromJill
from GAMSAuxFuncs import extract1dVarResultsFromGAMSModel,extract2dVarResultsIntoDictNoLA

# SET OPTIONS
warnings.filterwarnings("ignore")
pd.set_option('display.max_rows', 10)
pd.set_option('display.max_columns', 10)

# SCALARS
mwToGW = 1000
lbToShortTon = 2000

# Whether running on SC
RUNONSC = True

# ##############################################################################
# ##### UNIVERSAL PARAMETERS ###################################################
# ##############################################################################
def setKeyParameters():
    # ### START YEAR, END YEAR, AND STEPS FOR CE
    startYear, endYear, yearStepCE = 2045,2051, 5 #5

    # ### RE UPSAMPLING
    reDownFactor = .3                     # FRACTION of wind & solar sites per region dropped (0 = no sites dropped; .7 = 70% sites dropped); sites dropped with worst CFs
    
    # ### BUILD LIMITS
    yearIncDACS,yearIncHydrogen,yearIncCCS,yearIncNuclear = 2051,2051,2040,2050#2041,2041,2031,2031 #2040,2055,2035,2035  #year to allow investment in certain technologies

    # ### CE OPTIONS
    numBlocks, daysPerBlock, daysPerPeak = 4,7,3
    removeHydro = True                                  #whether to remove hydropower from fleet & subtract generation from demand, or to include hydro as dispatchable in CE w/ gen limit
    stoInCE,seasStoInCE = True,False                    # whether to allow new storage,new seasonal storage in CE model
    windGenFracOfDemand=0 #RPS wind carve-out
    prm=13 #PRM in excess of peak demand (%)

    return (startYear,endYear,yearStepCE,reDownFactor,yearIncDACS,yearIncHydrogen,yearIncCCS,yearIncNuclear,
        numBlocks,daysPerBlock,daysPerPeak,removeHydro,stoInCE,seasStoInCE,windGenFracOfDemand,prm)

def stateAssumptions(region,yearStepCE,reBuildRateMultiplier=1,thermalBuildRateMultiplier=4): 
    # ### MAX BUILDS
    #Max builds (MW) per region or, for wind & solar, grid point
    areaPerLatLongBox = 9745 #km^2 per degree lat x long (https://www.usgs.gov/faqs/how-much-distance-does-a-degree-minute-and-second-cover-your-maps?qt-news_science_products=0#qt-news_science_products)
    windDens,solarDens = .9,5.7 #W/m^2 (equiv to MW/km^2); https://www.seas.harvard.edu/news/2018/10/large-scale-wind-power-would-require-more-land-and-cause-more-environmental-impact
    reDensity = {'Onshore Wind':windDens,'Solar':solarDens,'Offshore Wind':windDens}
    #Set max capacity per tech. To turn off storage, use (stoInCE,seasStoInCE) flags above. 
    maxCapPerTech = {'Onshore Wind': areaPerLatLongBox * windDens, 'Offshore Wind': areaPerLatLongBox * windDens, 'Solar': areaPerLatLongBox * solarDens, 
        'Combined Cycle': 3000*thermalBuildRateMultiplier*yearStepCE,'Storage': 0, 'Dac': -0, 'CCS': 3000*thermalBuildRateMultiplier*yearStepCE, 
        'Nuclear': 3000*thermalBuildRateMultiplier*yearStepCE, 'Battery Storage': 99999, 'Hydrogen': 99999, 'Transmission': 99999} #across WECC, 2400 MW from 2019-2021 & 4100 MW from 2015-2021 of new NGCC + NGCT. 2020 annual max (1500 MW)
    #Max wind & solar builds (MW) per interconnection & region
    reCAGR = .3 #compounded annual growth rate for max RE builds - applied to histMaxREBuildRegion
    if region == 'WECC': #see SetupTransmissionAndZones for list of zones
        maxREPerZone = {'Onshore Wind': {'NWPP_NE':99999,'CAMX':99999,'Desert_Southwest':99999,'NWPP_Central':99999,'NWPP_NW':99999},
                        'Solar':{'NWPP_NE':99999,'CAMX':99999,'Desert_Southwest':99999,'NWPP_Central':99999,'NWPP_NW':99999}} 
        histMaxREBuildRegion = {'Wind':6700/2*yearStepCE*reBuildRateMultiplier,'Solar':7500/1*yearStepCE*reBuildRateMultiplier} #max build per interconnection per CE run. WECC 2022 WARA: 7.5 GW solar expected to be added in 2023; EIA 860: 6.7 GW & 5.8 GW wind & solar added in 2020 & 2021 across WECC
    elif region == 'EI':
        maxREPerZone = {'Onshore Wind': {'SERC':99999,'NY':99999,'NE':99999,'MISO':99999,'PJM':99999,'SPP':99999},
                        'Solar':{'SERC':99999,'NY':99999,'NE':99999,'MISO':99999,'PJM':99999,'SPP':99999}} 
        histMaxREBuildRegion = {'Wind':21000/3*reBuildRateMultiplier,'Solar':15000/3*reBuildRateMultiplier} #max build per interconnection per CE run. EI using 860 data: 15 GW solar 2019-2021; 21 GW 2019-2021 
    elif region == 'NY':
        maxREPerZone = {'Onshore Wind': {'p127':99999,'p128':99999},'Solar':{'p127':99999,'p128':99999},'Offshore Wind': {'p127':99999,'p128':99999}} 
        # histMaxREBuildRegion = {'Onshore Wind':567,'Solar':452,'Offshore Wind':567} #assume OSW = ONW. https://www.eia.gov/electricity/data.php#gencapacity, Existing capacity by energy source, by producer, by state back to 2000 (annual data from the EIA-860) 
        histMaxREBuildRegion = {'Onshore Wind':99999,'Solar':99999,'Offshore Wind':99999} #assume OSW = ONW. https://www.eia.gov/electricity/data.php#gencapacity, Existing capacity by energy source, by producer, by state back to 2000 (annual data from the EIA-860) 
        print('Running without limits on max builds')

    # ### CO2 EMISSION CAPS [https://www.eia.gov/environment/emissions/state/, table 3]
    if region == 'ERCOT': co2EmsInitial =  130594820     #METRIC TONS. Initial emission for ERCOT: 130594820.
    elif region == 'EI': 
        co2EmsInitial,wsGenFracOfDemandInitial,wGenFracOfDemandInitial =  1043526617,0,0 #1043526617 is from CEM output (2023 w/out investments); 1274060000 is from EIA data
        print('if running EI w/ renewable generation requirements, need to populate wsGenFracOfDemandInitial & wGenFracOfDemandInitial')
    elif region == 'WECC': co2EmsInitial,wsGenFracOfDemandInitial,wGenFracOfDemandInitial =  248800000,7,2   #2019 emissions: 248800000 METRIC TONS. wa,or,ca,nm,az,nv,ut,co,wy,id,mt
    elif region == 'NY': co2EmsInitial,wsGenFracOfDemandInitial,wGenFracOfDemandInitial =  30788000,5,3.6 #https://www.eia.gov/electricity/annual/html/epa_09_05.html - 2022, Table 9.5. Emissions from Energy Consumption; https://www.eia.gov/electricity/data.php#gencapacity, Detailed preliminary EIA-923 monthly and annual survey data (back to 1990)

    # ### CE AND UCED/ED OPTIONS
    balAuths = 'full'                                   # full: run for all BAs in interconn. TODO: add selection of a subset of BAs. [10/31 MC note: Papa Yaw has this code]
    compressFleet = True                                                # whether to compress fleet
    tzAnalysis = {'ERCOT':'CST','EI':'EST','WECC':'PST','NY':'EST'}[region]     # timezone for analysis
    fuelPrices = importFuelPrices('Reference case')                     # import fuel price time series
    transmissionEff = 0.95                                              # efficiency of transmission between zones (https://ars.els-cdn.com/content/image/1-s2.0-S2542435120305572-mmc1.pdf)
    greenField = False      #whether to run CE greenfield (no existing units) or not

    # ### CE OPTIONS
    runCE,ceOps = True,'ED'                           # ops are 'ED' or 'UC' (econ disp or unit comm constraints)
    includeRes = False                                  # whether to include reserves in CE & dispatch models (if False, multiplies reserve timeseries by 0)
    retireByAge = True                                  # whether to retire by age or not
    prmBasis = 'demand'                      # whether basis for planning reserve margin is peak demand ('demand') or net demand ('netdemand')

    retirementCFCutoff = .3                             # retire units w/ CF lower than given value
    ptEligRetCF = ['Coal Steam','O/G Steam','Combined Cycle'] #['Coal Steam']                        # which plant types retire based on capacity factor (economics)
    discountRate = 0.07 

    # ### DEMAND FLEXIBILITY PARAMETERS
    demandShifter = 0                                   # Percentage of hourly demand that can be shifted
    demandShiftingBlock = 4                             # moving shifting demand window (hours)

    # ### WARNINGS OR ERRORS
    if ceOps == 'UC': sys.exit('CEwithUC.gms needs to be updated for DACS operations - add DACS constraints and include gentechs set')

    return (balAuths,co2EmsInitial,compressFleet,tzAnalysis,fuelPrices,transmissionEff,
        runCE,ceOps,includeRes,retireByAge,prmBasis,retirementCFCutoff,ptEligRetCF,discountRate,maxCapPerTech,
        maxREPerZone,histMaxREBuildRegion,reCAGR,wsGenFracOfDemandInitial,wGenFracOfDemandInitial,reDensity,greenField,demandShifter,demandShiftingBlock)

def storageAssumptions():
    stoMkts = 'energy'                            # energy,res,energyAndRes - whether storage participates in energy, reserve, or energy and reserve markets
    stoFTLabels = ['Energy Storage','Pumped Storage']
    stoDuration = {'Energy Storage':'st','Hydrogen':'lt','Battery Storage':'st','Flywheels':'st','Batteries':'st','Pumped Storage':'st'} # mapping plant types to short-term (st) or long-term (lt) storage
    stoPTLabels = [pt for pt in stoDuration ]
    initSOCFraction = {pt:{'st':0,'lt':0}[dur] for pt,dur in stoDuration.items()} # get initial SOC fraction per st or lt storage
    stoMinSOC = 0     # min SOC
    stoEff = 0.81     #round-trip storage efficiency
    return stoMkts,stoFTLabels,stoDuration,stoPTLabels,initSOCFraction,stoMinSOC,stoEff

def importFuelPrices(fuelPriceScenario):
    fuelPrices = pd.read_csv(os.path.join('Data', 'Energy_Prices_Electric_Power.csv'), skiprows=4, index_col=0)
    fuelPrices = fuelPrices[[col for col in fuelPrices if fuelPriceScenario in col]]
    fuelPrices.columns = [col.split(':')[0] for col in fuelPrices.columns]
    return fuelPrices    
# ###############################################################################
# ###############################################################################
# ###############################################################################

# ###############################################################################
# ##### MASTER FUNCTION #########################################################
# ###############################################################################
#Main function to call. 
def macroCEM(region,climateScenario,wsGenFracOfDemand,wsGenYr,co2EndPercent,co2CapYr,forceBldgUpgrade,oswCapReq=9000,oswCapReqYr=2035):
    # Set key parameters
    (startYear,endYear,yearStepCE,reDownFactor,yearIncDACS,yearIncHydrogen,yearIncCCS,yearIncNuclear,
        numBlocks,daysPerBlock,daysPerPeak,removeHydro,stoInCE,seasStoInCE,windGenFracOfDemand,prm) = setKeyParameters()
    climateChange,nonCCReanalysis = True,False

    # Set assumptions
    (balAuths,co2EmsInitial,compressFleet,tzAnalysis,fuelPrices,transmissionEff,
        runCE,ceOps,includeRes,retireByAge,prmBasis,retirementCFCutoff,ptEligRetCF,discountRate,
        maxCapPerTechOrig,maxREPerZone,histMaxREBuildRegion,reCAGR,wsGenFracOfDemandInitial,wGenFracOfDemandInitial,
        reDensity,greenField,demandShifter,demandShiftingBlock) = stateAssumptions(region,yearStepCE)
    stoMkts,stoFTLabels,stoDuration,stoPTLabels,initSOCFraction,stoMinSOC,stoEff = storageAssumptions()
    (regLoadFrac, contLoadFrac, regErrorPercentile, flexErrorPercentile, regElig, contFlexInelig, regCostFrac,
        rrToRegTime, rrToFlexTime, rrToContTime) = defineReserveParameters(stoMkts, stoFTLabels)

    # Create results directory
    resultsDirAll = region+climateScenario+'RPS'+str(wsGenFracOfDemand)+'RY'+str(wsGenYr)+'CO2'+str(co2EndPercent)+'CY'+str(co2CapYr)+'Up'+str(forceBldgUpgrade)
    if not os.path.exists(resultsDirAll): os.makedirs(resultsDirAll)
    pd.Series(co2EmsInitial).to_csv(os.path.join(resultsDirAll,'initialCO2Ems.csv'))

    # Setup initial fleet and demand
    (genFleet, compressedGens, transRegions, pRegionShapes, lineLimits, lineDists, lineCosts) = getInitialFleetAndTransmission(startYear, 
        fuelPrices, compressFleet, resultsDirAll, regElig, regCostFrac, stoMinSOC, greenField, region, balAuths, 
        contFlexInelig, stoFTLabels, stoPTLabels, stoEff, stoInCE, climateScenario)

    # Run CE and/or ED/UCED
    for currYear in range(startYear, endYear, yearStepCE):
        # Get weather years
        weatherYears = list(range(currYear-yearStepCE+1,currYear+yearStepCE+1))
        if climateScenario == 'historical': weatherYears = [i-40 for i in weatherYears]
        print('Weather years:',weatherYears)

        # Set CO2 cap
        if currYear < co2CapYr:
            # currCo2Cap = co2EmsInitial + (co2EndPercent/100*co2EmsInitial - co2EmsInitial)/((co2CapYr-1) - startYear) * (currYear - startYear)
            currCo2Cap = co2EmsInitial + (co2EndPercent/100*co2EmsInitial - co2EmsInitial)/(co2CapYr - startYear) * (currYear - startYear)
        else:
            currCo2Cap = co2EndPercent/100*co2EmsInitial
        
        #Don't need this for Jill runs because the demand data she gave me has a combined weight of 1 across all years!
        # currCo2Cap *= len(weatherYears) #adjust annual cap if running more than 1 year!

        # Set RPS
        if currYear < wsGenYr:
            # currWSGenFracOfDemand = wsGenFracOfDemandInitial + wsGenFracOfDemand/(wsGenYr-1 - startYear) * (currYear - startYear)
            currWSGenFracOfDemand = wsGenFracOfDemand/(wsGenYr - startYear) * (currYear - startYear)
        else:
            currWSGenFracOfDemand = wsGenFracOfDemand
        currWindGenFracOfDemand = wGenFracOfDemandInitial + windGenFracOfDemand/(endYear-1 - startYear) * (currYear - startYear)
        
        # Set offshore wind requirement
        if currYear < oswCapReqYr:
            currOSWCapReq = oswCapReq/(oswCapReqYr - startYear) * (currYear - startYear)
        else:
            currOSWCapReq = oswCapReq

        print('Entering year ', currYear, ' with CO2 cap (million tons):', round(currCo2Cap/1e6),'\t and RE (%) & OSW (MW) requirement:',round(currWSGenFracOfDemand),round(currOSWCapReq))

        # Set maximum RE additions using annual CAGR from initial max historic builds while accounting for # years included in CE
        maxREInRegion = dict()
        for re in histMaxREBuildRegion: maxREInRegion[re] = sum([histMaxREBuildRegion[re] * (1 + reCAGR) ** (currYear - yr - startYear) for yr in range(yearStepCE)])

        # Create results directory
        resultsDir = os.path.join(resultsDirAll,str(currYear))
        if not os.path.exists(resultsDir): os.makedirs(resultsDir)
        
        # Get electricity demand profile
        demand,demandWeights,demandSavings,demandCosts,demandYearUpgrade0 = importDemand(currYear,weatherYears,climateScenario)
        demand.to_csv(os.path.join(resultsDir,'demandInitial'+str(currYear)+'.csv'))
        demandYearUpgrade0.to_csv(os.path.join(resultsDir,'demandFullYearUpgrade0'+str(currYear)+'.csv'))
        demandWeights.to_csv(os.path.join(resultsDir,'demandWeights'+str(currYear)+'.csv'))
        demandCosts.to_csv(os.path.join(resultsDir,'demandCosts'+str(currYear)+'.csv'))
        for upgrade in demandSavings: demandSavings[upgrade].to_csv(os.path.join(resultsDir,'demandSavingsUpgrade'+str(upgrade)+str(currYear)+'.csv'))

        # Get intertie capacities and prices from outside region
        interties = pd.read_csv(os.path.join('Data','NYInterties.csv'))
        interties.to_csv(os.path.join(resultsDir,'interties.csv'))

        # Run CE
        if currYear > startYear and runCE:
            #Initialize results & inputs
            if currYear == startYear + yearStepCE: priorCEModel, priorHoursCE, genFleetPriorCE = None, None, None
            (genFleet, genFleetPriorCE, lineLimits,priorCEModel, priorHoursCE) = runCapacityExpansion(genFleet, demand, demandYearUpgrade0, demandWeights,demandSavings,demandCosts, currYear, weatherYears, prm, prmBasis,
                                        discountRate, fuelPrices, currCo2Cap, numBlocks, daysPerBlock, daysPerPeak,
                                        retirementCFCutoff, retireByAge, tzAnalysis, resultsDir,
                                        maxCapPerTechOrig, regLoadFrac, contLoadFrac, regErrorPercentile, flexErrorPercentile,
                                        rrToRegTime, rrToFlexTime, rrToContTime, regElig, regCostFrac, ptEligRetCF,
                                        genFleetPriorCE, priorCEModel, priorHoursCE, stoInCE, seasStoInCE,
                                        ceOps, stoMkts, initSOCFraction, includeRes, reDownFactor, demandShifter,
                                        demandShiftingBlock, region, yearIncDACS, yearIncHydrogen,yearIncCCS, yearIncNuclear, transRegions,pRegionShapes,
                                        lineLimits, lineDists, lineCosts, contFlexInelig, 
                                        nonCCReanalysis, stoFTLabels, transmissionEff, removeHydro, climateChange, climateScenario,
                                        maxREPerZone, maxREInRegion, compressedGens, currWSGenFracOfDemand, currWindGenFracOfDemand, currOSWCapReq, reDensity,
                                        co2EmsInitial, resultsDirAll, interties, forceBldgUpgrade)
# ###############################################################################
# ###############################################################################
# ###############################################################################

# ###############################################################################
# ###### SET UP INITIAL FLEET AND DEMAND ########################################
# ###############################################################################
def getInitialFleetAndTransmission(startYear, fuelPrices, compressFleet, resultsDir, regElig, regCostFrac, 
        stoMinSOC, greenField, region, balAuths, contFlexInelig, stoFTLabels, stoPTLabels, stoEff, stoInCE, climateScenario):
    # GENERATORS
    genFleet = setupGeneratorFleet(region, startYear, fuelPrices, stoEff, stoMinSOC, stoFTLabels, stoInCE)
    #Modify retirement year for Diablo Canyon given recent extension
    if 'Diablo Canyon' in genFleet['Plant Name'].unique():
        print('Extending Diablo Canyon lifetime to 2035')
        genFleet.loc[genFleet['Plant Name']=='Diablo Canyon','Retirement Year'] = 2035

    # ADD COOLING TYPES TO GENERATORS TO CAPTURE DERATINGS WHEN RUNNING CLIMATE SCENARIOS
    genFleet = addCoolingTypes(genFleet, region)

    # DEFINE TRANSMISSION REGIONS
    transRegions = defineTransmissionRegions(region, balAuths)

    # TRANSMISSION
    genFleet, transRegions, limits, dists, costs, pRegionShapes = setupTransmissionAndZones(genFleet, transRegions, region, balAuths)
    for df, l in zip([limits, dists, costs],['Limits', 'Dists', 'Costs']): df.to_csv(os.path.join(resultsDir, 'transmission' + l + 'Initial.csv'))

    # ADD REGION-SPECIFIC GAS PRICES FROM JILL
    genFleet = addRegionalGasPricesFromJill(genFleet)
    genFleet.to_csv(os.path.join(resultsDir, 'genFleetInitialPreCompression.csv'))

    # COMBINE GENERATORS FOR SMALLER GEN FLEET AND ADD SIZE DEPENDENT PARAMS (COST, REG OFFERS, UC PARAMS)
    genFleet,compressedGens = compressAndAddSizeDependentParams(genFleet, compressFleet, regElig, contFlexInelig, regCostFrac, stoFTLabels, stoPTLabels)

    # IF GREENFIELD, ELIMINATE EXISTING GENERATORS EXCEPT TINY NG, WIND, & SOLAR PLANT (TO AVOID CRASH IN LATER FUNCTIONS)
    if greenField: genFleet = stripDownGenFleet(genFleet, greenField)

    # IF RUNNING CESM, CONVERT LONGITUDE TO 0-360 FROM -180-180
    # if climateScenario != None: genFleet['Longitude'] = genFleet['Longitude']%360 #do this for TGW and CESM

    # SAVE FILES
    genFleet.to_csv(os.path.join(resultsDir,'genFleetInitial.csv')),compressedGens.to_csv(os.path.join(resultsDir,'compressedUnitsFromGenFleet.csv'))
    return genFleet, compressedGens, transRegions, pRegionShapes, limits, dists, costs
# ###############################################################################
# ###############################################################################
# ###############################################################################

# ###############################################################################
# ###### RUN CAPACITY EXPANSION #################################################
# ###############################################################################
def runCapacityExpansion(genFleet, demand, demandYearUpgrade0, demandWeights,demandSavings,demandCosts, currYear, weatherYears, prm, prmBasis, discountRate, fuelPrices, currCo2Cap, numBlocks,
                         daysPerBlock, daysPerPeak, retirementCFCutoff, retireByAge, tzAnalysis, resultsDirOrig, maxCapPerTechOrig,
                         regLoadFrac,contLoadFrac, regErrorPercentile, flexErrorPercentile, rrToRegTime, rrToFlexTime,  rrToContTime,
                         regElig, regCostFrac, ptEligRetCF, genFleetPriorCE, priorCEModel, priorHoursCE, stoInCE, seasStoInCE,
                         ceOps, stoMkts, initSOCFraction, includeRes, reDownFactor, demandShifter, demandShiftingBlock, 
                         region, yearIncDACS, yearIncHydrogen,yearIncCCS, yearIncNuclear, transRegions, pRegionShapes, lineLimits, lineDists, lineCosts, contFlexInelig,
                         nonCCReanalysis, stoFTLabels, transmissionEff, removeHydro, climateChange, climateScenario,
                         maxREPerZone, maxREInRegion, compressedGens, currWSGenFracOfDemand, currWindGenFracOfDemand, currOSWCapReq, reDensity, co2EmsInitial,
                         resultsDirAll, interties, forceBldgUpgrade):
    # ###############CREATE RESULTS DIRECTORY FOR CE RUN AND SAVE INITIAL INPUTS
    resultsDir = os.path.join(resultsDirOrig, 'CE')
    if not os.path.exists(resultsDir): os.makedirs(resultsDir)
    print('Entering CE loop for year ' + str(currYear))
    lineLimits.to_csv(os.path.join(resultsDir,'lineLimitsForCE' + str(currYear) + '.csv'))
    pd.Series(currCo2Cap).to_csv(os.path.join(resultsDir,'co2CapCE' + str(currYear) + '.csv'))
    pd.Series(currWSGenFracOfDemand).to_csv(os.path.join(resultsDir,'wsGenFracOfDemandCE' + str(currYear) + '.csv'))
    pd.Series(currWindGenFracOfDemand).to_csv(os.path.join(resultsDir,'windGenFracOfDemandCE' + str(currYear) + '.csv'))
    pd.Series(currOSWCapReq).to_csv(os.path.join(resultsDir,'oswCapReqCE' + str(currYear) + '.csv'))
    pd.Series(maxREInRegion).to_csv(os.path.join(resultsDir,'maxREInRegion' + str(currYear) + '.csv'))
    maxCapPerTech = maxCapPerTechOrig.copy()

    # ###############PREPARE INPUTS FOR CEM
    # Whether IRA tax credits are still in effect
    ira = (extract0dVarResultsFromGAMSModel(priorCEModel,'vCO2emsannual') > co2EmsInitial*.25) if priorCEModel != None else True #IRA ends when emissions reach 75% reducion relative to 2022
    pd.Series(ira).to_csv(os.path.join(resultsDir,'iraInEffectCE' + str(currYear) + '.csv'))
    if priorCEModel != None:
        print(extract0dVarResultsFromGAMSModel(priorCEModel,'vCO2emsannual'),co2EmsInitial*.25,extract0dVarResultsFromGAMSModel(priorCEModel,'vCO2emsannual') > co2EmsInitial*.25)

    # Update new technology and fuel price data    
    newTechsCE = getNewTechs(regElig, regCostFrac, currYear, stoInCE, seasStoInCE,
            fuelPrices, yearIncDACS, yearIncHydrogen,yearIncCCS, yearIncNuclear, transRegions, contFlexInelig, weatherYears, climateScenario, genFleet, ira, lbToShortTon)
    genFleet = updateFuelPricesAndCosts(genFleet, currYear, fuelPrices, regCostFrac)

    # Retire units and create fleet for current CE loop
    if priorCEModel != None:                    # if not in first CE loop
        genFleet = retireUnitsCFPriorCE(genFleet, genFleetPriorCE, retirementCFCutoff,
            priorCEModel, priorHoursCE, ptEligRetCF, currYear)
    genFleet, genFleetForCE = createFleetForCurrentCELoop(genFleet, currYear, retireByAge)
    genFleet.to_csv(os.path.join(resultsDir, 'genFleetPreCEPostRetirements' + str(currYear) + '.csv'))
    genFleetForCE.to_csv(os.path.join(resultsDir, 'genFleetForCEPreRECombine' + str(currYear) + '.csv'))
    print('MODIFIED LIFETIMES SO ONLY THERMAL PLANTS RETIRE')

    # Combine wind, solar, and storage plants by region
    genFleetForCE = combineWindSolarStoPlants(genFleetForCE)
    
    # Get renewable CFs by plant and region and calculate net demand by region
    print('Loading RE data')
    windGen, offWindGen, solarGen, windGenRegion, solarGenRegion = getREGen(genFleet, 
        tzAnalysis, weatherYears, currYear, pRegionShapes, nonCCReanalysis, climateChange, climateScenario, region)

    # Demand values are missing some periods w/in the 10 years (EPWs don't have leap days), so trim RE to demand index
    windGen = windGen.loc[demandYearUpgrade0.index]
    if offWindGen is not None: offWindGen = offWindGen.loc[demandYearUpgrade0.index] #use 0 as a value when no offshore wind units exist
    solarGen = solarGen.loc[demandYearUpgrade0.index]
    windGenRegion = windGenRegion.loc[demandYearUpgrade0.index]
    solarGenRegion = solarGenRegion.loc[demandYearUpgrade0.index]

    # Trim RE to the hours in full demand timeseries (missing 1 hour @ start and last week of 10-year period)
    netDemand = demandYearUpgrade0 - windGenRegion - solarGenRegion

    # Set planning reserve w/ forced upgrade.
    if forceBldgUpgrade != 0:
        demandWithUpgrade = demand.sum(axis=1) + demandSavings[forceBldgUpgrade].sum(axis=1)
    else:
        demandWithUpgrade = demand.sum(axis=1)
    demandWithUpgrade.to_csv(os.path.join(resultsDir, 'demandUsedToCalculatePRM' + str(currYear) + '.csv'))
    planningReserve = demandWithUpgrade.max()*(1+prm/100)
    planningReserveHour = demandWithUpgrade.idxmax()

    # Remove hydropower generation from demand using net-demand-based heuristic
    genFleetForCE,hydroGen,demandYearUpgrade0,planningReserve = processHydro(genFleetForCE, demandYearUpgrade0, netDemand, weatherYears, removeHydro, planningReserve) 
    genFleetForCE.to_csv(os.path.join(resultsDir, 'genFleetForCE' + str(currYear) + '.csv'))

    # Get hours included in CE model (representative + special blocks)
    hoursForCE,blockNamesChronoList,blockWeights = createTimeBlocks(demand,demandWeights)
    pd.Series(blockWeights).to_csv(os.path.join(resultsDir,'blockWeightsCE' + str(currYear) + '.csv'))

    # (hoursForCE, planningReserve, blockWeights, socScalars, planningReserveHour, blockNamesChronoList, 
    #     lastRepBlockNames, specialBlocksPrior) = getHoursForCE(demand, netDemand, windGenRegion, solarGenRegion,
    #     daysPerBlock, daysPerPeak, currYear, resultsDir, numBlocks, prm, prmBasis, climateChange)
    
    # Get CFs for new wind and solar sites and add wind & solar sites to newTechs
    newCfs,maxCapPerTech = getNewRenewableCFs(genFleet, tzAnalysis, weatherYears, currYear, 
        pRegionShapes, nonCCReanalysis, climateChange, climateScenario, region, maxCapPerTech, reDensity)
    newTechsCE,newCfs,maxCapPerTech = addWSSitesToNewTechs(newCfs, newTechsCE, pRegionShapes, reDownFactor, maxCapPerTech)
    pd.Series(maxCapPerTech).to_csv(os.path.join(resultsDir,'buildLimitsForCE' + str(currYear) + '.csv'))

    # Demand values are missing some periods w/in the 10 years (EPWs don't have leap days), so trim RE to demand index
    newCfs = newCfs.loc[demandYearUpgrade0.index]
    
    # Calculating thermal power plant & thermal line capacity deratings, FORs, and capacity eligibilities towards PRM
    print('Calculating deratings and capacity adjustments')
    capDerates,capDeratesTechs,newTechsCE = calculatePlantCapacityDerates(genFleetForCE, newTechsCE, demandYearUpgrade0, weatherYears, climateScenario, compressedGens, tzAnalysis)
    lineDerates = calculateLineThermalDerates(lineLimits, demandYearUpgrade0) 
    fors,windFOR,solarFOR,forsTechs,prmEligWindSolar = importPRMCapacityAdjustments(genFleetForCE, 
                            newTechsCE, demandYearUpgrade0, prmBasis, region, nonCCReanalysis, weatherYears, compressedGens, climateScenario, tzAnalysis)

    # Demand values are missing some periods w/in the 10 years (EPWs don't have leap days), so trim derates and fors to demand index
    capDerates = capDerates.loc[demandYearUpgrade0.index]
    capDeratesTechs = capDeratesTechs.loc[demandYearUpgrade0.index]
    fors = fors.loc[demandYearUpgrade0.index]
    forsTechs = forsTechs.loc[demandYearUpgrade0.index]

    # Set reserves for existing and incremental reserves for new generators
    print('Calculating reserves')
    if includeRes:
        cont, regUp, flex, regDemand, regUpSolar, regUpWind, flexSolar, flexWind = calcWWSISReserves(windGenRegion, solarGenRegion, demand, regLoadFrac,
                                                                                                     contLoadFrac, regErrorPercentile, flexErrorPercentile)
        regUpInc, flexInc = getIncResForAddedRE(newCfs, regErrorPercentile, flexErrorPercentile)
    else:
        cont, regUp, flex, regDemand, regUpSolar, regUpWind, flexSolar, flexWind, regUpInc, flexInc = createEmptyReserveDfs(windGenRegion, newCfs)

    # Get timeseries hours for CE (demand, wind, solar, new wind, new solar, reserves) & save dfs
    (demandCE, windGenCE, solarGenCE, newCfsCE, contCE, regUpCE, flexCE, regUpIncCE, 
        flexIncCE, forsCE, forsTechsCE, capDeratesCE, capDeratesTechsCE, lineDeratesCE) = isolateDataInCEHours(hoursForCE, 
        demandYearUpgrade0, windGenRegion, solarGenRegion, newCfs, cont, regUp, flex, regUpInc, flexInc, fors, forsTechs, capDerates, capDeratesTechs, lineDerates)
    
    # Get total hydropower generation potential by block for CE. If running CESM, get monthly hydropower generation to populate monthly generation constraints.
    [hydroGenCE] = isolateDataInCEBlocks(hoursForCE,hydroGen)
    if climateScenario is not None and 'rcp' not in climateScenario: #running CESM 
        hydroGenMonthlyCE = hydroGen.groupby(pd.Grouper(freq="M")).sum()    
        hydroGenMonthlyCE['dt'] = hydroGenMonthlyCE.index
        hydroGenMonthlyCE.index = ['month'+str(c)+'h' for c in range(hydroGenMonthlyCE.shape[0])]
        hydroGenMonthlyCE.to_csv(os.path.join(resultsDir,'hydroGenMonthlyCE'+str(currYear)+'.csv'))
    else: #not running CESM
        hydroGenMonthlyCE = 'NA'

    # Save CE inputs
    for df, n in zip([windGen, solarGen, windGenRegion, solarGenRegion, newCfs, demandYearUpgrade0, netDemand, cont, regUp, flex, regUpInc, flexInc, regDemand, regUpSolar, regUpWind, flexSolar, flexWind, hydroGen, fors, forsTechs, capDerates, capDeratesTechs, lineDerates],
                     ['windGen','solarGen','windGenRegion','solarGenRegion','windSolarNewCFs','demand','netDemand','contRes','regUpRes','flexRes','regUpInc','flexInc','regUpDemComp','regUpSolComp','regUpWinComp','flexSolComp','flexWinComp','hydroGen','fors','forTechs','capDerates','capDeratesTechs','lineDerates']):
        df.to_csv(os.path.join(resultsDir, n + 'FullYr' + str(currYear) + '.csv'))
    if offWindGen is not None: offWindGen.to_csv(os.path.join(resultsDir, 'windOffshoreGenFullYr' + str(currYear) + '.csv'))
    for df, n in zip([demandCE, windGenCE, solarGenCE, newCfsCE, newTechsCE, contCE, regUpCE, flexCE, regUpIncCE, flexIncCE, hydroGenCE, forsCE, forsTechsCE, capDeratesCE, capDeratesTechsCE, lineDeratesCE, hoursForCE],
                     ['demand', 'windGen', 'solarGen','windAndSolarNewCFs','newTechs','contRes','regUpRes','flexRes','regUpInc','flexInc','hydroGen','fors','forTechs','capDerates','capDeratesTechs','lineDerates','hoursByBlock']):
        df.to_csv(os.path.join(resultsDir, n + 'CE' + str(currYear) + '.csv'))
    for scalar, n in zip([windFOR,solarFOR,prmEligWindSolar,planningReserve,planningReserveHour],['windFOR','solarFOR','windSolarPRMElig','planningReserveCE','planningReserveHour']): 
        pd.Series(scalar).to_csv(os.path.join(resultsDir,n + str(currYear) + '.csv'))   
    # pd.DataFrame([[k, v] for k, v in socScalars.items()],columns=['block','scalar']).to_csv(os.path.join(resultsDir,'socScalarsCE' + str(currYear) + '.csv'))

    # ###############SET UP CAPACITY EXPANSION
    # Create GAMS workspace and database. Parameters and sets are put into database below.
    ws, db, gamsFileDir = createGAMSWorkspaceAndDatabase()

    # Write .gms files for inclusion in CE. These files vary based on CE parameters, e.g. treatment of time.
    # writeTimeDependentConstraints(blockNamesChronoList, stoInCE, seasStoInCE, gamsFileDir, ceOps, lastRepBlockNames, specialBlocksPrior, removeHydro, hydroGenMonthlyCE)
    writeTimeDependentConstraints(blockNamesChronoList, stoInCE, seasStoInCE, gamsFileDir, ceOps, removeHydro, hydroGenMonthlyCE)
    writeBuildVariable(ceOps, gamsFileDir)
    writeBuildingUpgradeFix(forceBldgUpgrade,gamsFileDir)

    # Enter sets and parameters into database
    genSet, hourSet, hourSymbols, zoneOrder, lineSet, zoneSet = edAndUCSharedFeatures(db, genFleetForCE, hoursForCE, demandCE, contCE,regUpCE,flexCE,
                                                                             demandShifter, demandShiftingBlock, rrToRegTime, rrToFlexTime, rrToContTime,
                                                                             solarGenCE, windGenCE, transRegions, lineLimits, transmissionEff, capDeratesCE,
                                                                             lineDeratesCE)  
    stoGenSet, stoGenSymbols = storageSetsParamsVariables(db, genFleetForCE, stoMkts, stoFTLabels)
    stoTechSet, stoTechSymbols = ceSharedFeatures(db, planningReserveHour, genFleetForCE, newTechsCE, planningReserve, discountRate, currCo2Cap,
                                      genSet, hourSet, hourSymbols, newCfsCE, maxCapPerTech, maxREPerZone, maxREInRegion, regUpIncCE, flexIncCE, 
                                      stoMkts, lineDists, lineCosts, lineSet, zoneOrder, ceOps, region, stoFTLabels, zoneSet, 
                                      forsCE, forsTechsCE, capDeratesTechsCE, windFOR, solarFOR, prmEligWindSolar, currWSGenFracOfDemand, 
                                      currWindGenFracOfDemand, demandSavings, demandCosts, interties, currOSWCapReq)
    if ceOps == 'UC': ucFeatures(db, genFleetForCE, genSet),
    ceTimeDependentConstraints(db, hoursForCE, blockWeights, ceOps,  
                genSet, genFleetForCE, stoGenSet,stoGenSymbols, newTechsCE, stoTechSet, stoTechSymbols, initSOCFraction,
                hydroGenCE, zoneSet, hydroGenMonthlyCE)

    # Run CE model
    print('Running CE for ' + str(currYear))
    capacExpModel, ms, ss = runGAMS('CEWith{o}.gms'.format(o=ceOps), ws, db)

    # ########## SAVE AND PROCESS CE RESULTS
    pd.Series([ms,ss],index=['ms','ss']).to_csv(os.path.join(resultsDir, 'msAndSsCE' + str(currYear) + '.csv'))
    saveCapacExpOperationalData(capacExpModel, genFleetForCE, newTechsCE, hoursForCE, transRegions, lineLimits, resultsDir, 'CE', currYear)
    newGens,newStoECap,newStoPCap,newLines = saveCEBuilds(capacExpModel, resultsDir, currYear)
    genFleet = addNewGensToFleet(genFleet, newGens, newStoECap, newStoPCap, newTechsCE, currYear)
    # genFleet = retireUnitsByCF(genFleet,hoursForCE,capacExpModel,currYear,planningReserve,planningReserveHour)
    lineLimits = addNewLineCapToLimits(lineLimits, newLines)
    genFleet.to_csv(os.path.join(resultsDir, 'genFleetAfterCE' + str(currYear) + '.csv'))
    lineLimits.to_csv(os.path.join(resultsDir, 'lineLimitsAfterCE' + str(currYear) + '.csv'))

    bldgUpgrade = extract1dVarResultsFromGAMSModel(capacExpModel,'vNBldgUpgrade')
    pd.Series(bldgUpgrade).to_csv(os.path.join(resultsDir, 'vNBldgUpgrade' + str(currYear) + '.csv'))
    importsDf = pd.DataFrame(columns=interties.index,index=hoursForCE.index)
    for rec in capacExpModel.out_db['vImport']: importsDf.loc[rec.key(1),rec.key(0)] = rec.level
    importsDf.to_csv(os.path.join(resultsDir, 'vImport' + str(currYear) + '.csv'))

    return (genFleet, genFleetForCE, lineLimits, capacExpModel, hoursForCE)

# ###############################################################################
# ###############################################################################
# ###############################################################################

# ###############################################################################
# ################## GAMS FUNCTIONS #############################################
# ###############################################################################
def createGAMSWorkspaceAndDatabase():
    gamsFileDir = 'GAMS'
    gamsSysDir = '/home/mtcraig/gams40_3'
    # gamsFileDir = 'C:\\Users\\mtcraig\\Desktop\\Research\\Models\\MacroCEMJill\\GAMS'
    # gamsSysDir = 'C:\\GAMS\\43'
    ws = GamsWorkspace(working_directory=gamsFileDir, system_directory=gamsSysDir)
    db = ws.add_database()
    return ws, db, gamsFileDir

def runGAMS(gamsFilename, ws, db):
    t0 = time.time()
    model = ws.add_job_from_file(gamsFilename)
    opts = GamsOptions(ws)
    opts.defines['gdxincname'] = db.name
    model.run(opts, databases=db)
    ms, ss = model.out_db['pModelstat'].find_record().value, model.out_db['pSolvestat'].find_record().value
    if (int(ms) != 8 and int(ms) != 1) or int(ss) != 1: print('***********************Modelstat & solvestat:', ms, ' & ', ss, ' (ms1 global opt, ms8 int soln, ss1 normal)')
    print('Time (mins) for GAMS run: ' + str(round((time.time()-t0)/60)))
    return model, ms, ss

def edAndUCSharedFeatures(db, genFleet, hours, demand, contRes, regUpRes, flexRes, demandShifter, demandShiftingBlock, rrToRegTime, rrToFlexTime,
                          rrToContTime, hourlySolarGen, hourlyWindGen, transRegions, lineLimits, transmissionEff, capDeratesCE, lineDeratesCE, cnse=10000, co2Price=0):
    # SETS
    genSet = addGeneratorSets(db, genFleet)
    hourSet, hourSymbols = addHourSet(db, hours)
    zoneSet,zoneSymbols,zoneOrder = addZoneSet(db, transRegions)
    lineSet,lineSymbols = addLineSet(db, lineLimits)

    # PARAMETERS
    # Demand and reserves
    addDemandParam(db, demand, hourSet, zoneSet, demandShifter, demandShiftingBlock, mwToGW)
    addReserveParameters(db, contRes, regUpRes, flexRes, rrToRegTime, rrToFlexTime, rrToContTime, hourSet, zoneSet, mwToGW)

    # CO2 cap or price
    addCo2Price(db, co2Price)

    # Generators
    addGenParams(db, genFleet, genSet, mwToGW, lbToShortTon, zoneOrder)
    addCapacityDerates(db, genSet, hourSet, capDeratesCE)
    addExistingRenewableMaxGenParams(db, hourSet, zoneSet, hourlySolarGen, hourlyWindGen, mwToGW)
    addSpinReserveEligibility(db, genFleet, genSet)
    addCostNonservedEnergy(db, cnse)

    # Transmission lines
    addLineParams(db,lineLimits, transmissionEff, lineSet, zoneOrder, mwToGW)
    addLineDerates(db, lineSet, hourSet, lineDeratesCE)
    return genSet, hourSet, hourSymbols, zoneOrder, lineSet, zoneSet

def storageSetsParamsVariables(db, genFleet, stoMkts, stoFTLabels):
    (stoGenSet, stoGenSymbols) = addStoGenSets(db, genFleet, stoFTLabels)
    addStorageParams(db, genFleet, stoGenSet, stoGenSymbols, mwToGW, stoMkts)
    return stoGenSet, stoGenSymbols

def ed(db, socInitial, stoGenSet):
    addStorageInitSOC(db, socInitial, stoGenSet, mwToGW)

def ucFeatures(db, genFleet, genSet):
    addGenUCParams(db, genFleet, genSet, mwToGW)
    
def uc(db, stoGenSet, genSet, socInitial, onOffInitial, genAboveMinInitial, mdtCarriedInitial):
    addStorageInitSOC(db, socInitial, stoGenSet, mwToGW)
    addEguInitialConditions(db, genSet, onOffInitial, genAboveMinInitial, mdtCarriedInitial, mwToGW)

def ceSharedFeatures(db, planningReserveHour, genFleet, newTechs, planningReserve, discountRate, co2Cap, 
        genSet, hourSet, hourSymbols, newCfs, maxCapPerTech, maxREPerZone, maxREInRegion, regUpInc, 
        flexInc, stoMkts, lineDists, lineCosts, lineSet, zoneOrder, ceOps, region, stoFTLabels, zoneSet,
        forsCE, forsTechsCE, capDeratesTechsCE, windFOR, solarFOR, prmEligWindSolar, currWSGenFracOfDemand, currWindGenFracOfDemand,
        demandSavings, demandCosts, interties, currOSWCapReq):
    addBuildingSetsAndParameters(db,demandSavings,demandCosts,zoneSet,hourSet,mwToGW)
    addIntertieSetsAndParameters(db,interties,zoneSet,mwToGW,zoneOrder)

    # Sets
    addPlanningReserveHourSubset(db, planningReserveHour)
    addStorageSubsets(db, genFleet, stoFTLabels)
    (techSet, renewTechSet, stoTechSet, stoTechSymbols, thermalSet, dacsSet, CCSSet) = addNewTechsSets(db, newTechs)

    # Long-term planning parameters
    addPlanningReserveParam(db, planningReserve, mwToGW)
    addDiscountRateParam(db, discountRate)
    addCO2Cap(db, co2Cap)

    # New tech parameters
    addGenParams(db, newTechs, techSet, mwToGW, lbToShortTon, zoneOrder, True)
    addCapacityDerates(db, techSet, hourSet, capDeratesTechsCE, True)
    addHourlyFORs(db, forsCE, genSet, hourSet)
    addHourlyFORs(db, forsTechsCE, techSet, hourSet, True)
    addFORScalars(db, windFOR, solarFOR, prmEligWindSolar)
    addTechCostParams(db, newTechs, techSet, stoTechSet, mwToGW)
    addRenewTechCFParams(db, renewTechSet, hourSet, newCfs)
    addMaxNewBuilds(db, newTechs, thermalSet, renewTechSet, stoTechSet, dacsSet, CCSSet, zoneSet, maxCapPerTech, maxREPerZone, maxREInRegion, mwToGW)
    addWSMinGen(db, currWSGenFracOfDemand, currWindGenFracOfDemand)
    addOSWCapReq(db, currOSWCapReq, genFleet)
    if ceOps == 'UC': addGenUCParams(db, newTechs, techSet, mwToGW, True)
    addResIncParams(db, regUpInc, flexInc, renewTechSet, hourSet)
    addSpinReserveEligibility(db, newTechs, techSet, True)
    addStorageParams(db, newTechs, stoTechSet, stoTechSymbols, mwToGW, stoMkts, True)
    addNewLineParams(db, lineDists, lineCosts, lineSet, maxCapPerTech, zoneOrder, region, mwToGW)
    return stoTechSet, stoTechSymbols

def ceTimeDependentConstraints(db, hoursForCE, blockWeights, ceOps, 
        genSet, genFleet, stoGenSet, stoGenSymbols, newTechs, stoTechSet, stoTechSymbols, 
        initSOCFraction, hydroGenCE, zoneSet, hydroGenMonthlyCE):
    addHourSubsets(db, hoursForCE, hydroGenMonthlyCE)
    addSeasonDemandWeights(db, blockWeights)
    # addBlockSOCScalars(db, socScalars)
    addStoInitSOCCE(db, genFleet, stoGenSet, stoGenSymbols, mwToGW, initSOCFraction)
    addStoInitSOCCE(db, newTechs, stoTechSet, stoTechSymbols, mwToGW, initSOCFraction, True)
    addHydroGenLimits(db, hydroGenCE, hydroGenMonthlyCE, zoneSet, mwToGW)

# ###############################################################################
# ###############################################################################
# ###############################################################################
