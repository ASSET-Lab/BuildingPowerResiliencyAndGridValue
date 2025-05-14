#Michael Craig
#October 12, 2016
#Remove hydro (normal + pumped storage) units from fleet and subtract their monthly average generation
#from demand profile.

import copy, operator, os, numpy as np, pandas as pd, calendar
from CombinePlants import combinePlantsByRegion

#Inputs: gen fleet, demand and net demand dataframes
#Outputs: demand & net demand minus hydropower generation
def processHydro(genFleet, demand, netDemand, weatherYears, removeHydro, planningReserve):
    months,regions = pd.date_range(start='1/1/'+str(weatherYears[0]),end='12/31/'+str(weatherYears[-1]),freq='M'),netDemand.columns

    # print('in processhydro!')
    # print(months)
    if 'Hydro' in genFleet['FuelType'].unique():
        #Get installed capacity by region
        capByRegion,plantCapsByRegion = getHydroCapacByRegion(genFleet,regions)

        #Load historic or future data & (if EIA923) take extra steps to get monthly generation by region
        gen,netGenCols = import923() 
        monthGenByRegion = getMonthlyGenByRegion(capByRegion,plantCapsByRegion,gen,months,netGenCols) 
        # print(monthGenByRegion)
        monthGenByRegion.index = monthGenByRegion.index.map(lambda x: x + pd.DateOffset(years=demand.index.year.unique()[0] - weatherYears[0]))
        # print(monthGenByRegion)

        #Convert EIA data from monthly to hourly or daily
        genPerStepAll = convertRegionalMonthlyToSubMonthlyGen(monthGenByRegion,capByRegion,netDemand)

        #If removing hydro gen from demand, remove hydro from fleet; otherwise aggregate hydro by zone
        if removeHydro:
            demand -= genPerStepAll
            planningReserve -= genFleet.loc[genFleet['FuelType'] == 'Hydro']['Capacity (MW)'].sum()
            genFleet = genFleet.loc[genFleet['FuelType'] != 'Hydro']
        else: 
            for r in genFleet['region'].unique(): genFleet = combinePlantsByRegion(genFleet,'FuelType','Hydro',r)
    else: #create dummy copy of hourGenAll for later functions
        genPerStepAll = demand.copy()
        genPerStepAll *= 0

    # print('exiting processhydro')
    return genFleet, genPerStepAll, demand, planningReserve

def getHydroCapacByRegion(genFleet,regions):
    #Initialize df
    capByRegion,plantCapsByRegion = pd.Series(0,index=regions),dict()
    #Pull out hydro units
    hydroUnits = genFleet.loc[genFleet['FuelType'] == 'Hydro']
    for region in regions:
        #Aggregate hydro generators to plants, then get regional total capacity
        hydroRegion = hydroUnits.loc[hydroUnits['region']==region]
        initCap = hydroRegion['Capacity (MW)'].sum()

        plantCaps = hydroRegion.groupby('ORIS Plant Code')['Capacity (MW)'].apply(lambda x: np.sum(x.astype(float))).reset_index()
        plantCaps.index = plantCaps['ORIS Plant Code']
        plantCapsByRegion[region] = plantCaps

        capByRegion[region] = plantCaps['Capacity (MW)'].sum()

        if initCap>1: assert((initCap-capByRegion[region])<.01*initCap)
    return capByRegion,plantCapsByRegion

#Get monthly generation and capacity by region
def getMonthlyGenByRegion(capByRegion,plantCapsByRegion,gen,months,netGenCols):
    #Initialize dfs
    monthGenByRegion = pd.DataFrame(0,index=months,columns=capByRegion.index)
    
    #Pull out hydro units
    for region in capByRegion.index:
        #Match generators from fleet to 923 data
        genRegion = plantCapsByRegion[region].merge(gen,how='left',left_index=True,right_index=True)

        #Get total month generation
        for monthIdx in months:
            monthGenByRegion.loc[monthIdx,region] = genRegion[netGenCols[monthIdx.month-1]].sum() 
    return monthGenByRegion

#Convert hydropower generation from monthly to submonthly values
def convertRegionalMonthlyToSubMonthlyGen(monthGenByRegion,capByRegion,netDemand):
    #Determine time step (hourly or daily)
    hrsPerStep = 1 #if run total daily demand, use following code to scale hourly to daily capacity: # hrsPerStep = 24 if (netDemand.index[1]-netDemand.index[0]).days==1 else 1. Note that CC data is average, so don't scale now.
    #Pull out hydro units
    genPerStepAll = pd.DataFrame(index=netDemand.index,columns=netDemand.columns)

    for region in netDemand.columns:
        totalCapac = capByRegion[region]
        for dtIdx in monthGenByRegion.index:
            #Get monthly generation value
            monthGen = monthGenByRegion.loc[dtIdx,region]

            #Get net demand for month
            monthNetDemand = netDemand.loc[(netDemand.index.month==dtIdx.month) & (netDemand.index.year==dtIdx.year),region]
                
            #Only run below code if monthGen>0
            if monthGen > 1e-3: 
                #Calculate normalized weights so sum = 1, avoiding special case of all hours having negative net demand
                monthNetDemand.loc[monthNetDemand<0] = 0 #sometimes negative, which messes w/ weights; so set these to zero                
                if int(monthNetDemand.max()) != 0: 
                    wt = monthNetDemand/(monthNetDemand.max())
                    wt = wt/wt.sum() 
                else: 
                    wt = pd.Series(1/monthNetDemand.shape[0],index=monthNetDemand.index,name=region)

                #Estimate hourly generation using weights
                genPerStep = monthGen * wt * hrsPerStep
                assert((genPerStep.sum()-monthGen)<.01*monthGen)

                #If generation exceeds capacity in hours, reallocate that generation surplus to other hours
                hoursAboveCap,hoursBelowCap = genPerStep.loc[genPerStep>=(totalCapac*hrsPerStep)],genPerStep.loc[genPerStep<(totalCapac*hrsPerStep)]
                surplus = (hoursAboveCap - totalCapac).sum()
                while surplus > 0 and hoursBelowCap.shape[0] > 0:
                    #Evenly spread surplus generation across hours below capacity
                    genPerStep[hoursBelowCap.index] += surplus/hoursBelowCap.shape[0]
                    #Cap generation at capacity
                    genPerStep[hoursAboveCap.index] = totalCapac
                    #Reassign hours and recalculate surplus
                    hoursAboveCap,hoursBelowCap = genPerStep.loc[genPerStep>=(totalCapac*hrsPerStep)],genPerStep.loc[genPerStep<(totalCapac*hrsPerStep)]
                    surplus = (hoursAboveCap - totalCapac).sum()

                #Might exit while loop w/ remaining surplus - if so, all hours are full, so truncate
                if surplus > 0: genPerStep.loc[genPerStep>=(totalCapac*hrsPerStep)] = totalCapac*hrsPerStep

                #Place gen for month & region into df w/ all months & regions
                genPerStepAll.loc[genPerStep.index,region] = genPerStep
            
                #Make sure all generation was allocated (or dropped as surplus)
                assert((monthGen - (genPerStep.sum()+surplus))<.0001*monthGen)
            else:
                genPerStepAll.loc[monthNetDemand.index,region] = 0
    return genPerStepAll

def import923(metYear=2019):
    print('Using 2019 hydropower data - could improve!')
    #Column labels vary slightly b/wn years; handle that here
    if metYear == 2012: netGenCols,rftCol = ['Netgen_'+calendar.month_name[i][:3] for i in range(1,13)],'Reported Fuel Type Code'
    elif metYear == 2019: netGenCols,rftCol = ['Netgen\n'+calendar.month_name[i] for i in range(1,13)],'Reported\nFuel Type Code'
    else: netGenCols,rftCol = ['Netgen\r\n'+calendar.month_name[i] for i in range(1,13)],'Reported\r\nFuel Type Code'
        
    #Import, skipping empty top rows
    yrGen = pd.read_csv(os.path.join('Data','EIA923','gen' + str(metYear) + '.csv'),skiprows=5,header=0,thousands=',')
    yrGen = yrGen[['Plant Id',rftCol]+netGenCols]

    #Data very rarely has a missing value - replace with a zero for simplicity
    yrGen = yrGen.replace('.',0)

    #Slim down to hydro facilities
    yrGen = yrGen.loc[yrGen[rftCol] == 'WAT']
    yrGen.drop([rftCol],axis=1,inplace=True)

    #Get rid of , text and convert to float (thousands=',' doesn't work above)
    for lbl in netGenCols:
        yrGen[lbl] = yrGen[lbl].astype(str).str.replace(',','')
        yrGen[lbl] = yrGen[lbl].astype(float)

    #Aggregate unit to plant level
    yrGen = yrGen.groupby('Plant Id').apply(lambda x: np.sum(x.astype(float))).reset_index(drop=True)

    #Reindex
    yrGen['Plant Id'] = yrGen['Plant Id'].astype(int)
    yrGen.index = yrGen['Plant Id']
    yrGen.drop(['Plant Id'],axis=1,inplace=True)

    return yrGen,netGenCols

#CESM hydropower generation data is a monthly value that sums HOURLY operations, so for a 30
#day month, it is the sum of 30*24 periods of generation. But all other values in CESM data
#(inc. demand) are daily average values. So need hydropower generation that corresponds
#to a monthly value that sumds DAILY operations. So divide hydropower generation by 24!
def importCESMHydro(weatherYears,climateScenario,hourlyToDailyGen=24,climDir=os.path.join('Data','CESM')):
    if len(climateScenario) > 1: sys.exit('Hydropower code not built for more than 1 CESM member!')

    #Load hydro data
    gen = pd.read_csv(os.path.join(climDir,'hydro_monthly_'+climateScenario[0]+'.csv'),header=0,index_col=0,parse_dates=True)

    #Scale hydropower to daily
    gen /= hourlyToDailyGen

    #Rename columns
    gen.rename(columns={'Desert Southwest':'Desert_Southwest','NWPP Central':'NWPP_Central','NWPP Northeast':'NWPP_NE','NWPP Northwest':'NWPP_NW'},inplace=True)

    #Select weather years
    gen = gen[gen.index.year.isin(weatherYears)]

    return gen


# def convert923ToSubMonthlyGen(yrGen,genFleet,netDemand,netGenCols):
#     #Determine time step (hourly or daily)
#     hrsPerStep = 1 #if run total daily demand, use following code to scale hourly to daily capacity: # hrsPerStep = 24 if (netDemand.index[1]-netDemand.index[0]).days==1 else 1. Note that CC data is average, so don't scale now.
#     #Pull out hydro units
#     hydroUnits = genFleet.loc[genFleet['FuelType'] == 'Hydro']
#     genPerStepAll = pd.DataFrame(index=netDemand.index,columns=netDemand.columns)
#     for region in hydroUnits['region'].unique():
#         #Aggregate hydro generators to plants, then get regional total capacity
#         hydroRegion = hydroUnits.loc[hydroUnits['region']==region]
#         initCap = hydroRegion['Capacity (MW)'].sum()
#         capac = hydroRegion.groupby('ORIS Plant Code')['Capacity (MW)'].apply(lambda x: np.sum(x.astype(float))).reset_index()
#         capac.index = capac['ORIS Plant Code']
#         totalCapac = capac['Capacity (MW)'].sum()
#         assert((initCap-totalCapac)<.01*initCap)
        
#         #Match EIA data to hydro units
#         genRegion = capac.merge(yrGen,how='left',left_index=True,right_index=True)

#         #Get hourly generation based on regional net demand
#         for mnth in range(1,13):
#             #Get total month generation
#             monthGen = genRegion[netGenCols[mnth-1]].sum()

#             #Calculate hourly weights to determine hourly generation
#             monthNetDemand = netDemand.loc[netDemand.index.month==mnth,region]
#             monthNetDemand.loc[monthNetDemand<0] = 0 #sometimes negative, which messes w/ weights; so set these to zero

#             #Calculate normalized weights so sum = 1, avoiding special case of all hours having negative net demand
#             if int(monthNetDemand.max()) != 0: 
#                 wt = monthNetDemand/(monthNetDemand.max())
#                 wt = wt/wt.sum() 
#             else: 
#                 wt = pd.Series(1/monthNetDemand.shape[0],index=monthNetDemand.index,name=region)

#             #Estimate hourly generation using weights
#             genPerStep = monthGen * wt * hrsPerStep
#             assert((genPerStep.sum()-monthGen)<.01*monthGen)

#             #If generation exceeds capacity in hours, reallocate that generation surplus to other hours
#             hoursAboveCap,hoursBelowCap = genPerStep.loc[genPerStep>=(totalCapac*hrsPerStep)],genPerStep.loc[genPerStep<(totalCapac*hrsPerStep)]
#             surplus = (hoursAboveCap - totalCapac).sum()
#             while surplus > 0 and hoursBelowCap.shape[0] > 0:
#                 #Evenly spread surplus generation across hours below capacity
#                 genPerStep[hoursBelowCap.index] += surplus/hoursBelowCap.shape[0]
#                 #Cap generation at capacity
#                 genPerStep[hoursAboveCap.index] = totalCapac
#                 #Reassign hours and recalculate surplus
#                 hoursAboveCap,hoursBelowCap = genPerStep.loc[genPerStep>=(totalCapac*hrsPerStep)],genPerStep.loc[genPerStep<(totalCapac*hrsPerStep)]
#                 surplus = (hoursAboveCap - totalCapac).sum()

#             #Might exit while loop w/ remaining surplus - if so, all hours are full, so truncate
#             if surplus > 0: genPerStep.loc[genPerStep>=(totalCapac*hrsPerStep)] = totalCapac*hrsPerStep

#             #Place gen for month & region into df w/ all months & regions
#             genPerStepAll.loc[genPerStep.index,region] = genPerStep
        
#             #Make sure all generation was allocated (or dropped as surplus)
#             assert((monthGen - (genPerStep.sum()+surplus))<.0001*monthGen)
#     return genPerStepAll
