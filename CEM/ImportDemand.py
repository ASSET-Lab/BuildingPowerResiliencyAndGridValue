import operator, os, sys, pandas as pd

def importDemand(currYear,weatherYears,climateScenario,tzShift=-5):
    # dataDir='/nfs/turbo/seas-mtcraig/jmoraski/CEM_inputs/'
    dataDir=f'T:\\jmoraski\\CEM_inputs'

    #Representative demand time periods
    repPeriodsYear = (currYear - 40) if climateScenario == 'historical' else currYear
    demandRep = pd.read_csv(os.path.join(dataDir,'rep_days_csv',climateScenario + '_representative_periods_' + str(repPeriodsYear) + '.csv'),header=0,index_col=1,parse_dates=True)
    #Costs of upgrade scenarios
    cost = pd.read_csv(os.path.join(dataDir,'total_annualized_cost.csv'),header=0,index_col=0)

    #Add 'p' in front of region #s
    demandRep['reeds_balancing_area'] = demandRep['reeds_balancing_area'].replace({127:'p127',128:'p128'})

    #Convert rep period tz
    # demandRep.index = demandRep.index.shift(tzShift,freq='H') #shift for tz
  
    #Extract load and weights for rep periods
    upgrade0 = demandRep.loc[demandRep['upgrade']==0]
    baseDemand = upgrade0.pivot(columns='reeds_balancing_area',values='total_load')
    baseWeight = upgrade0.loc[upgrade0['reeds_balancing_area']=='p127','weight']

    #Calculate savings for rep periods 
    demandSavings = dict()
    for upgrade in demandRep['upgrade'].unique():
        if upgrade != 0: 
            upgradeDemand = demandRep.loc[demandRep['upgrade']==upgrade]
            upgradeDemand = upgradeDemand.pivot(columns='reeds_balancing_area',values='total_load')
            demandSavings[upgrade] = upgradeDemand - baseDemand

    #Full demand timeseries - 1 for each year
    demandFiles = list()
    for wy in weatherYears:
        demandFiles.append(pd.read_csv(os.path.join(dataDir,'full_timeseries_'+climateScenario,'full_timeseries_output_' + str(wy) + '.csv'),header=0,index_col=1,parse_dates=True))
    demand = pd.concat(demandFiles)

    #Add p infront of region #s
    demand['reeds_balancing_area'] = demand['reeds_balancing_area'].replace({127:'p127',128:'p128'})
    
    #Isolate upgrade 0 base demand for full year of demand
    demand = demand.loc[demand['upgrade']==0]
    demand = demand[['reeds_balancing_area','total_load']]
    # demand.index = pd.DatetimeIndex(demand['datetime_utc'])
    demand = demand.pivot(columns='reeds_balancing_area',values='total_load')

    #Shift tzs to local
    # demand.index = demand.index.shift(tzShift,freq='H') #shift for tz

    return baseDemand,baseWeight,demandSavings,cost,demand

def createTimeBlocks(demand,demandWeights):
    hoursForCE,blockNamesChronoList,blockWeights = pd.Series(9999,index=demand.index),list(),dict()
    blockCtr = 0

    for i in range(demand.index.shape[0]):
        if i>0: 
            shift = demand.index[i] - demand.index[i-1]
            if shift > pd.Timedelta(1,unit='h'):
                blockNamesChronoList.append(blockCtr)

                totalBlockWeight = demandWeights[demand.index[blockIdxStart:i]].sum()
                blockWeights[blockCtr] = totalBlockWeight

                blockCtr += 1
                blockIdxStart = i

            hoursForCE.iloc[i] = blockCtr

            if i == (demand.index.shape[0]-1):
                blockNamesChronoList.append(blockCtr)

                totalBlockWeight = demandWeights[demand.index[blockIdxStart:i]].sum()
                blockWeights[blockCtr] = totalBlockWeight
                                
        else:
            hoursForCE.iloc[i] = blockCtr
            blockIdxStart = 0

    return hoursForCE,blockNamesChronoList,blockWeights



    # if weatherYears == [2012]: demand = importEFSDemand(currYear,transRegions,demandScen) #2012
    # elif climateScenario != None and 'rcp' not in climateScenario[0]: 
    #     demand = importCESMDemand(weatherYears,climateScenario) #future years so use CC demand
    # elif climateScenario != None and 'rcp' in climateScenario[0]:

    # else: demand = importERA5Demand(weatherYears) #historic non 2012 data so use reanalysis


# def importEFSDemand(currYear,transRegions,demandScen):
#     print('Importing EFS demand')
#     #Initialize df
#     if currYear > 2050: currYear = 2050
#     dates = pd.date_range('1/1/'+str(currYear)+' 0:00','12/31/' + str(currYear) + ' 23:00',freq='H')
#     dates = dates[~((dates.month == 2) & (dates.day == 29))] #ditch leap day
#     demand = pd.DataFrame(index=dates)
    
#     #Read EFS data
#     filename = 'EP' + demandScen + '_FlexNONEload_hourly.csv'
#     rawDemand = pd.read_csv(os.path.join('Data','REEDS', filename), delimiter=',',header=0)
#     rawDemand = rawDemand.loc[rawDemand['year']==currYear]

#     #Iterate through dict of zone:p regions (from REEDS) & aggregate demand for p-regions
#     for zone,pRegions in transRegions.items():
#         for p in pRegions:
#             pDemand = rawDemand[p]
#             if zone in demand.columns:
#                 demand[zone] += pDemand.values
#             else:
#                 demand[zone] = pDemand.values
#     return demand

# def importCESMDemand(weatherYears,climateScenario):
#     print('Importing future climate demand')
#     demands = list()
#     #Import each ensemble member
#     for c in climateScenario:
#         filename = 'demand_' + c + '.csv'
#         demand = pd.read_csv(os.path.join('Data','CESM', filename), delimiter=',',header=0,index_col=0,parse_dates=True)

#         #Slim down to current years
#         demand = demand.loc[(demand.index.year >= weatherYears[0]) & (demand.index.year <= weatherYears[-1])]

#         #If multiple ensemble members, create MultiIndex with 2 levels (for ensemble member & region)        
#         if len(climateScenario)>1: 
#             idx = pd.MultiIndex.from_product([[c], demand.columns], names=['ensembleMember', 'region'])
#             demand.columns = idx

#         demands.append(demand)

#     #Combine into 1 array
#     demand = pd.concat(demands,axis=1)

#     return demand

# def importERA5Demand(weatherYears):
#     print('Importing historic demand')
#     # demand = pd.read_csv(os.path.join('Data','ERA5','hourlydemand_subregions_historic.csv'), delimiter=',',header=0,index_col=0,parse_dates=True)
#     demand = pd.read_csv('/nfs/turbo/seas-mtcraig/data_sharing/hari_paper1_ERA5/hourlydemand_subregions_historic2016to2021.csv', delimiter=',',header=0,index_col=0,parse_dates=True)
    
#     #Slim down to current years
#     demand = demand.loc[(demand.index.year >= weatherYears[0]) & (demand.index.year <= weatherYears[-1])]

#     return demand