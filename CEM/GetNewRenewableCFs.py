import os, copy, datetime, pandas as pd, datetime as dt, numpy as np, geopandas as gpd, shapely
from os import path
from GetRenewableCFs import *

#Output: dfs of wind and solar generation (8760 dt rows, arbitrary cols)
def getNewRenewableCFs(genFleet, tgtTz, weatherYears, currYear, pRegionShapes, 
            nonCCReanalysis, climateChange, climateScenario, region, maxCapPerTech, reDensity):
    print('Getting new renewable CFs!')
    #Importing single RE timeseries
    # print(climateScenario)
    newCfs,maxCapPerTech = getSingleNewRenewableCFsTimeseries(genFleet, tgtTz, weatherYears, currYear, 
            pRegionShapes, nonCCReanalysis, climateChange, climateScenario, region, maxCapPerTech, reDensity)
    
    # if not climateChange or len(climateScenario)==1: 
    #     newCfs,maxCapPerTech = getSingleNewRenewableCFsTimeseries(genFleet, tgtTz, weatherYears, currYear, 
    #         pRegionShapes, nonCCReanalysis, climateChange, climateScenario, region, maxCapPerTech, reDensity)
    # #Importing multiple RE timeseries for different ensemble members
    # else: 
    #     newCfsAll = list()
    #     for climateMember in climateScenario:
    #         newCfs,maxCapPerTech = getSingleNewRenewableCFsTimeseries(genFleet, tgtTz, weatherYears, currYear, 
    #             pRegionShapes, nonCCReanalysis, climateChange, climateMember, region, maxCapPerTech, reDensity)
    #         newCfs.columns = pd.MultiIndex.from_product([[climateMember], newCfs.columns], names=['ensembleMember', 'locs'])
    #         newCfsAll.append(newCfs)
    #     #Combine into 1 array
    #     newCfs = pd.concat(newCfsAll,axis=1)
    return newCfs,maxCapPerTech

#Output: dfs of wind and solar generation (8760 dt rows, arbitrary cols)
def getSingleNewRenewableCFsTimeseries(genFleet, tgtTz, weatherYears, currYear, 
        pRegionShapes, nonCCReanalysis, climateChange, climateMember, region, maxCapPerTech, reDensity):
    if currYear > 2050 and climateChange == False: currYear = 2050
    
    #Isolate wind & solar units
    windUnits,offWindUnits,solarUnits = genFleet.loc[genFleet['PlantType'] == 'Onshore Wind'], genFleet.loc[genFleet['PlantType'] == 'Offshore Wind'], genFleet.loc[genFleet['FuelType'] == 'Solar']
    
    #Get list of wind / solar sites in region.
    lats,lons,latsOn,lonsOn,latsOff,lonsOff,cf = loadData(weatherYears,pRegionShapes,climateMember,region,climateChange,nonCCReanalysis) #latlonRegion - output
    # print('******HEY')
    # print(cf)

    #Import available land per grid cell
    wArea,offwArea,sArea = importWSLandAvailable()

    #Match existing gens to CFs (use this for figuring our spare capacity available at each coordinate given existing generators)
    getCFIndexCC(windUnits,latsOn,lonsOn,cf,'wind'),getCFIndexCC(solarUnits,latsOn,lonsOn,cf,'solar')
    if offWindUnits.shape[0] > 0: getCFIndexCC(offWindUnits,latsOff,lonsOff,cf,'offwind')
    # if not climateChange: get_cf_index(windUnits, lats, lons),get_cf_index(solarUnits, lats, lons)
    # else: getCFIndexCC(windUnits,lats,lons,cf),getCFIndexCC(solarUnits,lats,lons,cf)

    #Calculate new CFs for given met year, but (in some cases) setting dt index to currYear
    yrForCFs = weatherYears# if (climateChange or (nonCCReanalysis and region == 'WECC')) else [currYear] #if not CC, relabel fixed met year to future year; if CC, have many years
    
    #DO I EVEN NEED THIS!!!!!!!!!! since shapefile from jill will constrain area
    #Filter out new CF locations to states being analyzed
    # stateBounds = latlonRegion.reset_index(drop=True)
    # stateBounds.columns = range(stateBounds.columns.size)
    # cf = enforceStateBounds(cf, stateBounds)

    #Calculate new CF sites    
    windCfs,maxNewWindCaps = calcNewCfs(windUnits, latsOn, lonsOn, cf, 'wind', 'Onshore Wind', yrForCFs, wArea, reDensity)
    solarCfs,maxNewSolarCaps = calcNewCfs(solarUnits, latsOn, lonsOn, cf, 'solar', 'Solar', yrForCFs, sArea, reDensity)
    offwindCfs,maxNewOffWindCaps = calcNewCfs(offWindUnits, latsOff,lonsOff, cf, 'offwind', 'Offshore Wind', yrForCFs, offwArea, reDensity)

    #Modify max capacity of wind and solar using cell-specific available capacity
    for re in ['Onshore Wind','Offshore Wind','Solar']: maxCapPerTech.pop(re)
    maxCaps = pd.concat([maxNewWindCaps,maxNewSolarCaps,maxNewOffWindCaps])
    for i in maxCaps.index: maxCapPerTech[i] = maxCaps[i]

    #Shift into right tz
    windCfs, offwindCfs, solarCfs = shiftTz(windCfs, tgtTz, weatherYears, 'wind'),shiftTz(offwindCfs, tgtTz, weatherYears, 'wind'),shiftTz(solarCfs, tgtTz, weatherYears, 'solar')
    # if not (climateChange or (region=='WECC' and nonCCReanalysis)): 
    #     windCfs, solarCfs = shiftTz(windCfs, tgtTz, currYear, 'wind'), shiftTz(solarCfs, tgtTz, currYear, 'solar')
    # else:
    #     print('Not shifting RE timezone - either climate data or Zihan runs, both of which are already time shifted to local.')

    return pd.concat([windCfs, solarCfs, offwindCfs], axis=1),maxCapPerTech

#Import available land area per grid cell based on Grace WU SL work (see gis-work channel; https://www.pnas.org/doi/10.1073/pnas.2204098120)
def importWSLandAvailable(shapeDir=os.path.join('Data','RenewableShapefilesNYC')):
    #Import shapefile w/ available land
    wShape,sShape = gpd.read_file(os.path.join(shapeDir,'onshore_wind.shp')),gpd.read_file(os.path.join(shapeDir,'solar.shp'))
    offwShape = gpd.read_file(os.path.join(shapeDir,'offshore_wind.shp'))

    #Filter based on suitability value
    wShape,sShape = wShape[wShape['value'] >= 70],sShape[sShape['value'] >= 70]

    return wShape,offwShape,sShape

# #Import available land area per grid cell based on Grace WU SL work (see gis-work channel; https://www.pnas.org/doi/10.1073/pnas.2204098120)
# def importWSLandAvailable(region,climateChange,wsSLScenario='sl1'):
#     #Import available area per grid cell based on RE dataset and spatial resolution (WECC and CONUS are the same in WECC except at edges of WECC)
#     # #Used the following two lines for Hari paper 2 runs around July 29, 2023
#     # wArea = pd.read_csv(os.path.join('Data','WindSolarSiting',wsSLScenario + '_wind_cesm2_areaassessment_reproject_hariRuns.csv'),index_col=0,header=0)
#     # sArea = pd.read_csv(os.path.join('Data','WindSolarSiting',wsSLScenario + '_solar_cesm2_areaassessment_reproject_hariRuns.csv'),index_col=0,header=0)
#     if climateChange:
#         if region == 'WECC': 
#             wArea = pd.read_csv(os.path.join('Data','WindSolarSiting',wsSLScenario + '_wind_cesm2_wecc.csv'),index_col=0,header=0)
#             sArea = pd.read_csv(os.path.join('Data','WindSolarSiting',wsSLScenario + '_solar_cesm2_wecc.csv'),index_col=0,header=0)
#         else:
#             wArea = pd.read_csv(os.path.join('Data','WindSolarSiting',wsSLScenario + '_wind_cesm2_conus.csv'),index_col=0,header=0)
#             sArea = pd.read_csv(os.path.join('Data','WindSolarSiting',wsSLScenario + '_solar_cesm2_conus.csv'),index_col=0,header=0)
#     else:
#         if region == 'WECC': 
#             wArea = pd.read_csv(os.path.join('Data','WindSolarSiting',wsSLScenario + '_wind_era5_wecc.csv'),index_col=0,header=0)
#             sArea = pd.read_csv(os.path.join('Data','WindSolarSiting',wsSLScenario + '_solar_era5_wecc.csv'),index_col=0,header=0)
#         else:
#             wArea = pd.read_csv(os.path.join('Data','WindSolarSiting',wsSLScenario + '_wind_era5_conus.csv'),index_col=0,header=0)
#             sArea = pd.read_csv(os.path.join('Data','WindSolarSiting',wsSLScenario + '_solar_era5_conus.csv'),index_col=0,header=0)

#     #Fill NAs with 0
#     wArea['ia_area'].fillna(0,inplace=True),sArea['ia_area'].fillna(0,inplace=True)

#     #Isolate lat & lon
#     locs = wArea['centroid_text'].str.split(' ')
#     lons,lats = [float(v[0].split('(')[1]) for v in locs.values],[float(v[1].split(')')[0]) for v in locs.values]
#     lats,lons = pd.Series(lats,index=locs.index),pd.Series(lons,index=locs.index)
#     wArea['Latitude'],wArea['Longitude'] = lats,lons+360

#     locs = sArea['centroid_text'].str.split(' ')
#     lons,lats = [float(v[0].split('(')[1]) for v in locs.values],[float(v[1].split(')')[0]) for v in locs.values]
#     lats,lons = pd.Series(lats,index=locs.index),pd.Series(lons,index=locs.index)
#     sArea['Latitude'],sArea['Longitude'] = lats,lons+360

#     return wArea,sArea

def getSquareAroundPoint(lon,lat, delta_size=0.1): #same delta size as when initially create grid - SEE GETRENEWABLECFS
    
    point_coords = np.array([lon,lat])
    
    c1 = point_coords + [-delta_size,-delta_size]
    c2 = point_coords + [-delta_size,+delta_size]
    c3 = point_coords + [+delta_size,+delta_size]
    c4 = point_coords + [+delta_size,-delta_size]
    
    square_geom = shapely.geometry.Polygon([c1,c2,c3,c4])
    
    return square_geom


def calcNewCfs(existingGens, lats, lons, cf, re, reDensLabel, yrForCFs, availableArea, reDensity, maxAreaDir=os.path.join('Data','RenewableShapefilesNYC'), ftToM=.0929): 
    print('using projection for NY state only! in GetNewRenewableCFs')
    print('using .1 delta size - if changed in GetRenewableCFs, change in GetNewRenewableCFs')    
    print('Calculating new CFs for ' + re)

    #Pull number of time steps from CF array (indexed by lat/lon/time, so time is idx 2)
    tSteps = cf[re].shape[2] 
    f = 'H' if tSteps >=8760 else 'D' #f = 'H' if tSteps == 8760 else 'D' #2/1/23
    #For each lat/lon, check existing capacity and, if spare room for more renewables, add CFs
    cfs,maxNewCaps,maxedLocs,maxAreas = dict(),dict(),list(),dict()
    latDiff,lonDiff = abs(lats[1]-lats[0]),abs(lons[1]-lons[0])

    #Check if file already exists
    areasFile = False
    if os.path.exists(os.path.join(maxAreaDir,re+'MaxArea.csv')): 
        print('Loading max areas file')
        maxAreas = pd.read_csv(os.path.join(maxAreaDir,re+'MaxArea.csv'),index_col=[0,1]) #0,1 to form MultiIndex tuple
        areasFile = True

    for latIdx in range(len(lats)):
        for lonIdx in range(len(lons)):
            lat,lon = np.round(lats[latIdx],4), np.round(lons[lonIdx],4)
            assert(availableArea.crs == 'epsg:4326')

            #Get max available area for RE at location
            if not areasFile:
                cell = getSquareAroundPoint(lon,lat) #create square around point
                intersected = availableArea[availableArea.intersects(cell)] #intersect square with shapefile of available area
                intersected['geometry'] = intersected.geometry.to_crs("epsg:2263") #put into epsg that provides accurate area measures. https://epsg.io/2263 - this is for NEW YORK.             
                maxArea = intersected.geometry.area.sum() #per 2263 projection above, units are in ft! So area is square feet!
                maxAreas[(lat,lon)]=maxArea
            else:
                maxArea = maxAreas.loc[(lat,lon)].values[0] #extract value from MultiIndexed Series

            maxCapacity = maxArea*reDensity[reDensLabel]*ftToM/1e6 #MW [=ft^2 * W/m^2 * .0929 m^2/ft^2 * MW/1e6W] #capitalize

            #Get generators & total capacity at location
            if existingGens.shape[0] > 0: 
                gensAtLoc = existingGens.loc[(existingGens['lat idx'] == latIdx) & (existingGens['lon idx'] == lonIdx)]
                existingCap = gensAtLoc['Capacity (MW)'].astype(float).sum()
            else:
                existingCap = 0

            #Get CFs @ coordinates
            coordCfs = cf[re][latIdx, lonIdx, :]

            # print(np.isnan(coordCfs).all())

            #Filter out any coordinates with all NANs (these are sites not in EI)
            if not np.isnan(coordCfs).all():
                #Filter out coords w/ no gen
                if coordCfs.sum() > 0:
                    if existingCap > 0: print('Existing and max capacity at ' + str(lat) + str(lon) + ':' + str(existingCap) + ' ' + str(maxCapacity))
                    #If can fit more capacity in cell, save CFs and update max capacity
                    if existingCap < maxCapacity: 
                        cfs[re + 'lat' + str(lat) + 'lon' + str(lon)] = coordCfs
                        maxNewCaps[re + 'lat' + str(lat) + 'lon' + str(lon)] = maxCapacity - existingCap
                    #Zero out CFs if current capacity exceeds max possible capacity in that cell
                    else:
                        cfs[re + 'lat' + str(lat) + 'lon' + str(lon)] = 0
                        maxNewCaps[re + 'lat' + str(lat) + 'lon' + str(lon)] = 0
                        if maxCapacity > 1e-3: maxedLocs.append((lat,lon))
    print('No new ' + re + ' capacity at following sites due to existing RE capacity:', maxedLocs)

    # print(cfs)
    # print(maxNewCaps)
    
    #Create df of CFs and series of max capacities
    cfs,maxNewCaps = pd.DataFrame(cfs),pd.Series(maxNewCaps)
    #Create dt index
    idx = pd.date_range('1/1/' + str(yrForCFs[0]) + ' 0:00','12/31/' + str(yrForCFs[-1]) + ' 23:00', freq=f)
    #Drop leap year extra time steps from index and df
    # idx = idx.drop(idx[(idx.month == 2) & (idx.day == 29)]) #2/1/23
    # cfs = cfs.iloc[:len(idx)]
    #Add datetime index to cfs df
    cfs.index = idx
    #Save max areas
    if not areasFile: pd.Series(maxAreas,name='max areas (sq ft)').to_csv(os.path.join(maxAreaDir,re+'MaxArea.csv'))
    return cfs,maxNewCaps

def calcHaversine(lat1,lon1,lat2,lon2,R=6371): #R = Earth radius (km)
    lat1,lon1,lat2,lon2 = lat1*np.pi/180,lon1*np.pi/180,lat2*np.pi/180,lon2*np.pi/180
    a = np.sin((lat1-lat2)/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin((lon1-lon2)/2)**2
    c = 2 * np.arctan2(np.sqrt(a),np.sqrt(1-a))
    return R*c

def enforceStateBounds(cf, stateBounds):
   for re in cf:
        for row in stateBounds.index:
            for col in stateBounds.columns:
                cf[re][row,col] *= stateBounds.loc[row,col]
    # plotCFs(cf)
   return cf



# import matplotlib.pyplot as plt
# def plotCFs(cf):
#     avgCfs,lats,lons = np.zeros((23,23)),np.zeros(23),np.zeros(23)
#     for re in cf:
#         cfs = cf[re]
#         for lat in range(cfs.shape[0]):
#             for lon in range(cfs.shape[1]):
#                 avgCfs[lat,lon] = cfs[lat,lon].mean()
#                 # lats[lat] = cfs['lat'][lat]
#                 # lons[lon] = cfs['lon'][lon]

#         plt.figure()
#         ax = plt.subplot(111)
#         im = ax.contourf(avgCfs,cmap='plasma')#,extent = [np.min(lons),np.max(lons),np.min(lats),np.max(lats)])
#         cbar = ax.figure.colorbar(im, ax=ax)#, ticks=np.arange(vmin,vmax,int((vmax-vmin)/5)))
#         plt.title(re)
#     plt.show()
