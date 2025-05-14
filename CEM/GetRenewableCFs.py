import os, copy, datetime, pandas as pd, geopandas as gpd, datetime as dt, numpy as np, xarray as xr
from os import path
from statistics import mode
from Haochi_useful_script.util import *

def defineTransmissionRegions(region,balAuths='full'):
    transRegions = dict()
    if balAuths == 'full': #running full interconnection
        if region == 'ERCOT': transRegions = createERCOTGroupings(transRegions)
        elif region == 'EI': transRegions = createEIGroupings(transRegions)
        elif region == 'WECC': transRegions = createWECCGroupings(transRegions)
        elif region == 'NY': transRegions['p127'],transRegions['p128'] = ['p127'],['p128']
    else:
        sys.exit('defineTransmissionRegions not set up for non-full run!')
    return transRegions

#SEE BOTTOM OF SCRIPT FOR ALTERNATE STRUCTURE FOR WORKING WITH CLIMATE DATA

#Output: dfs of wind and solar generation (8760 dt rows, arbitrary cols)
def getREGen(genFleet, tgtTz, weatherYears, currYear, pRegionShapes, nonCCReanalysis, climateChange, climateScenario, region):
    #Importing single RE timeseries
    # if not climateChange or len(climateScenario)==1: 
    windGen, offWindGen, solarGen, windGenByRegion, solarGenByRegion = getSingleREGenTimeseries(genFleet,  
        tgtTz, weatherYears, currYear, pRegionShapes, nonCCReanalysis, climateChange, climateScenario, region) #, latlonRegion - output
    #Importing multiple RE timeseries for different ensemble members
    # else: 
    #     windGens,solarGens,windGenByRegions,solarGenByRegions = list(),list(),list(),list()
    #     for climateMember in climateScenario:
    #         windGen, solarGen, windGenByRegion, solarGenByRegion, latlonRegion = getSingleREGenTimeseries(genFleet, 
    #             tgtTz, weatherYears, currYear, pRegionShapes, nonCCReanalysis, climateChange, climateMember, region)

    #         windGen.columns = pd.MultiIndex.from_product([[climateMember], windGen.columns.to_list()], names=['ensembleMember', 'locs'])
    #         solarGen.columns = pd.MultiIndex.from_product([[climateMember], solarGen.columns.to_list()], names=['ensembleMember', 'locs'])
    #         windGenByRegion.columns = pd.MultiIndex.from_product([[climateMember], windGenByRegion.columns], names=['ensembleMember', 'region'])
    #         solarGenByRegion.columns = pd.MultiIndex.from_product([[climateMember], solarGenByRegion.columns], names=['ensembleMember', 'region'])

    #         windGens.append(windGen),solarGens.append(solarGen),windGenByRegions.append(windGenByRegion),solarGenByRegions.append(solarGenByRegion)

    #     #Combine into 1 array
    #     windGen,solarGen,windGenByRegion,solarGenByRegion = pd.concat(windGens,axis=1),pd.concat(solarGens,axis=1),pd.concat(windGenByRegions,axis=1),pd.concat(solarGenByRegions,axis=1)

    # print('********')
    # print(windGen)
    # print(offWindGen)
    # print(solarGen)
    # print(windGenByRegion)
    return windGen, offWindGen, solarGen, windGenByRegion, solarGenByRegion

def getSingleREGenTimeseries(genFleet, tgtTz, weatherYears, currYear, pRegionShapes, nonCCReanalysis, climateChange, climateMember, region):
    if currYear > 2050 and climateChange == False: currYear = 2050
    #Isolate wind & solar units
    windUnits, offWindUnits, solarUnits = genFleet.loc[genFleet['PlantType'] == 'Onshore Wind'], genFleet.loc[genFleet['PlantType'] == 'Offshore Wind'], genFleet.loc[genFleet['FuelType'] == 'Solar']
    #Get list of wind / solar sites in region and their CFs
    lats,lons,latsOn,lonsOn,latsOff,lonsOff,cf = loadData(weatherYears,pRegionShapes,climateMember,region,climateChange,nonCCReanalysis) #latlonRegion - output
    #Match to CFs
    getCFIndexCC(windUnits,latsOn,lonsOn,cf,'wind'),getCFIndexCC(solarUnits,latsOn,lonsOn,cf,'solar')
    if offWindUnits.shape[0]>0: getCFIndexCC(offWindUnits,latsOff,lonsOff,cf,'offwind')
    #Get hourly generation (8760 x n df, n = num generators). Use given met year data but set dt index to currYear.
    yrForCFs = weatherYears# if (climateChange or (nonCCReanalysis and region == 'WECC')) else [currYear] #if not CC, relabel fixed met year to future year; if CC, have many years
    windGen,solarGen = convertCFToGeneration(windUnits,cf['wind'],yrForCFs),convertCFToGeneration(solarUnits,cf['solar'],yrForCFs)
    if offWindUnits.shape[0]>0: offWindGen = convertCFToGeneration(offWindUnits,cf['offwind'],yrForCFs)
    else: offWindGen = None
    #Shift into right tz
    windGen, solarGen = shiftTz(windGen, tgtTz, weatherYears, 'wind'),shiftTz(solarGen, tgtTz, weatherYears, 'solar')
    if offWindUnits.shape[0]>0: offWindGen = shiftTz(offWindGen, tgtTz, weatherYears, 'wind')
    # if not (climateChange or (region=='WECC' and nonCCReanalysis)): 
    #     windGen, solarGen = shiftTz(windGen, tgtTz, currYear, 'wind'),shiftTz(solarGen, tgtTz, currYear, 'solar')
    # else:
    #     print('Not shifting RE timezone - either climate data or Zihan runs, both of which are already time shifted to local.')
    #Combine by region and fill in missing regions with zeros
    windGenByRegion = windGen.groupby(level='region', axis=1).sum() 
    if offWindUnits.shape[0]>0: windGenByRegion += offWindGen.groupby(level='region', axis=1).sum()
    solarGenByRegion = solarGen.groupby(level='region', axis=1).sum()
    for genByRegion in [windGenByRegion, solarGenByRegion]:
        regionsNoGen = [r for r in genFleet['region'].unique() if r not in genByRegion.columns]
        for r in regionsNoGen: genByRegion[r] = 0
    return windGen, offWindGen, solarGen, windGenByRegion, solarGenByRegion#, latlonRegion

#Get list of wind / solar sites in region and their CFs
def loadData(weatherYears,pRegionShapes,climateMember,region,climateChange,nonCCReanalysis):
    lats,lons,latsOn,lonsOn,latsOff,lonsOff, cf = loadTGWClimateData(weatherYears, pRegionShapes, climateMember) #latlonRegion - output
    # if climateChange: 
    #     if 'rcp' not in climateMember: lats, lons, cf, latlonRegion = loadCESMClimateData(weatherYears, pRegionShapes, climateMember)
    #     else: 
    # else: #if not climate change, only running single weather year, so take [0] index
    #     if not nonCCReanalysis: lats, lons, cf, latlonRegion = loadNRELData(weatherYears, pRegionShapes)
    #     else: lats, lons, cf, latlonRegion = loadERA5Data(weatherYears, pRegionShapes,region) 
    return lats,lons,latsOn,lonsOn,latsOff,lonsOff, cf#, latlonRegion

# Outputs: numpy arrays of lats and lons, and then a dictionary w/ wind and solar cfs
# as an np array of axbxc, where a/b/c = # lats/# longs/# hours in year


############ Haochi added function for TGW spatial projection handling #######
############ Haochi added function for TGW spatial projection handling #######
#Create function that finds locational index of nearest coords to given coordinates
def findNearestIndexes(coord,df):
    dfDiff = df.copy()
    lat,lon = coord[0],coord[1]
    dfDiff['diff'] = np.sqrt((df['xlats']-lat)**2+(df['xlongs']-lon)**2)
    idx = dfDiff['diff'].idxmin()
    return dfDiff.loc[idx,'south_north'],dfDiff.loc[idx,'west_east'],(dfDiff.loc[idx,'xlats'],dfDiff.loc[idx,'xlongs'])

# convert the index (x,y), to sliced index (x',y')
def conv_index_to_sliced_index(min_sn, max_sn, min_we, max_we, input_orig_sn, input_orig_we):
    # check input range
    if input_orig_sn < min_sn or input_orig_sn > max_sn:
        print(f'input_orig_sn{input_orig_sn}\n min_sn{min_sn} max_sn{max_sn}')
        raise ValueError("Original south_north index out of the sliced range.")
    if input_orig_we < min_we or input_orig_we > max_we:
        print(f'input_orig_we{input_orig_we}\n min_we{min_we} max_we{max_we}')
        raise ValueError("Original west_east index out of the sliced range.")

    # convert index
    output_sliced_sn = input_orig_sn - min_sn
    output_sliced_we = input_orig_we - min_we

    return output_sliced_sn, output_sliced_we

#from dask.diagnostics import ProgressBar
def get_spatial_conversion(region='NY',year='2003',File_prefix ='Solar_CF', Var_name='Solar_CF' ):
    '''
    This function is for conduct spatial projection conversion from TGW to lat/lon, for single years
    '''
    region = region
    transRegions = defineTransmissionRegions(region)
    pRegionShapes = loadRegions(transRegions)

    #read the spatial projection file
    tgwCoords = pd.read_csv(os.path.join('Data','TGW','tgwIdxsAndCoords.csv'),index_col=0,header=0)

    # get the boundary lat/lon, and add 0.2deg margin
    x1, y1, x2, y2 = pRegionShapes.total_bounds
    x1 = x1 -0.2
    y1 = y1 -0.2
    x2 = x2+ 0.2
    y2 = y2+ 0.2
    
    # read the temporal xarray file to get the lat/lon
    # ds= xr.open_dataset(f'S:\\TGW\\TGW_based_CFs\\TGW_Output_by_Year\\{File_prefix}_{year}.nc') 
    ds= xr.open_dataset(f'/nfs/turbo/seas-mtcraig-climate/TGW/TGW_based_CFs/TGW_Output_by_Year/{File_prefix}_{year}.nc')
    ds2 = ds.isel(Time=0)
    
    # find the range for isel
    masked_data = ds2.where((ds2.lat >= y1) & (ds2.lat <= y2) & (ds2.lon >= x1) & (ds2.lon <= x2), drop=False)
    masked_data[Var_name].values
    # check the null value
    non_null_positions = masked_data[Var_name].notnull()

    # find the min/max range in south_north index for none-null variable
    sn_non_null_indices = np.where(non_null_positions.any(dim='west_east'))[0]
    min_sn = sn_non_null_indices.min()
    max_sn = sn_non_null_indices.max()

    # find the min/max range in west_east index for none-null variable
    we_non_null_indices = np.where(non_null_positions.any(dim='south_north'))[0]
    min_we = we_non_null_indices.min()
    max_we = we_non_null_indices.max()

    # get the sliced .Xarray file
    ds = xr.open_dataset(f'S:\\TGW\\TGW_based_CFs\\TGW_Output_by_Year\\{File_prefix}_{year}.nc',chunks={}) #f'/nfs/turbo/seas-mtcraig-climate/TGW/TGW_based_CFs/TGW_Output_by_Year/{File_prefix}_{year}.nc'
    ds = ds[Var_name].reset_coords(['lat', 'lon'], drop=True)
    ds=ds.isel(south_north=slice(min_sn, max_sn + 1), west_east=slice(min_we, max_we + 1))
    
    # 3.create the nearest grid cell method
    new_lats = np.linspace(y1+0.1, y2-0.1, num=ds.shape[1])  # set the sample index is the original SN/WE resolution
    new_lons = np.linspace(x1+0.1, x2-0.1, num=ds.shape[2])  # set the sample index is the original SN/WE resolution

    new_long_grid, new_lat_grid = np.meshgrid(new_lons, new_lats, indexing='xy')  # establish the grid

    # initialize data and matrix
    interpolated_data = np.empty((len(ds.Time), len(new_lats), len(new_lons)))
    nearest_sn_indices = np.empty_like(new_lat_grid, dtype=int)
    nearest_we_indices = np.empty_like(new_long_grid, dtype=int)

    for i in range(new_lat_grid.shape[0]):
        for j in range(new_lat_grid.shape[1]):
            lat, lon = new_lat_grid[i, j], new_long_grid[i, j]
            sn, we, nearestCoord = findNearestIndexes((lat, lon), tgwCoords) 
            sn, we = conv_index_to_sliced_index(min_sn, max_sn, min_we, max_we, sn, we)
            nearest_sn_indices[i, j], nearest_we_indices[i, j] = sn, we
            
    # use xarray isel to chose data for each locations
    interpolated_da = ds.isel(south_north=xr.DataArray(nearest_sn_indices, dims=['latitude', 'longitude']),
                                       west_east=xr.DataArray(nearest_we_indices, dims=['latitude', 'longitude']))
    # re-define the coordinator
    interpolated_da = interpolated_da.assign_coords(latitude=new_lats, longitude=new_lons)
    interpolated_da = interpolated_da.rename({
        'latitude': 'lat',
        'longitude': 'lon',
        'Time': 'time'
    })

    return interpolated_da    

def get_spatial_conversion_offwind(region='NY',year='2003',File_prefix ='Solar_CF', Var_name='Solar_CF' ):
    '''
    This function is for conduct spatial projection conversion from TGW to lat/lon, for single years
    '''
    
    pRegionShapes = gpd.read_file(os.path.join('/nfs/turbo/seas-mtcraig/mtcraig/MacroCEMForHaochi/Python/Data','REEDS','Shapefiles','offshore_wind.shp'))

    #read the spatial projection file
    tgwCoords = pd.read_csv('tgwIdxsAndCoords.csv',index_col=0,header=0)

    # get the boundary lat/lon, and add 0.2deg margin
    x1, y1, x2, y2 = pRegionShapes.total_bounds
    x1 = x1 -0.2
    y1 = y1 -0.2
    x2 = x2+ 0.2
    y2 = y2+ 0.2
    
    # read the temporal xarray file to get the lat/lon
    ds= xr.open_dataset(f'/nfs/turbo/seas-mtcraig-climate/TGW/TGW_based_CFs/TGW_Output_by_Year/{File_prefix}_{year}.nc')
    ds2 = ds.isel(Time=0)
    
    # find the range for isel
    masked_data = ds2.where((ds2.lat >= y1) & (ds2.lat <= y2) & (ds2.lon >= x1) & (ds2.lon <= x2), drop=False)
    masked_data[Var_name].values
    # check the null value
    non_null_positions = masked_data[Var_name].notnull()

    # find the min/max range in south_north index for none-null variable
    sn_non_null_indices = np.where(non_null_positions.any(dim='west_east'))[0]
    min_sn = sn_non_null_indices.min()
    max_sn = sn_non_null_indices.max()

    # find the min/max range in west_east index for none-null variable
    we_non_null_indices = np.where(non_null_positions.any(dim='south_north'))[0]
    min_we = we_non_null_indices.min()
    max_we = we_non_null_indices.max()

    # get the sliced .Xarray file
    ds = xr.open_dataset(f'/nfs/turbo/seas-mtcraig-climate/TGW/TGW_based_CFs/TGW_Output_by_Year/{File_prefix}_{year}.nc',chunks={})
    ds = ds[Var_name].reset_coords(['lat', 'lon'], drop=True)
    ds=ds.isel(south_north=slice(min_sn, max_sn + 1), west_east=slice(min_we, max_we + 1))
    
    # 3.create the nearest grid cell method
    new_lats = np.linspace(y1+0.1, y2-0.1, num=ds.shape[1])  # set the sample index is the original SN/WE resolution
    new_lons = np.linspace(x1+0.1, x2-0.1, num=ds.shape[2])  # set the sample index is the original SN/WE resolution

    new_long_grid, new_lat_grid = np.meshgrid(new_lons, new_lats, indexing='xy')  # establish the grid

    # initialize data and matrix
    interpolated_data = np.empty((len(ds.Time), len(new_lats), len(new_lons)))
    nearest_sn_indices = np.empty_like(new_lat_grid, dtype=int)
    nearest_we_indices = np.empty_like(new_long_grid, dtype=int)

    for i in range(new_lat_grid.shape[0]):
        for j in range(new_lat_grid.shape[1]):
            lat, lon = new_lat_grid[i, j], new_long_grid[i, j]
            sn, we, nearestCoord = findNearestIndexes((lat, lon), tgwCoords) 
            sn, we = conv_index_to_sliced_index(min_sn, max_sn, min_we, max_we, sn, we)
            nearest_sn_indices[i, j], nearest_we_indices[i, j] = sn, we
            
    # use xarray isel to chose data for each locations
    interpolated_da = ds.isel(south_north=xr.DataArray(nearest_sn_indices, dims=['latitude', 'longitude']),
                                       west_east=xr.DataArray(nearest_we_indices, dims=['latitude', 'longitude']))
    # re-define the coordinator
    interpolated_da = interpolated_da.assign_coords(latitude=new_lats, longitude=new_lons)
    interpolated_da = interpolated_da.rename({
        'latitude': 'lat',
        'longitude': 'lon',
        'Time': 'time'
    })

    return interpolated_da    

############ Haochi added function for TGW spatial projection saving #######
############ Haochi added function for TGW spatial projection saving #######
def get_cf(climateScenario,region='NY',year='2003',File_prefix ='Solar_CF',Var_name='Solar_CF'):
    '''
    This function is for read single years CFs, 
    To check if the .NC file exists, if not, call 'get_spatial_conversion' function to finish the spatial projection conversion from TGW
    '''

    # base_dir = f'/nfs/turbo/seas-mtcraig/mtcraig/MacroCEMForHaochi/Python/Data/TGW'
    base_dir = f'T:\\mtcraig\\MacroCEMForHaochi\\Python\\Data\\TGW' 
    
    if climateScenario == 'historical': 
        climateLabel = climateScenario + '_1980-2020'
    else:
        climateLabel = climateScenario + '_2020-2060'
    base_dir = os.path.join(base_dir,climateLabel)

    file_name = f'{region}_{File_prefix}_{year}.nc'
    file_path = os.path.join(base_dir, file_name)  # Correct path joining
    
    if os.path.exists(file_path):  # Correct function for checking if the file exists
        print(f"existing {file_name}, read")
        ds_var = xr.open_dataset(file_path).load()  # Use the correct file path
        ds_var = ds_var[Var_name]
    else:
        print(f"NOT existing {file_name}, create")
        ds_var = get_spatial_conversion(region, year, File_prefix,Var_name)  # Make sure to assign the result of get_file
        
        # save files
        print(f'Saving Xarray {file_name}')
        # with ProgressBar():
        ds_var.to_netcdf(file_path)
    
    return ds_var

############ Haochi modified function for TGW  #######
############ Haochi modified function for TGW  #######
def loadTGWClimateData(weatherYears, pRegionShapes, tgwMember, dataDir=os.path.join('Data','TGW')):
    print('Using TGW met data!')
    
    weather_years = weatherYears
    wind_pow_gen,solar_pow_gen,offwind_pow_gen = list(),list(),list()
    
    # Concatenate multi years! 
    for year in weather_years:
        print(f'Handling {year}')
        wind_pow_gen.append(get_cf(tgwMember,year=year,File_prefix ='Wind_CF',Var_name='Wind_CF'))
        solar_pow_gen.append(get_cf(tgwMember,year=year,File_prefix ='Solar_CF',Var_name='Solar_CF'))
        offwind_pow_gen.append(get_cf(tgwMember,year=year,File_prefix ='offwind',Var_name='Wind_CF'))

    # Concatenate data along the time dimension
    windPowGen = xr.concat(wind_pow_gen, dim='time')
    solarPowGen = xr.concat(solar_pow_gen, dim='time')
    offwindPowGen = xr.concat(offwind_pow_gen, dim='time')
    
    # #NOTE: from here to below comment, code does not seem to do anything. 
    # #Get lat and lons for both datasets
    # latsPd,lonsPd = pd.DataFrame(solarPowGen['lat'], columns = ['lat']),pd.DataFrame(solarPowGen['lon'], columns=['lon'])
    
    # #Create df w/ lats & lons
    # latlonList = [(i, j) for i in latsPd.lat for j in lonsPd.lon]
    # latlonPd = pd.DataFrame(data=latlonList, columns=['lat', 'lon'])
    # latlonGpd = gpd.GeoDataFrame(latlonPd, geometry=gpd.points_from_xy(latlonPd.lon, latlonPd.lat))

    # #Join w/ p region shapes
    # latlonPshapeJoin = gpd.sjoin(latlonGpd, pRegionShapes, how="inner", op='intersects')
    # latlonPshapeJoin = latlonPshapeJoin.sort_values(by=['lat', 'lon'])
    # latlonPshapeJoin = latlonPshapeJoin.reset_index()

    # #Crosstab indicates which cells are in p-regions and which are not
    # latlonRegion = pd.crosstab(latlonPshapeJoin['lat'],latlonPshapeJoin['lon'])
    # #Correct longitudes back into 0 to 360 format
    # latlonRegion.columns = latlonRegion.columns + 360
    # #NOTE: from here to above comment, code does not seem to do anything. 

    #Store data
    cf = dict()
    cf["solar"] = np.array(solarPowGen[:])
    cf["wind"] = np.array(windPowGen[:])
    cf['offwind'] = np.array(offwindPowGen[:])
    cf['offwind'][cf['offwind']>1]=1 #truncate CFs to 1; there are some higher than 1 values    
    solarPowGen.close(), windPowGen.close(), offwindPowGen.close()

    #Reshape arrays so index is lat/lon/time instead of time/lat/lon
    if cf['solar'].shape[0] > 360: #0 idx is time dimension
        cf['solar'],cf['wind'],cf['offwind'] = np.swapaxes(cf['solar'],0,1),np.swapaxes(cf['wind'],0,1),np.swapaxes(cf['offwind'],0,1)
        cf['solar'],cf['wind'],cf['offwind'] = np.swapaxes(cf['solar'],1,2),np.swapaxes(cf['wind'],1,2),np.swapaxes(cf['offwind'],1,2)

    #Concat lats/lons for onshore & offshore locations - lats/lons for onshore& offshore are non-overlapping!
    latsOn,lonsOn = np.array(solarPowGen['lat'][:]),np.array(solarPowGen['lon'][:])
    latsOff,lonsOff = np.array(offwindPowGen['lat'][:]),np.array(offwindPowGen['lon'][:])
    allLats,allLons = np.concatenate([latsOn,latsOff]),np.concatenate([lonsOn,lonsOff])

    return allLats,allLons,latsOn,lonsOn,latsOff,lonsOff,cf#,latlonRegion

# Convert the latitude and longitude of the vg into indices for capacity factor matrix,
#then save that index into df
# More detail: The simulated capacity factor maps are of limited resolution. This function
#               identifies the nearest simulated location for renewable energy generators
#               and replaces those generators' latitudes and longitudes with indices for 
#               for the nearest simulated location in the capacity factor maps
def get_cf_index(RE_generators, powGen_lats, powGen_lons):
    RE_generators.loc[:,"lat idx"] = find_nearest_impl(RE_generators["Latitude"].astype(float), powGen_lats).astype(int)
    RE_generators.loc[:,"lon idx"] = find_nearest_impl(RE_generators["Longitude"].astype(float), powGen_lons).astype(int)

# Find index of nearest coordinate. Implementation of get_RE_index
def find_nearest_impl(actual_coordinates, discrete_coordinates):
    indices = []
    for coord in actual_coordinates:
        indices.append((np.abs(coord-discrete_coordinates)).argmin())
    return np.array(indices)

def getCFIndexCC(units,lats,lons,cf,re):
    centerLat,centerLon = (lats[0]+lats[-1])/2,(lons[0]+lons[-1])/2
    latStep,lonStep = lats[1]-lats[0],lons[1]-lons[0]
    # print(units)
    # print(units[['Latitude','Longitude']])
    # print(lats)
    # print(lons)
    # print(centerLat,centerLon,latStep,lonStep)
    for idx,unit in units.iterrows():
        lat,lon = unit['Latitude'],unit['Longitude']
        closestLat,closestLon = np.abs(lat-lats).argmin(),np.abs(lon-lons).argmin()
        # print(closestLat,closestLon,lats[closestLat],lons[closestLon],lat,lon,np.shape(lats),np.shape(lons))
        cfAtClosestPoint = cf[re][closestLat][closestLon]
        while np.isnan(cfAtClosestPoint[0]):
            lat,lon = lat + ((centerLat-lat)/abs(centerLat-lat))*latStep,lon + ((centerLon-lon)/abs(centerLon-lon))*lonStep
            closestLat,closestLon = np.abs(lat-lats).argmin(),np.abs(lon-lons).argmin()
            cfAtClosestPoint = cf[re][closestLat][closestLon] 
        units.loc[idx,'lat idx'],units.loc[idx,'lon idx'] = int(closestLat),int(closestLon)
    units['lat idx'] = units['lat idx'].astype(int)
    units['lon idx'] = units['lon idx'].astype(int)

# Find expected hourly capacity for RE generators. Of shape (8760 hrs, num generators)
# When running CC, using average demand, so don't scale up generation by 24
def convertCFToGeneration(RE_generators,cf,yrForCFs):
    #Pull number of time steps from CF array (indexed by lat/lon/time, so time is idx 2)
    tSteps = cf.shape[2] 
    f = 'H' if tSteps >=8760 else 'D' #f = 'H' if tSteps == 8760 else 'D' #2/1/23
    # repeat nameplate capacity in array
    RE_nameplate = np.tile(RE_generators["Capacity (MW)"].astype(float),(tSteps,1))
    # multiply by variable hourly capacity factor
    times = np.tile(np.arange(tSteps),(RE_generators['Capacity (MW)'].size,1)).T 
    RE_capacity = np.multiply(RE_nameplate, cf[RE_generators["lat idx"], RE_generators["lon idx"], times])
    # convert to datetime index
    idx = pd.date_range('1/1/'+ str(yrForCFs[0]) + ' 0:00','12/31/' + str(yrForCFs[-1]) + ' 23:00',freq=f)
    # create df with 2 sets of column labels
    reGen = pd.DataFrame(RE_capacity,columns=[RE_generators['GAMS Symbol'].values,RE_generators['region'].values])
    # drop excess time in leap year from datetime index and from CFs
    # if (len(idx) > 365 and f == 'D') or (len(idx) > 8760 and f == 'H'): idx = idx.drop(idx[(idx.month==2) & (idx.day ==29)]) #2/1/23
    # reGen = reGen.iloc[:len(idx)]
    # add dt index & relabel columns
    reGen.index = idx
    reGen.columns.names = ['GAMS Symbol','region']
    return reGen

#shift tz from UTC to local time
def shiftTz(reGen,tz,weatherYears,reType,tzOffsetDict = {'PST':-7,'CST': -6,'EST': -5}):
    origIdx = reGen.index
    reGen.index = reGen.index.shift(tzOffsetDict[tz],freq='H') #shift for tz
    reGen = reGen[reGen.index.year.isin(weatherYears)] #filter to year of interest
    reGen = reGen.append([reGen.iloc[-1]]*abs(tzOffsetDict[tz]),ignore_index=True) #replace dropped hours at end of year so have full year
    if reType=='solar': reGen.iloc[-5:] = 0 #set nighttime hours to 0 for newly added hours in prior row
    reGen.index=origIdx
    return reGen


#ALTERNATE STRUCTURE
# set year
# isolate W&S gens
# load data
# for w or solar:
#     for each gen in w or s:
#         get closest cell
#         get CF
#         gen = cap * CF
#         add to list
#     concat w or s gens
# concat for re type
# groupby for regional gen
# fill in missing regions
# return