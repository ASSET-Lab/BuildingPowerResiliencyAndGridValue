import os, xarray as xr
from Haochi_useful_script.util import *

def get_derate_file(climateScenario,region='NY',year='2003',File_prefix ='Derate'):
#This function is for read single years CFs, To check if the .NC file exists, 
#if not, call 'get_derate_spatial_conversion' function to finish the spatial projection conversion from TGW
#    from dask.diagnostics import ProgressBar
    # ProgressBar().register()

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
    else:
        print(f"NOT existing {file_name}, create")

        T2 = get_derate_spatial_conversion(year=year,File_prefix ='Derate',Var_name='T2')
        RH = get_derate_spatial_conversion(year=year,File_prefix ='Derate',Var_name='RH')
        PSFC = get_derate_spatial_conversion(year=year,File_prefix ='Derate',Var_name='PSFC')
        Q2 = get_derate_spatial_conversion(year=year,File_prefix ='Derate',Var_name='Q2')
        merged_dataset = xr.merge([T2, RH, PSFC, Q2])
        # keep varname as Q2
        merged_dataset = merged_dataset.rename({'T2': 'TREFHT', 'RH': 'RHREFHT', 'PSFC': 'PS'})
        merged_dataset['TREFHT'] = merged_dataset['TREFHT'] + 273.15
        ds_var = merged_dataset

        # save files
        print(f'Saving Xarray {file_name}')
        # with ProgressBar():
        ds_var.to_netcdf(file_path)
    
    return ds_var

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

#Import met vars for derates & TDFORs of thermal units (T, rh, pressure)
def importNonREMet(weatherYears,climateScenario,nonCCReanalysis):
    if climateScenario != None: 
        if 'rcp' not in climateScenario[0]: metVars = importCESMMetVars(weatherYears,climateScenario[0])
        else: metVars = importTGWMetVars(weatherYears,climateScenario[0])
    elif nonCCReanalysis: metVars = importERA5MetVars(weatherYears)
    else: metVars = None
    return metVars

def importERA5MetVars(weatherYears):
    #Load FORs file that contains temperatures (in C) and slice down to weather years
    # temps = xr.open_dataset(os.path.join('Data','ERA5','wecc_FOR_ERA5_hourly_PST2016to2021.nc'))
    temps = xr.open_dataset('/nfs/turbo/seas-mtcraig/data_sharing/hari_paper1_ERA5/wecc_FOR_ERA5_hourly_PST2016to2021.nc')
    temps = temps.sel(time=slice(str(weatherYears[0])+"-01-01", str(weatherYears[-1])+"-12-31"))
    return temps

def importCESMMetVars(weatherYears,cesmMember):
    print(cesmMember)
    temps = xr.open_dataset(os.path.join('Data','CESM','derate_fields_' + cesmMember + '.nc'))
    temps = temps.sel(time=slice(str(weatherYears[0])+"-01-01", str(weatherYears[-1])+"-12-31")) #T in K    
    return temps

def importTGWMetVars(weatherYears,tgwMember):
    # print(tgwMember)
    region = 'NY'
    weather_years = weatherYears
    derate_files = []
    print('derate climate scenario:',tgwMember)

    # Concatenate multi years! 
    for year in weather_years:
        print(f'Handling {year}_for_derate')
        derate_files.append(get_derate_file(tgwMember,region = region, year = year))

    # Concatenate data along the time dimension
    temps = xr.concat(derate_files, dim='time')
    
    return temps

### Haochi added function for spatial projection handling ## 
### Haochi added function for spatial projection handling ## 

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


def get_derate_spatial_conversion(region='NY',year='2003',File_prefix ='Derate', Var_name='T2' ):
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
    ds= xr.open_dataset(f'S:\\TGW\\TGW_based_CFs\\TGW_Output_by_Year\\{File_prefix}_{year}.nc') #f'/nfs/turbo/seas-mtcraig-climate/TGW/TGW_based_CFs/TGW_Output_by_Year/{File_prefix}_{year}.nc'
    # ds= xr.open_dataset(f'/nfs/turbo/seas-mtcraig-climate/TGW/TGW_based_CFs/TGW_Output_by_Year/{File_prefix}_{year}.nc')
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
    ds = xr.open_dataset(f'S:\\TGW\\TGW_based_CFs\\TGW_Output_by_Year\\{File_prefix}_{year}.nc',chunks={}) 
    # ds = xr.open_dataset(f'/nfs/turbo/seas-mtcraig-climate/TGW/TGW_based_CFs/TGW_Output_by_Year/{File_prefix}_{year}.nc',chunks={})
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



