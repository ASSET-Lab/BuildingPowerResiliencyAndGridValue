import pandas as pd, geopandas as gpd, numpy as np
from SetupTransmissionAndZones import assignGensToPRegions

def addWSSitesToNewTechs(newCfsOrig,newTechsCE,pRegionShapes,reDownFactor,maxCapPerTech):
    #Create copy for isolating sites
    newCfs = newCfsOrig.copy()

    sitesDfList = list()
    #For wind & solar, repeat tech row for each potential site, then remove original row
    for l,pt in zip(['wind','solar','offwind'],['Onshore Wind','Solar PV','Offshore Wind']): #ft ['Wind','Solar','Wind']
        #Get row of new techs CE w/ right fuel type
        # re = newTechsCE.loc[newTechsCE['FuelType']==ft]
        re = newTechsCE.loc[newTechsCE['PlantType']==pt]
        #Find new wind & solar sites w/ that fuel type
        sites = [c for c in newCfs if l in c]
        #Repeat tech row, then fill in plant specifics in each row
        sitesDf = pd.concat([re]*len(sites),ignore_index=True)
        sitesDf['SiteName'] = sites
        sitesDf['PlantType'] = pt
        #Get lat/lon
        txt = sitesDf['SiteName'].str.split('lat',expand=True)[1]
        sitesDf[['Latitude','Longitude']] = txt.str.split('lon',expand=True).astype(float)
        newTechsCE.drop(re.index,inplace=True)
        sitesDfList.append(sitesDf)

    #Combine wind & solar rows into df, then map to regions
    sitesDf = pd.concat(sitesDfList)
    sitesDf = sitesDf.drop('region',axis=1)
    if sitesDf['Longitude'].max() > 180: 
        sitesDf['Longitude'] -= 360
        sitesDf = assignGensToPRegions(sitesDf,pRegionShapes)
        sitesDf['Longitude'] += 360
    else:
        sitesDf = assignGensToPRegions(sitesDf,pRegionShapes)

    #Drop bottom percentile of wind and solar sites by annual average CF in each region
    allSitesToDrop = list()
    for region in sitesDf['region'].unique():
        regionGens = sitesDf.loc[sitesDf['region']==region]
        for re in ['wind','solar','offwind']:
            regionFTGens = regionGens.loc[regionGens['SiteName'].str.contains(re)]['SiteName'].values
            meanCfs = newCfs[regionFTGens].mean()
            meanCfs.sort_values(inplace=True) #ascending
            sitesToDrop = meanCfs.iloc[:int(meanCfs.shape[0]*reDownFactor)].index.values
            allSitesToDrop.extend(sitesToDrop)
    newCfsOrig.drop(allSitesToDrop,axis=1,inplace=True)
    sitesDf = sitesDf.loc[~sitesDf['SiteName'].isin(allSitesToDrop)]

    #Add remaining WS sites onto newTechsCE
    newTechsCE = pd.concat([newTechsCE,sitesDf],ignore_index=True)
    newTechsCE.reset_index(inplace=True,drop=True)

    #Create GAMS symbol as plant type + location + region
    newTechsCE['GAMS Symbol'] = newTechsCE['PlantType'] + 'lat' + np.round(newTechsCE['Latitude'],3).astype(str) + 'lon' + np.round(newTechsCE['Longitude'],3).astype(str) + newTechsCE['region']
    
    #Relabel new CF columns & maximum capacities to match GAMS symbols; use dictionary of original site name (lacking region) to GAMS symbol (includes region)
    reRows = newTechsCE.loc[newTechsCE['ThermalOrRenewableOrStorage']=='renewable']
    reRows.index = reRows['SiteName']
    gamsDict = reRows['GAMS Symbol'].to_dict()
    if type(newCfsOrig.columns)==pd.core.indexes.multi.MultiIndex:
        newCfsOrig = newCfsOrig.T.loc[:,[k for k in gamsDict],:].rename(gamsDict,level='locs')
        newCfsOrig = newCfsOrig.T
    else:
        newCfsOrig = newCfsOrig[[k for k in gamsDict]].rename(gamsDict,axis=1)
    
    #For each unique site (site name), attach max capacity to GAMS symbol and remove original site name
    for sn in gamsDict:
        maxCapPerTech[gamsDict[sn]] = maxCapPerTech[sn]
        maxCapPerTech.pop(sn)

    #Filter out wind and solar sites that were contained in original CFs but don't actually belong in region (e.g., offshore sites) and sites that were dropped due to low CFs
    toRemove = list()
    for k in maxCapPerTech:
        if 'wind' in k or 'solar' in k: #lower case is crucial here. Ex: wind site names in original CFs are 'windlatXXlonXX'; replaced by 'Onshore WindlatXXlonXX'
            toRemove.append(k)
    for k in toRemove: maxCapPerTech.pop(k)

    #Drop SiteName column
    newTechsCE.drop('SiteName',axis=1,inplace=True)

    return newTechsCE,newCfsOrig,maxCapPerTech
