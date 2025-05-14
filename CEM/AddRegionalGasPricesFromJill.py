import os, pandas as pd

def addRegionalGasPricesFromJill(genFleet):
	prices = pd.read_csv(os.path.join('Data','GasPricesJill.csv'))
	existingPrices = prices.loc[prices['Zone'].str.contains('Existing')]
	for r in ['p127','p128']:
		prices = existingPrices.loc[existingPrices['Zone'].str.contains(r[1:])] #region name is 127 or 128, not p127 or p128
		price = prices.loc[prices['Month'].isin([5,6,7,8,9])]['Avg Fuel Price 2010-2023 ($2023)'].mean() #use average summer price
		price = price*294/304.7 #update for CPI to go from 2023 to 2022 dollars
		genFleet.loc[(genFleet['FuelType']=='Natural Gas')&(genFleet['region']==r),'FuelPrice($/MMBtu)'] = price
	return genFleet
