#Michael Craig
#October 4, 2016
#Functions update input fleets' fuel prices for given year

from SetupGeneratorFleet import *
from AddRegionalGasPricesFromJill import addRegionalGasPricesFromJill

def updateFuelPricesAndCosts(fleet,currYear,fuelPrices,regCostFrac):
    if currYear > 2050: currYear = 2050
    fleet = addFuelPrices(fleet,currYear,fuelPrices)
    fleet = calcOpCost(fleet)
    fleet['RegOfferCost($/MW)'] = regCostFrac*fleet['OpCost($/MWh)']*fleet['RegOfferElig']

    fleet = addRegionalGasPricesFromJill(fleet) 
    return fleet