#Shell script for running macro CEM

import sys,os
from MacroCEM import macroCEM

#Set working directory to location of this script
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

#Process inputs and call master function
inputData = sys.argv[1:] #exclude 1st item (script name)
wsGenFracOfDemand = int(inputData[0]) #integer giving % of demand that must be supplied by wind and solar generation by final year
wsGenYr = int(inputData[1]) #year in which to fully enforce wsGenFracOfDemand
co2Cap = int(inputData[2]) #integer giving % of CO2 emissions in final year relative to first year
co2CapYr = int(inputData[3]) #year in which to fully enforce co2 cap
forceBldgUpgrade = int(inputData[4]) #which building upgrade package to force; if 0, no upgrade occurs
climateScenario=inputData[5]

macroCEM('NY',climateScenario,wsGenFracOfDemand,wsGenYr,co2Cap,co2CapYr,forceBldgUpgrade)