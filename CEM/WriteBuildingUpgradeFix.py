#Michael Craig
#April 24, 2020
#Write .gms file w/ CE text for time blocks and associated constraints

import os

def writeBuildingUpgradeFix(forceBldgUpgrade,gamsFileDir):
	if forceBldgUpgrade != 0: 
		txt = '\nvNBldgUpgrade.fx[\'upgrade'+str(forceBldgUpgrade)+'\']=1;\n'
	else:
		txt = '\n'
	g = open(os.path.join(gamsFileDir,'CEBldgUpgradeFix.gms'),'w')
	g.write(txt)
	g.close()