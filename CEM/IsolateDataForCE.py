##### ISOLATE DATA FOR USE IN CE
#Inputs: , any # of dfs
#Outputs: any # of dfs w/ just hours used in CE
def isolateDataInCEHours(hoursForCE,*args):
    return [df.loc[hoursForCE.index] for df in args]

    # [hydroGenCE] = isolateDataInCEBlocks(hoursForCE,hydroGen)

def isolateDataInCEBlocks(hoursForCE,*args):
    blockDataAll = list()
    for df in args:
    	dfInCE = df.loc[hoursForCE.index]
    	origTotalGen = dfInCE.sum().sum()
    	dfInCE['block'] = hoursForCE
    	blockData = dfInCE.groupby(['block']).sum()
    	if dfInCE.sum().sum()>1e-3: assert((origTotalGen - blockData.sum().sum())<=origTotalGen*.0001)
    	blockDataAll.append(blockData)
    return blockDataAll