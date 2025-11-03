import pandas as pd
from sys import argv
import os
import json

config = pd.read_table(argv[1],index_col=0,  names=['path','description'])[['path']].fillna('__NaNPATH__')
param = pd.read_table(config.at['param', 'path'], index_col=0 ,names=['value'])

for index in config.index:
    if config.at[index, 'path'] == '__NaNPATH__':
        if index == 'outdir':
            config.at[index, 'path'] = config.at['genecat', 'path'] + os.sep + 'Bin_SB' + 'intra_phylo'
        elif index == 'tempdir':
            config.at[index, 'path'] = config.at['outdir', 'path'] + os.sep + "MAGRecTK"
        elif index == 'mgslist':
            config.at[index, 'path'] = 'ALL'
    
    elif index == 'param':
        config.at[index, 'path'] = os.path.dirname(config.at[index, 'path']) + os.sep + os.path.splitext(os.path.basename(config.at[index, 'path']))[0] + '.json'
    elif index == 'tempdir':
        config.at[index, 'path'] = config.at[index, 'path'] + os.sep + "MAGRecTK"


outpath = config.at['outdir', 'path']

with open(f"{os.path.dirname(argv[1])}{os.sep}config.json", "w") as f:
    json.dump(dict(zip(config.index, config['path'])),f,indent=4)

with open(f"{os.path.dirname(argv[1])}{os.sep}param.json", "w") as f:
    json.dump(dict(zip(param.index, param['value'])),f,indent=4)

