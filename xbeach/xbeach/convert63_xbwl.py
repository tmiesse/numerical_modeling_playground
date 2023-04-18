
import pandas as pd
import netCDF4 as nc4
import numpy as np
import os
#import oct2py
#from oct2py import octave
import pathlib as pl
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt

root = pl.Path(r'/home/admin/neptune_work/realtime/xbeach/Magothy_Bay')
out = root / 'input'



f63 = nc4.Dataset(str(root / 'tide_from_arslaan' /'fort.63.nc'))
def timeseries(nc_data,station):
    data = nc_data['zeta'][:,station]
    table = pd.DataFrame({'water_level':data})
    return table

def write_tide(root_dir,time,front,back):
    root_dir = pl.Path(root_dir)
    file = root_dir / 'tide.txt'
    data = []
    for i in range(0,len(time)):
        #if 'nan' in float(front[i]):
         #   front[i] = 0
         #  back[i] = 0
        data.append('    {:.4e}    {:.4e}'.format(float(time[i]),float(front[i])*1.5) + '\n')
        #data.append('    {:.4e}    {:.4e}    {:.4e}'.format(float(time[i]),float(front[i])*1.5,float(back[i])*1.5) + '\n')
    with open(file,'w') as fin:
        fin.writelines(data)
    return

time = np.arange(0,302400,3600)
wl = timeseries(f63,286406)




write_tide(out,time,wl['water_level'])





