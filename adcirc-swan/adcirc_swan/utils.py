

from math import cos, sin, asin, sqrt, radians
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import os
#import adcirc
#import swan
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from IPython.display import HTML
import plotly.graph_objs as go
import plotly.offline as po
import pandas as pd

def output_files(root_dir):
    f61 = os.path.join(root_dir,'fort.61.nc')
    f62 = os.path.join(root_dir,'fort.62.nc')
    f63 = os.path.join(root_dir,'fort.63.nc')
    f64 = os.path.join(root_dir,'fort.64.nc')
    f73 = os.path.join(root_dir,'fort.73.nc')
    f74 = os.path.join(root_dir,'fort.74.nc')
    rad = os.path.join(root_dir,'rads.64.nc')
    mrad= os.path.join(root_dir,'maxrs.63.nc')
    max_f63 = os.path.join(root_dir,'maxele.63.nc')
    max_wv63 = os.path.join(root_dir,'maxwvel.63.nc')
    max_v63  = os.path.join(root_dir,'maxvel.63.nc')
    minpr63  = os.path.join(root_dir,'minpr.63.nc')
    swan_hs = os.path.join(root_dir,'swan_HS.63.nc')
    max_hs = os.path.join(root_dir,'swan_HS_max.63.nc')
    swan_dir = os.path.join(root_dir,'swan_DIR.63.nc')
    max_dir = os.path.join(root_dir,'swan_DIR_max.63.nc')
    files = [f61,f62,f63,f64,f73,f74,max_f63,max_wv63,max_v63,minpr63,swan_hs,max_hs,swan_dir,max_dir,rad,mrad]
    return files

def timeseries(nc_data,start,freq,station,name):
    data = nc_data['zeta'][:,station]*3.28084#106378]
    table = pd.DataFrame(data)
    date  = pd.date_range(start=start,periods=int(len(table)),freq=freq)
    table.insert(0,'Date Time',date)
    table = table.rename(columns={0:name})
    return table

def plot_timeseries(table,start,end,station):
    pred = pd.read_csv('https://tidesandcurrents.noaa.gov/api/datagetter?product=predictions&application=NOS.COOPS.TAC.WL&begin_date='+start+'&end_date='+end+'&datum=MSL&station='+station+'&time_zone=GMT&units=english&interval=6&format=csv')
    obs  = pd.read_csv('https://tidesandcurrents.noaa.gov/api/datagetter?product=water_level&application=NOS.COOPS.TAC.WL&begin_date='+start+'&end_date='+end+'&datum=MSL&station='+station+'&time_zone=GMT&units=english&format=csv')
    trace = go.Scatter(x = table['Date Time'],y = table[' Prediction'],
                name = 'ADCIRC',mode = 'lines',
                line = dict(
                    color = ('rgb(204, 0, 153)')))
    trace2 = go.Scatter(x = pred['Date Time'],y = pred[' Prediction'],
                    name = 'NOAA Prediction',mode = 'lines',visible = 'legendonly',
                    line = dict(
                        color = ('rgb(255, 153, 51)')))
    trace3 = go.Scatter(x = obs['Date Time'],y = obs[' Water Level'],
                    name = 'NOAA Observation',mode = 'lines',
                    line = dict(
                        color = ('rgb(100, 100, 53)')))
    data = [trace,trace2,trace3]
    return data
def plot_tidal(table,start,end,station):
    pred = pd.read_csv('https://tidesandcurrents.noaa.gov/api/datagetter?product=predictions&application=NOS.COOPS.TAC.WL&begin_date='+start+'&end_date='+end+'&datum=MSL&station='+station+'&time_zone=GMT&units=english&interval=6&format=csv')
    obs  = pd.read_csv('https://tidesandcurrents.noaa.gov/api/datagetter?product=water_level&application=NOS.COOPS.TAC.WL&begin_date='+start+'&end_date='+end+'&datum=MSL&station='+station+'&time_zone=GMT&units=english&format=csv')
    table = table.set_index(pd.DatetimeIndex(table['Date Time']))
    pred = pred.set_index(pd.DatetimeIndex(pred['Date Time']))
    trace = go.Scatter(x = table['Date Time'],y = table[' Prediction'],
                name = 'ADCIRC',mode = 'lines',
                line = dict(
                    color = ('rgb(204, 0, 153)')))
    trace2 = go.Scatter(x = pred['Date Time'],y = pred[' Prediction'],
                    name = 'NOAA Prediction',mode = 'lines',
                    line = dict(
                        color = ('rgb(100, 100, 153)')))
    trace3 = go.Scatter(x = obs['Date Time'],y = obs[' Water Level'],
                    name = 'NOAA Observation',mode = 'lines',
                    line = dict(
                        color = ('rgb(100, 100, 53)')))
    data = [trace,trace2]
    return data

def find_node_ak(fmax,obs_lat,obs_lon):
    min_distance = None
    best_index = 0
    x = fmax.variables['x'][:]
    y = fmax.variables['y'][:]
    for i in range(len(x)):
        current_distance = (float(y[i]) - obs_lat)**2 + (float(x[i]) - obs_lon)**2
        if min_distance is None or current_distance < min_distance:
            best_index = i
            min_distance = current_distance
    #print("best_index:{} ".format(best_index))
    return best_index

def calc_distance(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    km = 6371 * c
    return km

def layout(title,xaxis,yaxis):
    layout = go.Layout(dict(title=title),xaxis = dict(title = xaxis),
                   yaxis = dict(title = yaxis),legend= dict(orientation="h"),
                   font = dict(color = 'rgb(0,0,0)'),paper_bgcolor = 'rgb(255,255,255)',
                   plot_bgcolor = 'rgb(255,255,255)')
    return layout
    
def usgs_gauges(start,end,station,vdatum:str='NAVD88'):
    names=['owner','id','Date_time','zone','water level(ft) NAVD88','code']
    url = f'https://nwis.waterdata.usgs.gov/nwis/uv?cb_62620=on&format=rdb&site_no={station}&period=&begin_date={start}&end_date={end}'
    #site = requests.get(url = url).content.decode()
    usgs = pd.read_csv(url,sep='\t',skiprows=28,names=names,header=None)
    if 'NAVD88' not in vdatum:
        factor = 0.37
        usgs.insert(4,'water_level(ft) at MSL',usgs['water level(ft) NAVD88']-factor)
    return usgs

def plot_usgs(table,usgs):
    
    trace = go.Scatter(x = table['Date Time'],y = table[' Prediction'],
                name = 'ADCIRC',mode = 'lines',
                line = dict(color = ('rgb(204, 0, 153)')))

    trace2 = go.Scatter(x = pd.to_datetime(usgs['Date_time']),y = usgs['water_level(ft) at MSL'],
        name = 'USGS observations',mode = 'lines',
        line = dict(color = ('rgb(100, 100, 53)')))

    data = [trace,trace2]
    return data



def force2owi(file,dx,dy,dt,ilon,ilat,start,end,swlat,swlon,header,grib,lat1,lat2,lon1,lon2):
    if file == 'fort.222':
        u = grib.data_vars['u10'].values[:]
        v = grib.data_vars['v10'].values[:]
    elif file == 'fort.221':
        prmsl = grib.data_vars['msl'].values[:]
    elif file == 'fort.225':
        ice = grib.data_vars['siconc'].values[:]
    else:
        raise TypeError('not correct forcing format')
    
    
    lat = grib.coords['latitude'].values
    lon= grib.coords['longitude'].values
    idx = np.where((lon1 >= lon)&(lon<=lon2))[0]
    idy = np.where((lat1 >= lat)&(lat<=lat2))[0]
    with open(file,'w') as fin:
        fin.write(header)
        for t in range(len(dt[:])):
            param ='iLat= {}iLong= {}DX={:.4f}DY={:.4f}SWLat={:.4f}SWLon={:.4f}DT={} \n'.format(ilat,ilon,dx,dy,swlat,swlon,str(dt[t].strftime('%Y%m%d%H%M')))
            fin.write(param)
            i = 0
            data = []
            d1,d2 = [],[] 
            for y in idy[:]:
                for x in idx[:]:
                    if file == 'fort.222':
                        windx = u[t,y,x]
                        windy = v[t,y,x]
                        d1.append(' {:{width}.{prec}f}'.format(windx,width=9,prec=4))
                        d2.append(' {:{width}.{prec}f}'.format(windy,width=9,prec=4))
                    elif file == 'fort.221':
                        press = prmsl[t,y,x]
                        data.append(' {:{width}.{prec}f}'.format(press,width=9,prec=4))
                    elif file == 'fort.225':
                        sea = ice[t,y,x]
                        if np.isnan(sea):
                            sea = -1.0
                        else:
                            sea =sea * 100
                        data.append(' {:{width}.{prec}f}'.format(sea,width=9,prec=4))
            i = 0
            if file == 'fort.222':
                for d in range(len(d1)):
                    if (i == 7) :
                        fin.write(d1[d]+'\n')
                        i=0
                    elif (d+1>=len(d1)):
                        fin.write(d2[d]+'\n')
                        i = 0
                        for d in range(len(d2)):
                            if (i == 7) :
                                fin.write(d2[d]+'\n')
                                i=0
                            elif (d+1>=len(d2)):
                                fin.write(d2[d]+'\n')
                            else:
                                fin.write(d2[d])
                                i +=1
                    else:
                        fin.write(d1[d])
                        i +=1
            else:
                for d in range(len(data)):
                    if (i == 7) :
                        fin.write(data[d]+'\n')
                        i=0
                    elif (d+1>=len(data)):
                        fin.write(data[d]+'\n')
                    else:
                        fin.write(data[d])
                        i +=1
    return print('generated owi formatted forcings')