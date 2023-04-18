
from math import cos, sin, asin, sqrt, radians
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import os; import xarray as xr
import datetime
import matplotlib.pyplot as plt
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from IPython.display import HTML
import plotly.graph_objs as go
import plotly.offline as po
import pandas as pd;from pyproj import Proj, transform
import requests
import numpy as np
from scipy import ndimage as nd
from datetime import datetime
from bs4 import BeautifulSoup

def find_columns(data):
    data2 = []
    for f in data.split(' '):
        if f != '':
            data2.append(f)   
    return data2

def read_fort51(fort51:str,amplitude:dict,phase:dict):
    keys = list(amplitude.keys())
    with open(fort51,'r') as fin:
        lines0 = fin.readlines()    
        lines=lines0
        ke = 0
        for i in range(10,len(lines)):
            vals_1 = find_columns(lines[i].strip()) 
            if (lines[i].strip().isnumeric()):
                ke=1
            else:            
                amp1,pha1=float(find_columns(lines[i].strip())[0]),float(find_columns(lines[i].strip())[1])
                if ke:
                    amplitude[keys[ke-1]].append(amp1)
                    phase[keys[ke-1]].append(pha1)
                    ke+=1
    return amplitude,phase

def find_ahps(geo,obs_lat,obs_lon):
    min_distance = 30
    best = 0
    #x1,y1 = x.shape
    i = 0
    for g in geo:
        current_distance = (g.y - obs_lat)**2 + (g.x - obs_lon)**2
        if min_distance is None or current_distance < min_distance:
            best = i
            min_distance = current_distance
        i+=1
    #print("best_index:{} ".format(best_index))
    return best

def evaluate_harmonics(stations:pd.DataFrame,constituents:list,amplitudes,phas):
    adcirc_tidal = {}
    NOAA_Tidal_database = {}
    for station in stations[:]:
        url = r'https://tidesandcurrents.noaa.gov/harcon.html?unit=0&timezone=0&id={}'.format(station.split(':')[1])
        r = requests.get(url)
        data = r.text
        soup = BeautifulSoup(data, "lxml")
        #---Data
        nm = soup.find_all('h3')
        name_stn=str(nm)[5:-6].split(', ')[-1]
        data = soup.find_all('table')[0] 
        data_rows = data.find_all('tr')[1:]
        name1 = []
        amp1 = []
        phase1 = []
        for row in data_rows:
            d = row.find_all('td')
            if d[1].get_text() in constituents:
                name = d[1].get_text() 
                amp = d[2].get_text()
                phase = d[3].get_text()
                speed = d[4].get_text()
                fullname = d[5].get_text()
                name1.append(name)
                amp1.append(amp)
                phase1.append(phase)
        NOAA_Tidal_database.update({f'Station:{name_stn} ID:{station}':{'amp':{},'phas':{}}})
        for n in range(len(name1)):
            NOAA_Tidal_database[f'Station:{name_stn} ID:{station}']['amp'].update({name1[n]:float(amp1[n])})
            NOAA_Tidal_database[f'Station:{name_stn} ID:{station}']['phas'].update({name1[n]:float(phase1[n])})

        adcirc_tidal.update({f'{station}':{'amp':{},'phas':{}}})
        for j in constituents:
            adcirc_tidal[f'{station}']['amp'].update({j:amplitudes[j][stations.index(station)]})
            adcirc_tidal[f'{station}']['phas'].update({j:phas[j][stations.index(station)]})
    return adcirc_tidal,NOAA_Tidal_database

def noaa_harmonics(stations:pd.DataFrame,constituents:list):
    NOAA_Tidal_database = {}
    for station in stations[:]:
        url = r'https://tidesandcurrents.noaa.gov/harcon.html?unit=0&timezone=0&id={}'.format(station.split(':')[1])
        r = requests.get(url)
        data = r.text
        soup = BeautifulSoup(data, "lxml")
        #---Data
        nm = soup.find_all('h3')
        name_stn=str(nm)[5:-6].split(', ')[-1]
        data = soup.find_all('table')[0] 
        data_rows = data.find_all('tr')[1:]
        name1 = []
        amp1 = []
        phase1 = []
        for row in data_rows:
            d = row.find_all('td')
            if d[1].get_text() in constituents:
                name = d[1].get_text() 
                amp = d[2].get_text()
                phase = d[3].get_text()
                speed = d[4].get_text()
                fullname = d[5].get_text()
                name1.append(name)
                amp1.append(amp)
                phase1.append(phase)
        NOAA_Tidal_database.update({f'{station}':{'amp':{},'phas':{}}})
        for n in range(len(name1)):
            NOAA_Tidal_database[f'{station}']['amp'].update({name1[n]:float(amp1[n])})
            NOAA_Tidal_database[f'{station}']['phas'].update({name1[n]:float(phase1[n])})

    return NOAA_Tidal_database



def fill(data, invalid=None):
    """
    Replace the value of invalid 'data' cells (indicated by 'invalid') 
    by the value of the nearest valid data cell

    Input
        data:    numpy array of any dimension
        invalid: a binary array of same shape as 'data'. True cells set where data
                 value should be replaced.
                 If None (default), use: invalid  = np.isnan(data)

    Output: 
        Return a filled array. 
    """
    #import numpy as np
    #import scipy.ndimage as nd

    if invalid is None: invalid = np.isnan(data)

    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    return data[tuple(ind)]

def sat_seaice2owi(header:str,files:list,inProj:Proj,\
    outProj:Proj,lat1:float,lat2:float,lon1:float,lon2:float,\
    sea_ice:str='siconc',file:str='fort.225'):
    '''
    this function is a derivative of the force2owi function.
    it is expecting satellite sea ice data. meaning the timestep of the data is seperated by files
    '''

    with open(file,'w') as fin:
        fin.write(header)
        for f in files:
            # loop and read in data from each file
            data= xr.open_dataset(f)
            
            lat = data.coords['ygrid'].values
            lon= data.coords['xgrid'].values
            xgrid,ygrid = np.meshgrid(lon,lat)
            xt,yt = np.reshape(xgrid,xgrid.shape[0]*xgrid.shape[1]),np.reshape(ygrid,ygrid.shape[0]*ygrid.shape[1])
            x2,y2 = transform(inProj,outProj,xt,yt)
            #xgrid2,ygrid2 = np.meshgrid(x2,y2)
            ice = np.reshape(data.data_vars[sea_ice].values[0,:,:],x2.shape[0])
            dt = pd.to_datetime(f.name.split('_')[4])

            # finding all the points within the bounding box
            ind = np.where((lon1+360 <= x2)&(x2<=lon2+360)&(lat1 <= y2)&(y2<=lat2))[0]
            idx = np.where((lon1+360 <= x2)&(x2<=lon2+360))[0]
            idy = np.where((lat1 <= y2)&(y2<=lat2))[0]
            dx = np.abs(x2[0]-x2[1])
            dy = np.abs(y2[0]-y2[1])
            # writing parameters for each time step in file
            param ='iLat={:4d}iLong={:4d}DX={:6.4f}DY={:6.4f}SWLat={:8.5f}SWLon={:8.3f}DT={:12s}\n'.\
            format(len(idy),len(idx),dx,dy,lat1,lon1,str(dt.strftime('%Y%m%d%H%M')))#,width=8,prec=5)
            fin.write(param)
            i = 0
            content = []
            d1,d2 = [],[] 
            #ind = np.where(idy==idx)[0]
            for y in ind[:]:
                if file == 'fort.225':
                    if ice[y]>2:
                        sea = 0
                    else:
                        sea = ice[y]
                    if np.isnan(sea):
                        sea = -1.0
                    else:
                        sea =sea #* 100
                    content.append(' {:{width}.{prec}f}'.format(sea,width=9,prec=4))
            for d in range(len(content)):
                if (i == 7) : 
                    fin.write(content[d]+'\n')
                    i=0
                elif (d+1>=len(content)):
                    fin.write(content[d]+'\n')
                else:
                    fin.write(content[d])
                    i +=1

    print(f'amount of data {len(idx)*len(idy)} data is {len(data)}')
    return print('generated owi formatted forcings')


def force2owi(file:str,dx:float,dy:float,dt:np.array,ilon:int,ilat:int,\
                start:datetime,end:datetime,swlat:int,swlon:int,header:str,\
                grib:xr.core.dataset.Dataset,lat1:float,lat2:float,lon1:float,lon2:float,\
                uwind:str='u10',vwind:str='v10',psl:str='msl',sea_ice:str='siconc',\
                forcing:str='era5',interval:int=1,latcfs:np.array=None,loncfs:np.array=None):
    '''
    inputs:
        file= adcirc owi formatted file (e.g fort.221, fort.222, fort.225)
        dx = x resolution
        dy = y resolution
        dt = length of time
        ilon= length of array of longitudes
        ilat= length array of latitudes
        start= datetime of start
        end = datetime of end
        swlat = south west lat coord
        swlon = south west lon coord
        header = string for the header of owi files
        (e.g. 'Oceanweather WIN/PRE Format                            {:12s}   \
        {:12s}\n'.format(str(start.strftime('%Y%m%d%H')),str(end.strftime('%Y%m%d%H'))))
        lat1,lat2,lon1,lon2 = bounding box of the model domain
        uwind = name of uwind variable in grib file
        vwind = name of vwind variable in grib file
        psl = name of pressure at msl in grib file
        sea_ice = name of sea ice in grib file
        forcing = forcing type. current approved types era5 and cfs
        interval= time interval meant to change the hourly timesteps in cfs
        latcfs = when using cfs you need to pass the wind lats too so that the pressure can be interpolated to the same grid
        loncfs = when using cfs you need to pass the wind lons too so that the pressure can be interpolated to the same grid

    '''
    lat = grib.coords['latitude'].values
    lon= grib.coords['longitude'].values
    xgrid,ygrid = np.meshgrid(lon,lat)
    
    if file == 'fort.222':
        u = grib.data_vars[uwind].values[:]
        v = grib.data_vars[vwind].values[:]
        data = grib.data_vars[vwind].values[:]
    elif file == 'fort.221':  
        prmsl = grib.data_vars[psl].values[:]
        data  = grib.data_vars[psl].values[:]
    elif file == 'fort.225':
        ice = grib.data_vars[sea_ice].values[:]
        data= grib.data_vars[sea_ice].values[:]
    else:
        raise TypeError('not correct forcing format')

    idx = np.where((lon1+360 <= lon)&(lon<=lon2+360))[0]
    idy = np.where((lat1 <= lat)&(lat<=lat2))[0]
    if (len(data.shape) > 3):
        step = True
        print('cfs 4 array')
        dt2 = dt
        steps = 7
        dt = grib.coords['valid_time'].values[:len(dt)]
    
    '''
    CFS forcing assumed 
    arrays need to be [hour,date,lat,lon]
    '''
    xgrid,ygrid = np.meshgrid(lon,lat)
    
    with open(file,'w') as fin:
        fin.write(header)
        if forcing=='cfs':
            for count,t in enumerate(dt):
                iterate = 0
                ts = pd.to_datetime(dt[count])

                for tt in np.arange(0,len(t)-1,interval):
                    delta = pd.to_datetime(t[tt])
                    param ='iLat={:4d}iLong={:4d}DX={:6.4f}DY={:6.4f}SWLat={:8.5f}SWLon={:8.3f}DT={:12s}\n'.\
                    format(ilat,ilon,dx,dy,swlat,swlon,str(delta.strftime('%Y%m%d%H%M')))#,width=8,prec=5)
                    fin.write(param)
                    i = 0
                    data = []
                    d1,d2 = [],[]
                    if file == 'fort.221':
                        if len(latcfs)<2:
                            print('Need wind lat/lons to interpolate pressure to same grid')
                            break
                        else:
                            idx = np.where((lon1+360 <= loncfs)&(loncfs<=lon2+360))[0]
                            idy = np.where((lat1 <= latcfs)&(latcfs<=lat2))[0]
                            xgrid2,ygrid2 = np.meshgrid(loncfs,latcfs)
                            press = fill(grib.data_vars[psl].values[count,tt,:,:])
                            latr = np.reshape(ygrid,(ygrid.shape[0]*ygrid.shape[1]))
                            lonr = np.reshape(xgrid,(xgrid.shape[0]*xgrid.shape[1]))
                            press2= np.reshape(press,(press.shape[0]*press.shape[1]))
                            press3= scipy.interpolate.griddata((lonr,latr),press2,(xgrid2,ygrid2),method='nearest')
                    for y in reversed(idy[:]):
                        for x in idx[:]:
                            if file == 'fort.222':
                                windx = u[count,tt,y,x]
                                windy = v[count,tt,y,x]
                                d1.append(' {:{width}.{prec}f}'.format(windx,width=9,prec=4))
                                d2.append(' {:{width}.{prec}f}'.format(windy,width=9,prec=4))
                            elif file == 'fort.221':
                                prmsl = press3[y,x]/100
                                data.append(' {:{width}.{prec}f}'.format(prmsl,width=9,prec=4))
                            elif file == 'fort.225':
                                sea = ice[count,tt,y,x]
                                if np.isnan(sea):
                                    sea = -1.0
                                else:
                                    sea =sea #* 100
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
                    iterate += 1 
                                
        else:
            
            for t in range(len(dt)):
                param ='iLat={:4d}iLong={:4d}DX={:6.4f}DY={:6.4f}SWLat={:8.5f}SWLon={:8.3f}DT={:12s}\n'.\
                format(ilat,ilon,dx,dy,swlat,swlon,str(dt[t].strftime('%Y%m%d%H%M')))#,width=8,prec=5)
                fin.write(param)
                i = 0
                data = []
                d1,d2 = [],[] 
                for y in reversed(idy[:]):
                    for x in idx[:]:
                        if file == 'fort.222':
                            windx = u[t,y,x]
                            windy = v[t,y,x]
                            d1.append(' {:{width}.{prec}f}'.format(windx,width=9,prec=4))
                            d2.append(' {:{width}.{prec}f}'.format(windy,width=9,prec=4))
                        elif file == 'fort.221':
                            press = prmsl[t,y,x]/100
                            data.append(' {:{width}.{prec}f}'.format(press,width=9,prec=4))
                        elif file == 'fort.225':
                            sea = ice[t,y,x]
                            if np.isnan(sea):
                                sea = -1.0
                            else:
                                sea =sea #* 100
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
        #fin.write('\n')
    
    print(f'amount of data {len(idx)*len(idy)} data is {len(data)}')
    return print('generated owi formatted forcings')

'''
for modifying arctic meshes
with open(root / 'fort_v2.14','w') as fout:
    with open(root / 'fort.14','r') as fin:
        lines = fin.readlines()
        llimit = int(lines[1].strip().split(' ')[-1])
        track = 0
        for line in lines[:]:
            if (track >= 2) & (track <= llimit):
                data = line.split(' ')
                t = 0
                for i in range(len(data)):

                    try:
                        num = float(data[i])
                        if (num > 0) and (t==1):
                            num = num-360
                            data[i] = str(num)
                            print(data)
                        t+=1
                    except:
                        pass
                line = ' '.join(data)
            else:
                pass
            track+=1
            fout.write(line)
            
        
        #try:
        #    digits = [i if float(i) for i in data]



# wavewatch fort file
root = pl.Path('/Users/tmiesse/work/FHRL/arctic/ww3')
with open(root / 'global_mesh2' /'fort.msh','w') as fout:
    with open(root / 'global_mesh' / 'fort.msh','r') as fin:
        lines = fin.readlines()
        llimit = int(lines[4].strip().split(' ')[-1])
        track = 0
        for line in lines[:]:
            if (track >= 5) & (track <= llimit+4):
                data = line.split(' ')
                t = 0
                for i in range(len(data)):
                    try:
                        num = float(data[i])
                        if (num > 0) and (t==1):
                            num = num-360
                            data[i] = str(num)
                            print(data)
                        t+=1
                    except:
                        pass
                line = ' '.join(data)
            else:
                pass
            track+=1
            fout.write(line)
'''