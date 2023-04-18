%%      
% This code generates parameters for xbeach in 1d. There are 4 transects to
% choose from that represent the real transects out in the field.
% It is similar to the 2d parameter codes created and it uses similar
% settings such as:
%           -   waves
%           -   tide
%           -   vegetation
%           -   wind


run '/Users/tmiesse/work/libraries/delft/oetsettings.m';
clear all; clc;

paths = dir('/Users/tmiesse/work/FHRL/eeslr/modelling/xbeach/dan_study/new_grid/response2reviews/seasonality3/spart/');

outpath = '/Users/tmiesse/work/FHRL/eeslr/modelling/xbeach/dan_study/new_grid/response2reviews/seasonality3/spart/output/';

%%      Define Path Directory
% This helps with keeping everything in order instead of writing the file
% path multiple times

%path = '/Users/tmiesse/work/FHRL/eeslr/modelling/xbeach/dan_study/new_grid/high_water/spart';
%names = {'/'};
root = '/Users/tmiesse/work/FHRL/eeslr/modelling/xbeach/dan_study/new_grid/response2reviews/review1/';
xbo  =xb_read_output(root);%,'vars', {'H_mean'});
xb_dat2nc(root,[root,'event1_response.nc'],'vars',{'H','H_mean','zs_mean',...
                'zb_mean','zs','zb'});
            
%%
for i=1:length(paths)
    root = [paths(i).folder,'/',paths(i).name];
    subpath = dir(root);
    for ii=1:length(subpath)
        root2 = [subpath(ii).folder,'/',subpath(ii).name];
        if isfile([root2,'/','dims.dat'])
            %if isfile([root,[outpath,paths(i).name,'.nc']])
            %    ;
            %else
                xbo  =xb_read_output(root2);%,'vars', {'H_mean'});
                xb_dat2nc(root2,[outpath,paths(i).name,'_',subpath(ii).name,'.nc'],'vars',{'H','H_mean','zs_mean',...
                    'zb_mean','zs','zb'});
            %end
        end
    end
end

%%
% B. Model
dxmin           = 1;
dymin           = 1;
outputformat    = 'netcdf';


%%      Generate the waves 
path = '/Users/tmiesse/work/FHRL/eeslr/modelling/xbeach/assateague/spectrum/SPECTRUM_RESULTS/IRENE/';
sp2 = xb_sp22xb([path,'ASSA.txt']);

% xb_wave = xb_generate_waves('Hm0',0.73,'gammajsp',3.3,'s',10,...
%                             'mainang',270,'fnyq',.45,'Tp',4.45);   
    
%%    Generate the tide

 %xb_tide=xb_generate_tide('time',tide_time,'front',tide_height);    

%%    Generate the settings for xbeach to follow

xb_set=xb_generate_settings('outputformat',outputformat,... 
        'thetamin',180,'thetamax',360,'dtheta',90,'dtheta_s',5,...
        'instat','swan','morfac', 1,'posdwn',-1,'avalanching',0,...
        'morstart', 0,'CFL', 0.9,'front', 'abs_2d','random',0,...
        'back', 'abs_2d','left','neumann','right','neumann','mpiboundary','auto',...
        'thetanaut', 1,'zs0',0.6,'single_dir',0,...
        'tstop', 3600,'tstart', 0,...
        'tint', 3600,'tintm',360,'tintg',360,'epsi',-1,'facua',0.300,'bedfriction', 'manning',...
        'meanvar',{'zb','zs', 'H','u','v','sigm','Cd'} ,... 
        'globalvar',{'zb', 'zs','H','u','v','sigm','Cd'});


%%      Generate the bathy for Tansect 

%load '/Users/tmiesse/work/dewberry/vb_terrace/modeling/gis/delft_inputs/dem_pts/vb_terrace.mat';


path = '/Users/tmiesse/Downloads/';
file = 'temp.csv';
%path = '/Users/tmiesse/work/FHRL/eeslr/field/adcp/adcp22/';
%file = 'ADCP202.csv';
x=importdata([path,file],',');


dx = x.data(:,8);

dz = x.data(:,7);

figure;pcolor(x,y,z);caxis([-2 2]);shading interp; colormap jet

%%
xb2 =xb.data(1).value;
y = xb2.data(17).value;
x = xb2.data(16).value;
z = squeeze(xb.data(3).value(1,:,:));
figure;pcolor(x,y,z);caxis([-2 2]);shading interp; colormap jet
%plot(x,z,'k'); hold on
dxmin = 1;
dymin = 10;
%xm = xbo.data(1).value.data(13).value;
%zm = squeeze(xbo.data(9).value(1,1,:));
%file = importdata('/Users/tmiesse/work/FHRL/eeslr/modelling/xbeach/erratum/eslr_study/es_xs.csv',',');
%x = file.data(:,2);
%z = file.data(:,3);
[xx,zz] = xb_grid_xgrid(dx,dz,'dxmin',dxmin);

xb_bathy=xb_generate_bathy('x',x,'y',y,'z',z,'optimize',true,...
    'crop',false,'xgrid',{'dxmin',5},'ygrid',{'dymin',5},...
    'world_coordinates',true,'rotate',false,'finalise', {'zmin',8.5});
% xb_bathy=xb_generate_bathy('x',x,'y',y,'z',z,'optimize',true,...
%     'crop',false,'xgrid',{'dxmin',dxmin},'ygrid',{'dymin',dymin},...
%     'world_coordinates',true,'rotate',false,'finalise', {'zmin',8.5});

ygrid                   = xs_get(xb_bathy,'yfile.yfile');
xgrid                   = xs_get(xb_bathy,'xfile.xfile');
zgrid                   = xs_get(xb_bathy,'depfile.depfile');
figure;pcolor(xgrid,ygrid,zgrid);shading interp;caxis([-2 2]); colormap jet
%figure;plot(xgrid,zgrid);
%%      Fine tune the grid 
ygrid                   = xs_get(xb_bathy,'yfile.yfile');
xgrid                   = xs_get(xb_bathy,'xfile.xfile');
zgrid                   = xs_get(xb_bathy,'depfile.depfile');
figure;pcolor(xgrid,ygrid,zgrid);caxis([-2 2]);shading interp; colormap jet

% xb_bathy                         = xs_set(xb_bathy, 'yfile.yfile', y);
% xb_bathy                         = xs_set(xb_bathy, 'xfile.xfile', x);
% xb_bathy                         = xs_set(xb_bathy, 'depfile.depfile', z);
[xc yc, dim dir idx] = xb_get_coastline(xgrid, ygrid, zgrid);
figure;pcolor(xgrid,ygrid,zgrid);caxis([-2 2]);shading interp; colormap jet
figure;pcolor(zgrid);caxis([-2 2]);shading interp; colormap jet
pcolor(idx);shading interp; colormap jet
%%      Create Vegetation Map
[r,c] = size(xgrid);
xveg(1:r,1:c) = 0;
veg_h(1:r,1:c) = 0;
for i=1:r
    for ii=1:c
        if (0.9>zgrid(i,ii)) && (zgrid(i,ii)>=-0.015)
            xveg(i,ii)=1;
            veg_h(i,ii)=.4319;
%         elseif (0.25>zgrid(i,ii)) && (zgrid(i,ii)>=0.1)
%             xveg(i,ii)=2;
%             veg_h(i,ii)=.6319;   
%         elseif (zgrid(i,ii)>=0.25)
%             xveg(i,ii)=3;
%             veg_h(i,ii)=.8319;   
        else
            xveg(i,ii) = 0;
            veg_h(i,ii) = 0;
        end
    end
end

nsec  =   [1];               % number of vertical sections (only for mangrroves)
ah    =  [0.8];    % vegetation height of each section (m) - [roots -> trunk -> canopy]
bv    =  [0.025]; % stem diameter / blade width (m)
N     =  [250];      % density (units/m2)
Cd    =  [4];


xb_veg=xb_generate_settings('vegetation',1,'veggiefile','vegetation.txt','veggiemapfile','spartina_map.txt');

%figure;
%pcolor(xgrid,ygrid,xveg);shading interp; colormap jet;%caxis([-1 2])
%%      create bedfriction map

bedfrict(1:r,1:c)=0;
for i=1:r
    for ii=1:c
        if zgrid(i,ii)<0
            bedfrict(i,ii)=.02;
        end
        if zgrid(i,ii)>=0
            bedfrict(i,ii)=.045;
        end
    end
end

xb_set                     = xs_set(xb_set, 'bedfricfile', xs_set([], 'bedfricfile', bedfrict)); 

%% Create the Parameter file and any other necessary file to run xbeach
% *reminder a linux system cannot read a map file created by using the
% fprint function instead use the dlmwrite with a precision of 3.
path = '/Users/tmiesse/work/dewberry/vb_terrace/modeling/xbeach/calib2/';
xb_tot=xs_join(xb_bathy,xb_set);
xb_write_input([path, 'params.txt'], xb_tot)

fid = fopen([path,'vegetation.txt'],'w');
    fprintf(fid,'%s\n',[path,'spartina.txt']);
fclose(fid);
fid = fopen([path, 'spartina.txt'],'w');
     fprintf(fid,'%s\n', ['nsec = ', num2str(nsec)]);
     fprintf(fid,'%s\n', ['ah = ',num2str(ah)]);
     fprintf(fid,'%s\n', ['bv = ',num2str(bv)]);
     fprintf(fid,'%s\n', ['N  = ',num2str(N)]);
     fprintf(fid,'%s\n', ['Cd = ',num2str(Cd)]);
fclose(fid);

dlmwrite([path, 'spartina_map.txt'],xveg,'precision',3);

bedfric                   = xs_get(xb_tot,'bedfricfile.bedfricfile');
dlmwrite([path, 'bedfricfile.txt'],bedfric,'precision',3);


