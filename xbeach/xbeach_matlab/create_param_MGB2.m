%% Code to make Parameter file for xbeach in 2d with multiple vegetation
% Tyler Miesse

% Currently (8/29) this code is meant to generate magothy bay from the gps
% data. 
%   The settings include:
%       -Tide which uses data from noaa.
%       -waves using the jonswap setting in xbeach
%       -vegetation
%       -maybe wind that comes from noaa

clear all; close all; clc;

%%      Define Path Directory
% This helps with keeping everything in order instead of writing the file
% path multiple times

destin          = 'Z:\Project_NFWF\3_Modeling\2_xBeach\magothy_bay_attempt\';
destin2         = 'Z:\Project_NFWF\3_Modeling\2_xBeach\magothy_bay_attempt\MGB_files_to_create_bathy\';
destinbathy     = 'Z:\Project_NFWF\3_Modeling\2_xBeach\magothy_bay_attempt\MGB_bathy.mat_files\';
destout         = 'Z:\Project_NFWF\3_Modeling\2_xBeach\magothy_bay_attempt\MGB_sensitivity_test_inputs\2d_test\';
mkdir(destout);
%desttide        = 'Z:\Project_NFWF\19_WORKSPACE\Tyler\xbeach\xbeach_example\magothy_bay_attempt\Bday_tide_and_wind\';
destswan        = 'Z:\Project_NFWF\3_Modeling\2_xBeach\magothy_bay_attempt\juan_swan_outputs\';
% B. Model
dxmin           = 5;
dymin           = 5;
outputformat    = 'fortran';
        
%%      Call the topo/bathy file 

cd(destinbathy)
load 'mgb_bathy_standard.mat'
x=fliplr(output_X2);
y=output_Y2;
z=fliplr(output_Z_CRM2);
%figure;pcolor(x,y,z); shading flat
%%      Determine swan files if being used
% *reminder (8/29) currently coupling swan with xbeach is to much for stampede
%  be aware that it uses to much resources and it will give the error that
%  there is not enough memory to store the output

% cd(destswan)
% load 'swan_inc_output.mat';
% s=S(31,:).';
% swan=xb_swan_read([destswan,'stations.spc']);
% swan2=xb_swan_split(swan,'location');
% waves=xb_swan_write([destout,'\wave.txt'],swan2);
% xb_swan=xb_generate_settings('instat','swan','bcfile','filelist.txt','rt',36000,'dtbc',360);

%%      determine the tide either use data files or xbeach default

tide=importdata([destin2,'WLFort63.txt']);

tide_h(:,1)=tide.textdata(4:3:end,2);
tide_t(:,1)=tide.textdata(3:3:end,2);

tide_height=str2double(tide_h(:,1));
tide_time=str2double(tide_t(:,1));
tide_time(1:end,1)=3600;
for i=1:length(tide_time(1:(end-1),1))
    tide_time(i+1,1)=tide_time(i,1)+3600;
end

%%      read wind file if being used in model
% *reminder using wind is similar to using tide but it needs to be in a
%  different format and xbeach does not have its own function to read wind.
%  This means that generate settings functions needs to be used.

% wind=xlsread([desttide,'wind_20170620.xlsx']);
% wind_speed=wind(1:220,6);
% wind_dir=wind(1:220,7);
% wind_time(1:220,1)=0;
% for i=1:219
% wind_time(1+i,1)=wind_time(i,1)+360; %3 hours of tide in seconds
% end

%%      Creating grid from the bathymetry defined earlier

xb_bathy=xb_generate_bathy('x',x,'y',y,'z',z,...
    'xgrid',{'dxmin',dxmin},'ygrid',{'dymin',dymin},'crop',false,...
    'world_coordinates',true,'rotate',false,'finalise', {'zmin',1.6});
% % create the xbeach bathy again or us the saved generated bathy. It should
% % be preferred to load the save file instead of re-generating the bathy
% % unless there was a change.

%load(['Z:\Project_NFWF\3_Modeling\2_xBeach\magothy_bay_attempt\MGB_generated_bathy\','generated_bathy_for_MGB']);
%load(['Z:\Project_NFWF\3_Modeling\2_xBeach\magothy_bay_attempt\MGB_generated_bathy\','generated_bathy_for_MGB2']);
%%





%%      Generate the waves 

xb_wave=xb_generate_waves('Hm0', 0.34525, ... 
        'mainang', 100,'Tp',4,...
        'gammajsp',3.3,'s',10,'fnyq',0.45);
    
%%      Generate the tide

  xb_tide=xb_generate_tide('time',tide_time,'front',tide_height);    

%%      Generate settings for xbeach to follow

xb_set=xb_generate_settings('outputformat',outputformat,... 
        'thetamin',0,'thetamax',180,'dtheta',90,'dtheta_s',10,...
        'morfac', 1,'posdwn',-1,'instat','jons',...
        'morstart', 0,'CFL', 0.75,'front', 'abs_2d','dtbc',2, ...
        'back', 'abs_2d','left','neumann','right','neumann','mpiboundary','auto',...
        'thetanaut', 1,'zs0',.45,'random',0,'taper',1,...
        'tstop', 86401,'tstart', 0,...
        'tint', 86400,'tintm',1200,'tintg',36000,'epsi',-1,'facua',0.30,'bedfriction', 'manning',...
        'meanvar',{'zb', 'zs', 'H','urms','Cdrag'} ,...
        'globalvar',{'zb', 'zs','H','urms','Cdrag'});
    
%%      Call the grid files from xb_bathy
%*reminder the x and z grid need to be flipped for xbeach to properly read
%the files. This make sure the grid is in the proper direction.

xgrid                   = xs_get(xb_bathy,'xfile.xfile');
ygrid                   = xs_get(xb_bathy,'yfile.yfile');
zgrid                   = xs_get(xb_bathy,'depfile.depfile');
%%
ygrid = ygrid(5:91,1:510);
xgrid = xgrid(5:91,1:510);
zgrid = zgrid(5:91,1:510);
%%
xgrid=fliplr(xgrid);
zgrid=fliplr(zgrid);
xb_bathy                         = xs_set(xb_bathy, 'xfile.xfile', xgrid);
xb_bathy                         = xs_set(xb_bathy, 'depfile.depfile', zgrid);

%%      Clean errors from the elevation file
% set the max and min that the elevation can be

id1 = find(zgrid > 1.75);
for i = 1:length(id1)
    zgrid(id1(i)) = 1.75;
end
% 
% id2 = find(zgrid < -5);
% for i = 1:length(id2)
%     zgrid(id2(i)) = -5;
% end
% %%
% C. Straight boundaries
[nx, ny]                            = size(zgrid);
roundnumber                         = 1;
first                               = 1;
second                              = roundnumber;
third                               = nx-(roundnumber-1);
four                             	= nx;
five                                = ny - (roundnumber-1);
six                                 = ny;
zgrid([first:second],:)             = repmat(zgrid(second,:),[roundnumber,1]); 
zgrid([(third:four)],:)             = repmat(zgrid((third),:),[roundnumber,1]);
zgrid(:,[five:six])                 = repmat(zgrid(:,six),[1,roundnumber]);
% for i=822:827
%     zgrid(1:182,i)=1.75;
% end
xb_bathy                         = xs_set(xb_bathy, 'depfile.depfile', zgrid);
%%
figure;
pcolor(xgrid,ygrid,zgrid); shading flat
%%      Create Vegetation Map
% currently the vege map is being created by using the elvation. This works
% for single and double vegetation but it can not be used for more than
% that. In the future to keep up with accuracy come up with new code for
% more than 2 vege and it goes by area.

xveg=zeros(size(zgrid));

for i=1:56
    for ii=1:376
        if zgrid(i,ii)>.25
            xveg(i,ii)=1;
            veg_h(i,ii)=.25;
        end
        if zgrid(i,ii)>0 && zgrid(i,ii)<.25
            xveg(i,ii)=2;
            veg_h(i,ii)=.65;
        end
        if zgrid(i,ii)<=0
            xveg(i,ii)=0;
            veg_h(i,ii)=0;
        end
    end
end
nsec  = 1;               % number of vertical sections (only for mangrroves)
ah    = .25;    % vegetation height of each section (m) - [roots -> trunk -> canopy]
bv    = .02; % stem diameter / blade width (m)
N     = 200;      % density (units/m2)
Cd    = -3;
nsec2 = 1;               % number of vertical sections (only for mangrroves)
ah2   = 1.25 ;    % vegetation height of each section (m) - [roots -> trunk -> canopy]
bv2   = .03; % stem diameter / blade width (m)
N2    =  200;      % density (units/m2)
Cd2   = -3; 

xb_veg=xb_generate_settings('vegetation',1,'veggiefile','vegetation.txt','veggiemapfile','multi_spartina_map.txt');

%% create bedfriction map

for i=1:56
    for ii=1:376
        if xveg(i,ii)==1
            bedfrict(i,ii)=.02;
        else
            bedfrict(i,ii)=0.0225;
        end
    end
end

xb_set                     = xs_set(xb_set, 'bedfricfile', xs_set([], 'bedfricfile', bedfrict)); 


%%      Generate wind settings

% xb_wind=xb_generate_settings('wind',1,'windfile','wind.txt');
% 
% fid = fopen([destout,'\wind.txt'],'w');
%     fprintf(fid,'%s\n','wind_details.txt');
% fclose(fid);
% 
% wind_time=wind_time.';
% wind_speed=wind_speed.';
% wind_dir=wind_dir.';
% 
% fid = fopen([destout, '\wind_details.txt'],'w');
%     %fprintf(fid,'%s\n', ['time = ', num2str(wind_time)]);
%     fprintf(fid,'%s\n', ['windth = ', num2str(wind_dir)]);
%     fprintf(fid,'%s\n', ['windv = ', num2str(wind_speed)]);


%%      Plot initial grid from xbeach and the topo/bathy from gis

figure;
pcolor(xgrid,ygrid,zgrid);hold on
pcolor(xgrid,ygrid,veg_h+zgrid);
% for j=1:182
%     for jj = 1:827
%     pcolor([xgrid(j,jj),xgrid(j,jj)],[zgrid(j,jj) zgrid(j,jj)+veg_h(j,jj)],'Color',[76 153 0]/255);
%     end
% end
shading flat
colorbar
title('Bathymetry of Magothy Bay');
% [x1,y1]=ginput(1);
% close;
% [num, idx]=min(abs(xgrid(1,:)-x1));
% [num2, idy]=min(abs(ygrid(:,2)-y1)); 

%%      Create the Parameter file and any other necessary file to run xbeach
% *reminder a linux system cannot read a map file created by using the
% fprint function instead use the dlmwrite with a precision of 3.

cd(destout);
xb_tot=xs_join(xb_bathy,xb_veg,xb_set,xb_tide,xb_wave);
%xb_write_input([destin '\params.txt'],xb_tot);

xb_write_input([destout, '\params.txt'], xb_tot)
% 
% fid = fopen([destout,'\filelist.txt'], 'w');
%     fprintf(fid,'%s\n', 'FILELIST');
%         fprintf(fid, '%10i%10.4f%50s\n', 36000, 36000, ('wave001.txt'));      
% fclose(fid);

fid = fopen([destout,'\vegetation.txt'],'w');
    fprintf(fid,'%s\n','spartina.txt');
    fprintf(fid,'%s\n','spartina2.txt');
fclose(fid);

fid = fopen([destout, '\spartina.txt'],'w');
     fprintf(fid,'%s\n', ['nsec = ', num2str(nsec)]);
     fprintf(fid,'%s\n', ['ah = ',num2str(ah)]);
     fprintf(fid,'%s\n', ['bv = ',num2str(bv)]);
     fprintf(fid,'%s\n', ['N  = ',num2str(N)]);
     fprintf(fid,'%s\n', ['Cd = ',num2str(Cd)]);
fclose(fid);
fid = fopen([destout, '\spartina2.txt'],'w');
     fprintf(fid,'%s\n', ['nsec = ', num2str(nsec2)]);
     fprintf(fid,'%s\n', ['ah = ',num2str(ah2)]);
     fprintf(fid,'%s\n', ['bv = ',num2str(bv2)]);
     fprintf(fid,'%s\n', ['N  = ',num2str(N2)]);
     fprintf(fid,'%s\n', ['Cd = ',num2str(Cd2)]);
fclose(fid);
dlmwrite([destout, '\multi_spartina_map.txt'],xveg,'precision',3);

bedfric                   = xs_get(xb_tot,'bedfricfile.bedfricfile');
dlmwrite([destout, '\bedfricfile.txt'],bedfric,'precision',6);

