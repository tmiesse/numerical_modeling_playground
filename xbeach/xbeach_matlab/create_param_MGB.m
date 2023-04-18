%% Code to make Parameter file for xbeach in 2d with single vegetation
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

destin          = 'Z:\Project_NFWF\6_WORKSPACE\Tyler\xbeach\xbeach_example\magothy_bay_attempt\';
destin2          = 'Z:\Project_NFWF\6_WORKSPACE\Tyler\xbeach\xbeach_example\magothy_bay_attempt\MGB_files_to_create_bathy\';
destinbathy     = 'Z:\Project_NFWF\6_WORKSPACE\Tyler\xbeach\xbeach_example\magothy_bay_attempt\MGB_bathy.mat_files\';
destout         = 'Z:\Project_NFWF\6_WORKSPACE\Tyler\xbeach\xbeach_example\magothy_bay_attempt\mgb_grid_final_input_NO_OUTPUT_FILES_ALLOWED\working_on_the_tide\';
mkdir(destout);
desttide          = 'Z:\Project_NFWF\6_WORKSPACE\Tyler\xbeach\xbeach_example\magothy_bay_attempt\Bday_tide_and_wind\';
destswan        = 'Z:\Project_NFWF\6_WORKSPACE\Tyler\xbeach\xbeach_example\magothy_bay_attempt\juan_swan_outputs\';

% determine the grid size of the topo/bathy
dxmin               = 1;
dymin               = 1;
outputformat        = 'fortran';
        
%%      Call the topo/bathy file 

cd(destinbathy)
load 'mgb_bathy_standard.mat'
x=fliplr(output_X2);
y=output_Y2;
z=fliplr(output_Z_CRM2);
figure;pcolor(x,y,z); shading flat
%%      Determine swan files if being used
% *reminder (8/29) currently coupling swan with xbeach is to much for stampede
%  be aware that it uses to much resources and it will give the error that
%  there is not enough memory to store the output

% cd(destswan)
% load 'swan_inc_output.mat';
% s=S(31,:).';
% swan=xb_swan_read([destswan,'stations.spc']);
% swan2=xb_swan_split(swan,'location');
% waves=xb_swan_write([destout,'wave1'],swan2(1));
% waves=xb_swan_write([destout,'wave2'],swan2(2));
% waves=xb_swan_write([destout,'wave3'],swan2(3));
% waves=xb_swan_write([destout,'wave4'],swan2(4));
% waves=xb_swan_write([destout,'wave5'],swan2(5));
% waves=xb_swan_write([destout,'wave6'],swan2(6));
% waves=xb_swan_write([destout,'wave7'],swan2(7));

%%      determine the tide either use data files or xbeach default

tide=xlsread([destin2,'hurricane_matthew_try']);

tide_height(1:220,:)=tide(1:220,6);

tide_time(1:220,1)=0;
for i=1:219
tide_time(1+i,1)=tide_time(i,1)+360; %3 hours of tide in seconds
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
    'world_coordinates',true,'rotate',false,'finalise', {'zmin',1.7});

%%      Generate the waves 

xb_wave=xb_generate_waves('Hm0', 1, ... 
        'Tp', 4,'mainang', 90, ...
        'gammajsp',3.3,'s',10,'fnyq',.5);
    
%%      Generate the tide

  xb_tide=xb_generate_tide('time',tide_time,'front',tide_height);    

%%      Generate the settings for xbeach

xb_set=xb_generate_settings('outputformat',outputformat,... 
        'thetamin',45,'thetamax',135,'dtheta',90,'dtheta_s',10,...
        'instat','jons','morfac', 10,'posdwn',-1,...
        'morstart', 0,'CFL', 0.7,'front', 'abs_2d', ...
        'back', 'abs_2d','left','neumann','right','neumann','mpiboundary','auto',...
        'thetanaut', 1,'zs0',1,...
        'tstop', 36001,'tstart', 0,'bcfile','filelist.txt',...
        'tint', 36000,'tintm',600,'tintg',600,'epsi',-1,'facua',0.10,'bedfriction', 'manning',...
        'meanvar',{'zb', 'zs', 'H','u','v','sedero'} ,...
        'globalvar',{'zb', 'zs','H','u','v','sedero'});
    
%%      Call the grid files from xb_bathy
%*reminder the x and z grid need to be flipped for xbeach to properly read
%the files. This make sure the grid is in the proper direction.
  
xgrid                   = xs_get(xb_bathy,'xfile.xfile');
ygrid                   = xs_get(xb_bathy,'yfile.yfile');
zgrid                   = xs_get(xb_bathy,'depfile.depfile');
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
for i=822:827
    zgrid(1:182,i)=1.75;
end
xb_bathy                         = xs_set(xb_bathy, 'depfile.depfile', zgrid);

%%      Create Vegetation Map
% currently the vege map is being created by using the elvation. This works
% for single and double vegetation but it can not be used for more than
% that. In the future to keep up with accuracy come up with new code for
% more than 2 vege and it goes by area.

xveg=zeros(size(zgrid));

for i=1:182
    for ii=1:827
        if zgrid(i,ii)>.1
            xveg(i,ii)=1;
        else
            xveg(i,ii)=0;
        end
    end
end

npts = 2;               % number of vertical sections (only for mangrroves) * it mmust be 2 or greater
zv   = [.005 .75];      % vegetation height of each section (m) - [roots -> trunk -> canopy]
bv   = [.03 .02];       % stem diameter / blade width (m)
N    = [200 200];        % density (units/m2)
Cd   = [1 1];   

xb_veg=xb_generate_settings('vegetation',1,'veggiefile','vegetation.txt','veggiemapfile','spartina_map.txt');

%%      create bedfriction map

for i=1:182
    for ii=1:827
        if xveg(i,ii)==1
            bedfrict(i,ii)=.012;
        else
            bedfrict(i,ii)=0.02;
        end
    end
end

xb_set                     = xs_set(xb_set, 'bedfricfile', xs_set([], 'bedfricfile', bedfrict)); 


%%      Generate the wind settings

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


%% Plot initial grid from xbeach and the topo/bathy from gis

figure;pcolor(xgrid,ygrid,zgrid); shading flat
figure;pcolor(x,y,z); shading flat

%% Create the Parameter file and any other necessary file to run xbeach
% *reminder a linux system cannot read a map file created by using the
% fprint function instead use the dlmwrite with a precision of 3.

cd(destout);
xb_tot=xs_join(xb_bathy,xb_wave,xb_veg,xb_set,xb_tide);
%xb_write_input([destin '\params.txt'],xb_tot);

xb_write_input([destout, '\params.txt'], xb_tot)
fid = fopen([destout,'\filelist.txt'], 'w');
    fprintf(fid, '%s\n','wave1');
    fprintf(fid, '%s\n','wave2');
    fprintf(fid, '%s\n','wave3');
    fprintf(fid, '%s\n','wave4');
    fprintf(fid, '%s\n','wave5');
    fprintf(fid, '%s\n','wave6');
    fprintf(fid, '%s\n','wave7');
fclose(fid);

fid = fopen([destout,'\vegetation.txt'],'w');
    fprintf(fid,'%s\n','spartina.txt');
fclose(fid);

fid = fopen([destout, '\spartina.txt'],'w');
     fprintf(fid,'%s\n', ['npts = ', num2str(npts)]);
     fprintf(fid,'%s\n', ['zv = ',num2str(zv)]);
     fprintf(fid,'%s\n', ['bv = ',num2str(bv)]);
     fprintf(fid,'%s\n', ['N  = ',num2str(N)]);
     fprintf(fid,'%s\n', ['Cd = ',num2str(Cd)]);
fclose(fid);

dlmwrite([destout, '\spartina_map.txt'],xveg,'precision',3);

bedfric                   = xs_get(xb_tot,'bedfricfile.bedfricfile');
save('bedfricfile.txt', 'bedfric', '-ascii')

