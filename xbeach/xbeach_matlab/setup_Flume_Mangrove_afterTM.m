%% Set up for OSU Mangrove wave attenuation experiments
%
% This script generates a simple idealized 1D XBeach model(surfbeat mode)
% for a case with mangrove trees in In the Oregon State Wave Flume and is
% meant to model planned physical tests that will be on going during 
% the summer of 2020.It is a modification of a script developed for a 2-D 
% XBeach mangrove model and is based on the work by Phan et al. (2015, and
% makes use of a number of open source MATLAB functions that are part of
% the Open Earth initiative and can be download freely via
% https://publicwiki.deltares.nl/display/OET/OpenEarth
%
% For more information on the source files, please refer to:
%
%https://xbeach.readthedocs.io/en/latest/tutorials/vegetation_field2d/
% 
%  Last Modified: June, 2020 By Kiernan Kelty
%%


clear;clc;close all

%% Initialize generating auto script
% initialize constant variables (e.g. path to generate files, wave height)
pathout = '/Users/tmiesse/work/libraries/xbeach/xbeach_matlab';
Hm0   = 0.3;
Tp    = 1;
zs0   = 2.5;
dxmin = 0.1;
Lh    = 18;          % length over which vegetation is present
xveg  = [ 125 125+Lh];  % location of vegetation
nsec  = 3;             % number of vertical sections (only for mangroves)
ah    = [1.26 0.8 1.6]; % vegetation height of each section(m) [roots -> trunk -> canopy]
bv    = [0.029 0.11 0.041];   % stem diameter / blade width (m)
N     = [4.4 0.37 4.4];    % density (units/m2)
drag  = [1 1 1];         % drag coefficient (-)
% Define profile(based on 1:20 beach slope)
x     = [274 241.4481 241.4481 215.6363 25.0316 17.7661 0];
z     = [4 4 2.9665 0.8668 0.8668 0 0];
x     = fliplr(x);
z     = fliplr(z);

% ignore the naming scheme haha
% i did this for simulating xbeach 1000 times

lazy = fopen([pathout, '/copy.bat'],'wt');
lazier = fopen([pathout, '/auto_simulate.bat'],'wt');
%fprintf(lazy,'%s\n','#!bin/bash');
%fprintf(lazier,'%s\n','#!bin/bash');

%%
zs = [1,1.5,2,2.5,3];
for ii=1:2
    
    % do something like this for the values you are changing
    zs0 = zs(ii);
    
    [xx,zz] = xb_grid_xgrid(x,z,'dxmin',dxmin,'Tm',Tp/1.2,'wl',zs0);
    xb_bathy = xb_generate_bathy('x', xx, 'z', zz,...
    'crop',false,'finalise',false,'world_coordinates',false,...
    'rotate',false,'optimize',false);
     
    str = ([pathout,'grid',num2str(length(xx)),'_zs',num2str(zs0*10)]);
    %str2= (['test','/','drag',num2str(Cd2),'_waves',num2str(wave_h2)]);
    mkdir(str);
    
    xb_wave  = xb_generate_waves('Hm0',Hm0,'Tp',Tp,'duration',1200,'mainang',270,...
    's', 100000,'fnyq',2,'gammajsp',7);
    
    xb_set = xb_generate_settings('sedtrans',0,'morphology',0,'roller',1,'CFL',0.9,...
    'order',1,'ARC',0,'rt',600,'cf',0.001,'zs0',zs0,...
    'front','waveflume','left','wall','right','wall',...
    'back','abs_2d','posdwn',-1,'depthscale',50,'rho',1000,'random',0,...
    'thetamin',-90,'thetamax',90,'dtheta',180,...
    'tstart',600,'tstop',1200,'tintg',0.1,'tintm',600,...
    'outputformat','netcdf',...
    'nglobalvar',{'H','zs','u','zb','k','sigm'},...
    'nmeanvar',{'H','zs','k','sigm'});

    veg = zeros(size(xx));
    veg(xx >= xveg(1) & xx < xveg(2)) = 1;

    % Save vegetaion in put in XBeach structure
    xb_veg = xb_generate_settings('vegetation',1,...
    'veggiefile','vegetation.txt','veggiemapfile','mangrovebed.txt');
    
    xb_tot=xs_join(xb_bathy,xb_veg,xb_set,xb_wave);
    xb_write_input([str, '\params.txt'], xb_tot)

    fid = fopen([str,'/vegetation.txt'],'w');
    fprintf(fid,'%s\n','spartina.txt');
    fclose(fid);

    fid = fopen([str, '/spartina.txt'],'w');
        fprintf(fid,'%s\n', ['nsec = ', num2str(nsec)]);
        fprintf(fid,'%s\n', ['ah = ',num2str(ah)]);
        fprintf(fid,'%s\n', ['bv = ',num2str(bv)]);
        fprintf(fid,'%s\n', ['N  = ',num2str(N)]);
        fprintf(fid,'%s\n', ['Cd = ',num2str(drag)]);
    fclose(fid);
    
    

end
% Save wave input in XBeach structure

%% 2) Mean water level
zs0  = 2.5; % mean water level [m]

%% 3) Grid and Bathymetry profile 
dxmin = 0.1;% minimum grid size (near shore line)


% Use xb_grid_xgrid tool to get optimized cross-shore grid


% Save grid/bathy in XBeach structure



%% 4) Vegetation input


% Find grid points with vegetation


%% 5) General Model settings
% Directly save model settings in XBeach structure


%% 6) Plot experimental setup
figure;

% Plot profile
figure(1)
plot(xx,zz,'k');hold on
plot([274 0],[zs0 zs0],'b--');hold on
plot([274 0],[2.5 2.5],'b-.');hold on
ylim([0 5])
xlim([0 270])
% Plot vegetation
x_mang = [ 125 143 143 125 125];
z_mang_roots = [0.8668  0.8668 2.1268 2.1268 0.8668];
z_mang_trunk = [2.1268  2.1268 2.9268 2.9268 2.1268];
z_mang_canopy = [2.9268  2.9268 3.9148 3.9148 2.9268];

plot(x_mang,z_mang_roots,'r')
hold on 
plot(x_mang,z_mang_trunk,'b')
hold on 
plot(x_mang,z_mang_canopy,'g')

    
% Figure labels
xlabel('x [m]')
ylabel('z [m]')
legend('Cross Shore Profile','Water Level, h = 1.77m','Water Level, h = 2.5m','Mangrove roots','Mangrove Trunk','Mangrove Canopy','location','Best')

% Title
title({'Experimental setup: Mangrove Wave','Attenuation Experiments, OSU Wave Flume'})





% %% 7) Save model setup
% 
% % Define location / path of XBeach executable
% binDir = [pwd filesep 'XBeachExecutable'];
% 
% % 7.1) Save model setup for test 1 (no mangroves)
% dirName1 = 'testb4.01_kelty';mkdir(dirName1)
% 
% % Combine XBeach structures
% xb_tot = xs_join(xb_bat,xb_wav,xb_set);
% xb_write_input([dirName1 '\params.txt'],xb_tot);
% 
% 
% % Define location / path of XBeach executable
% binDir = [pwd filesep 'XBeachExecutable'];
% 
% % 7.2) Save model setup for test 2 (with mangroves)
% dirName2 = 'testm4.01_kelty';mkdir(dirName2)
% xb_tot = xs_join(xb_bat,xb_wav,xb_veg,xb_set);
% xb_write_input([dirName2 '\params.txt'],xb_tot);
% 
% 
% % Save additional vegetation input files (3x) Veggiefile (lists all
% % vegetation types in model):
% fid = fopen([dirName2 '\vegetation.txt'],'w');
%     fprintf(fid,'%s','mangroves.txt');
% fclose(fid);
% 
% % Veggiechars (lists all parameters for specified vegetation type):
% fid = fopen([dirName2,'\mangroves.txt'],'w');
%     fprintf(fid,'%s\n', ['nsec = ', num2str(nsec)]);
%     fprintf(fid,'%s\n', ['ah = ',num2str(ah)]);
%     fprintf(fid,'%s\n', ['bv = ',num2str(bv)]);
%     fprintf(fid,'%s\n', ['N  = ',num2str(Nv)]);
%     fprintf(fid,'%s\n', ['Cd = ',num2str(Cd)]);
% fclose(fid);
% 
% % Veggiemapfile (location of vegetation types):
% save([dirName2 '\mangrovebed.txt'],'veg','-ascii');% vegetation characteristitcs
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
