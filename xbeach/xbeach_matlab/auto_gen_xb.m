% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
clear;  clc;
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% Define paths 
% this can all be done on one or two paths
run '/Users/tmiesse/work/libraries/delft/oetsettings.m';

destin          = 'Z:\Project_NFWF\3_Modeling\2_xBeach\magothy_bay_attempt\';
destinbathy     = 'Z:\Project_NFWF\3_Modeling\2_xBeach\magothy_bay_attempt\MGB_bathy.mat_files\';
dest_tide       = 'Z:\Project_NFWF\3_Modeling\2_xBeach\real_time\tide_from_arslaan\';
destout         = 'Z:\Project_NFWF\3_Modeling\2_xBeach\magothy_bay_attempt\validation_inputs\Jose_v2\';
destout2         = 'Z:\Project_NFWF\3_Modeling\2_xBeach\magothy_bay_attempt\validation_inputs\';
dxmin           = 2.5;
dymin           = 2.5;
outputformat    = 'netcdf';

%% generate cross-section inputs
% this part is most likely what needs to be changed
path = [destinbathy,'mgb_bathy_standard_real.mat'];
load path;
x=(output_X2);
z=(output_Z_CRM2);
y=output_Y2;

[xgrid,zgrid]=xb_grid_xgrid(x,z,'dxmin',dxmin);

xb_bathy=xb_generate_bathy('x',xgrid,'z',zgrid,...
    'xgrid',{'dxmin',dxmin},'crop',false,'optimize',false,...
    'world_coordinates',true,'rotate',false); 

xgrid2                   = xs_get(xb_bathy,'xfile.xfile');
zgrid2                   = xs_get(xb_bathy,'depfile.depfile');
xb_bathy                         = xs_set(xb_bathy, 'xfile.xfile', xgrid2);
xb_bathy                         = xs_set(xb_bathy, 'depfile.depfile', zgrid2);

%% generate bedfriction and vegetation inputs
xveg=zeros(size(zgrid));
bedfrict=zeros(size(zgrid));
r,c = size(zgrid);
% vege types here are based on elevation
for i=1:r
    for ii=1:c
        if ii >= 315
            if zgrid2(i,ii)>.05
                xveg(i,ii)=1;
                veg_h(i,ii)=.25;
            end
            if zgrid(i,ii)>0 && zgrid(i,ii)<.25
                xveg(i,ii)=1;
                veg_h(i,ii)=.5;
            end
            if zgrid2(i,ii)<=0
                xveg(i,ii)=0;
                veg_h(i,ii)=0;
            end
            if zgrid2(i,ii)>1.25
                xveg(i,ii)=0;
                veg_h(i,ii)=.1;
            end
        else
            xveg(i,ii) = 0;
        end
    end
end

for i=1:r
    for ii=1:c
        if zgrid2(i,ii)<0
            bedfrict(i,ii)=.02;
        end
        if zgrid2(i,ii)>=0
            bedfrict(i,ii)=0.021;
        end
    end
end


%% Initialize generating auto script
cd(destout)
wave_h=0;
tide  = 1.0043;
N = 366.97;
Cd = 1;
ah = 0.476;
bv    = 0.005587;
% ignore the naming scheme haha
% i did this for simulating xbeach 1000 times

lazy = fopen([destout2, '\copy.bat'],'wt');
lazier = fopen([destout2, '\auto_simulate.bat'],'wt');
%fprintf(lazy,'%s\n','#!bin/bash');
%fprintf(lazier,'%s\n','#!bin/bash');
%% 
%for i=1:4
% for this automation I was changing the drag on vegetation

    for ii=1:6
        %tide=tide+.25;
        wave_h=0.4290;
    for iii=1:2

        wave_h2= wave_h*100;
        Cd2 = Cd*-1;
        str = ([destout,'drag',num2str(Cd2),'_waves',num2str(wave_h2)]);
        str2= (['test','/','drag',num2str(Cd2),'_waves',num2str(wave_h2)]);
        mkdir(str);
%     cd(str)
        xb_wave=xb_generate_waves('Hm0',wave_h,'gammajsp',3.3,'s',10,'Tp',2.828,...
                            'mainang',90,'fnyq',.45);   

        xb_set=xb_generate_settings('outputformat',outputformat,... 
        'thetamin',0,'thetamax',180,'dtheta',180,'dtheta_s',5,'xori',0,...
        'instat','jons','taper',1,'morfac', 1,'posdwn',-1,'avalanching',0,...
        'morstart', 0,'CFL', 0.2,'front', 'abs_1d','random',0,...
        'back', 'abs_1d','left','neumann','right','neumann','mpiboundary','auto',...
        'thetanaut', 1,'zs0',tide,'single_dir',0,...
        'tstop', 3601,'tstart', 0,'alpha',0.005,...
        'tintm',1200,'tintg',1200,'epsi',-1,'facua',0.300,'bedfriction', 'manning',...
        'meanvar',{'zb', 'zs', 'H','urms','Cdrag','Qb','E','sigm','Dveg','Df'} ,...
        'globalvar',{'zb', 'zs','H','urms','Cdrag'});

        xb_set                     = xs_set(xb_set, 'bedfricfile', xs_set([], 'bedfricfile', bedfrict)); 

        nsec  =  1;               % number of vertical sections (only for mangrroves)
                                  % vegetation height of each section (m) - [roots -> trunk -> canopy]
                                  % stem diameter / blade width (m)
                                  % density (units/m2)
        drag  = Cd;
        xb_veg=xb_generate_settings('vegetation',1,'veggiefile','vegetation.txt','veggiemapfile','spartina_tran.txt');


        xb_tot=xs_join(xb_bathy,xb_veg,xb_set,xb_wave);
        xb_write_input([str, '\params.txt'], xb_tot)

        fid = fopen([str,'\vegetation.txt'],'w');
        fprintf(fid,'%s\n','spartina.txt');
        fclose(fid);

        fid = fopen([str, '\spartina.txt'],'w');
            fprintf(fid,'%s\n', ['nsec = ', num2str(nsec)]);
            fprintf(fid,'%s\n', ['ah = ',num2str(ah)]);
            fprintf(fid,'%s\n', ['bv = ',num2str(bv)]);
            fprintf(fid,'%s\n', ['N  = ',num2str(N)]);
            fprintf(fid,'%s\n', ['Cd = ',num2str(drag)]);
        fclose(fid);

            fprintf(lazy,'%10s\n',['cp /home/vse/Neptune_work_folder/XBeach/validation/Jose xbeach ',
                '/home/vse/Neptune_work_folder/XBeach/validation/Jose/',str2]);
            fprintf(lazier,'%10s\n',['cd /home/vse/Neptune_work_folder/XBeach/validation/Jose/',str2]);
            fprintf(lazier,'%10s\n','./xbeach');
    
            dlmwrite([str, '\spartina_tran.txt'],xveg,'precision',3);

            bedfric                   = xs_get(xb_tot,'bedfricfile.bedfricfile');
            dlmwrite([str, '\bedfricfile.txt'],bedfric,'precision',6);
            if wave_h >= 2
                break
            end
        wave_h = wave_h+0.001;
   end
    Cd = Cd + 1;
    end
    %bv = bv + (.1/3);
%end
%%
fclose(lazy);
fclose(lazier);



