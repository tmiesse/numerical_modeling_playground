%% XBeach JIP tutorial: how to set-up a storm-impact model? 
% v0.1 Nov-15
clear all, close all, clc


%% 0. Define
destin = '/Users/tmiesse/work/FHRL/dunex/modeling/gis/dem_con_ascii/';

% Determine area of interest (not boundaries from XBeach grid, but close)
output_lon_w        = -75.658619;
output_lon_e        = -75.653557;
output_lat_s        = 35.221179;
output_lat_n        = 35.225334; 
output_delta_x      = 1.004953200000000013e-05;
output_delta_y      = 1.004953200000000013e-05;
% Determine size of grid cells
% output_delta_x      = 1;
% output_delta_y      = 1;
[output_X,output_Y] = meshgrid(output_lon_w:output_delta_x:output_lon_e, output_lat_n:-output_delta_y:output_lat_s);

%% 1. Load CRM

[input_Z_CRM, R_CRM]            = arcgridread_v2([destin,'usgs_lidar.asc']); % 

% Get grid information from input dataset     
[input_n_y, input_n_x]           = size(input_Z_CRM);
input_delta_x                   = R_CRM(2,1);
input_delta_y                   = R_CRM(1,2);
input_lon_w                     = R_CRM(3,1);
input_lon_e                     = R_CRM(3,1)+input_delta_x*(input_n_x-1);
input_lat_n                     = R_CRM(3,2);
input_lat_s                     = R_CRM(3,2)+input_delta_y*(input_n_y-1);
[input_X_CRM, input_Y_CRM]          = meshgrid(input_lon_w:R_CRM(2,1):input_lon_e,input_lat_n:R_CRM(1,2):input_lat_s);
% Griddata to output grid
output_Z_CRM                     = griddata(input_X_CRM, input_Y_CRM,input_Z_CRM, output_X, output_Y);


%%
[input_Z_CRM2, R_CRM2]            = arcgridread_v2([destin, 'VBDEM10ft_newsurvey.asc']); % 
 
% Get grid information from input dataset     
[input_n_y2, input_n_x2]           = size(input_Z_CRM2);
input_delta_x                   = R_CRM2(2,1);
input_delta_y                   = R_CRM2(1,2);
input_lon_w2                     = R_CRM2(3,1);
input_lon_e2                     = R_CRM2(3,1)+input_delta_x*(input_n_x2-1);
input_lat_n2                     = R_CRM2(3,2);
input_lat_s2                     = R_CRM2(3,2)+input_delta_y*(input_n_y2-1);
[input_X_CRM2, input_Y_CRM2]          = meshgrid(input_lon_w2:R_CRM2(2,1):input_lon_e2,input_lat_n2:R_CRM2(1,2):input_lat_s2);
% Griddata to output grid

output_Z_CRM2                     = griddata(input_X_CRM2, input_Y_CRM2,input_Z_CRM2, output_X, output_Y);
%%
idn = find(output_Z_CRM<=-3.40282306073710e+10);
output_Z_CRM(idn) = NaN;

idcrm = (isnan(output_Z_CRM) | isnan(output_Z_CRM2));
output_Z_CRM(idcrm) = output_Z_CRM2(idcrm);



%%
%figure; pcolor(input_X_CRM ,input_Y_CRM , input_Z_CRM); shading interp;  

figure;pcolor(output_X,output_Y,output_Z_CRM);caxis([-2 2]);shading interp; colormap jet;colorbar






%%
save(['/Users/tmiesse/work/FHRL/dunex/modeling/xbeach/frisco_field_site/dem_mat/nc_site.mat'],('output_X'),('output_Y'),('output_Z_CRM'));

%%
