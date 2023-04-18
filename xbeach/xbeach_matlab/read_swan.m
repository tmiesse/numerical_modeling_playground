%%      



run '/Users/tmiesse/work/libraries/delft/oetsettings.m';
clear all; clc;

%%      Define Path Directory
% This helps with keeping everything in order instead of writing the file
% path multiple times

path = '/Users/tmiesse/work/FHRL/arctic/model/no_ice_storm2011/';

file = swan_io_spectrum([path,'2DSpecOut.nc']);


%%      
swan_io_spectrum2nc(file,[path,'2Dspecout.nc']);
%file.
%pcolor_spectral(file.frequency,file.directions)%,file.VaDens)

%%