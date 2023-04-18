%%

run '/Users/tmiesse/work/libraries/delft/oetsettings.m';
clear all; clc;


%% Load model results

root = '/Users/tmiesse/Downloads';


xbo=xb_read_output(root);
xbo2 = xbo.data(1).value;

%%

t  = xbo2.data(24).value;
H  = xbo.data(17).value;
zs = xbo.data(13).value;
zb = xbo.data(9).value;
x  = xbo2.data(23).value;
y  = xbo2.data(24).value;
any = max(H);
%% Only for 1D model 
baseline=-8;
	index=1:499;
    brown= [0.8 0.8 .6];
    blue2=[0 0.3 1];           
    zs1 = zs(5,:,:);
    zb1 = zb(1,:,:);
    A=figure();
        waves = zs(10,1,:) + H(10,1,:);
        A(1)=plot(x,waves,'c'); hold on
        A(2)=plot(x,zs(10,1,:),'b'); hold on
        A(3)=plot(x,zb(10,1,:),'k'); hold on

        A(6)=fill(x(index([ 1 1:end end])),...
               [baseline waves(index) baseline],...
                blue2,'EdgeColor','none'); alpha(.5)                     
        A(7)=fill(x(index([ 1 1:end end])),...
               [baseline zs1(index) baseline],...
               'b','EdgeColor','none'); alpha(.5)                
        A(8)=fill(x(index([1 1:end end])),...
               [baseline zb1(index) baseline],...
               brown,'EdgeColor','none');



%% For 2D model


figure;

c(2) = pcolor(x,y,squeeze(zb(2,:,:))); shading interp;hold off


colormap('jet'); colorbar




%%










%%


%%



%%