
clear;
close all;

% run(['/home/admin/neptune_work/realtime/xbeach/matlab_for_rt/matlab/','oetsettings.m']);
xbo=xb_read_output('Z:\Users\Tyler\projects\forecast\xbeach\croatan\');
xbo2=xbo.data(1).value;

date1=datetime('now');
delta_t=datenum(date1);
delta_t=delta_t-hours(8)-minutes(30);
delta_t2(1:234)=delta_t+minutes(1:20:4680);
delta_t3 = datestr(delta_t2);
delta_t3 = datetime(delta_t3); 

%%
 for ii=1
     H          =squeeze(xbo.data(17).value(ii,:));
     zs         =squeeze(xbo.data(13).value(ii,:));
     zb         =squeeze(xbo.data(9).value(ii,:));
     x          =xbo2.data(23).value';
     y          =squeeze(xbo2.data(24).value);
     t          =squeeze(xbo2.data(19).value);
 end
t_end = length(t);

%%






%%

vid_path=('Z:\Users\Tyler\projects\forecast\xbeach\croatan\');
cd(vid_path);
t_end=length(t);
x_begin2 = 75;
x_begin     =30;
H3          =squeeze(xbo.data(17).value(:,x_begin2))/.707;
zs2         =squeeze(xbo.data(13).value(:,x_begin));
zs3	    =round(zs2,4);
H4  	    =round(H3,4);

%zs4 	    =round(zs2*3.28084,3);
%H5 	    =round(H3*3.28084,3);

delta_t3.Format = 'yyyy-MM-dd HH:mm';
T = table(delta_t3(1:t_end,1),zs3(1:t_end,1),H4(1:t_end,1));
writetable(T,'tides.csv');
fence(1,1:486)=0;
fence(1,89:3:136)=3.25;
fence(1,89:3:101)=fence(1,92:3:104)-(zb(1,92:3:104))*.15;
fence(1,107:3:136)=fence(1,110:3:139)-(zb(1,110:3:139))*.275;
fence(1,104)=0;
veg_h(1,1:486)=0;
veg_h(1,127:300)=.15;

%%
% vid_path=('/home/admin/neptune_work/realtime/xbeach/croatan/figures/');
cd(vid_path);
video=VideoWriter('xbeach');

video.FrameRate=4;
open(video);
t_end=length(t);

for i=1:2%t_end
         if i<t_end
            H1          =squeeze(xbo.data(17).value(i,:));
            zs          =squeeze(xbo.data(13).value(i,:));
            zb1         =squeeze(xbo.data(9).value(i,:));
            x1          =xbo2.data(23).value';
            y1          =squeeze(xbo2.data(24).value);
	        sigma1      =squeeze(xbo.data(21).value(i,:));
	    
            Tp = (2*pi)./sigma1;
            alph = 0.0013;
            Depth       = zs-zb1;
            L = disper(Depth,Tp);
            Irib = tan(alph)./sqrt(H1/L);
 
            R = ((1.86*Irib)^0.71)*H1;
            zs1 = zs + R;
	        waves 	= (H1+zs1);
	baseline=-8;
	index=1:486;
    brown= [0.8 0.8 .6];
    blue2=[0 0.3 1];           

    A=figure('units','normalized','outerposition',[0 0 1 .7]);

        A(1)=plot(x1,waves,'c'); hold on
        A(2)=plot(x1,zs1,'b'); hold on
        A(3)=plot(x1,zb1,'k'); hold on

        A(4)=plot(x(1,89:3:101), zb(1,89:3:101)+fence(1,89:3:101)-1.25, '-o k');
        A(4)=plot(x(1,89:3:101), zb(1,89:3:101)+fence(1,89:3:101)-1.25, '-o k');
        A(4)=plot(x(1,89:3:101), zb(1,89:3:101)+fence(1,89:3:101)-1.5, '-o k');
        A(4)=plot(x(1,89:3:101), zb(1,89:3:101)+fence(1,89:3:101)-1.75, '-o k');
        A(4)=plot(x(1,89:3:101), zb(1,89:3:101)+fence(1,89:3:101)-2, '-o k');
        A(4)=plot(x(1,89:3:101), zb(1,89:3:101)+fence(1,89:3:101)-2.25, '-o k');
        A(4)=plot(x(1,89:3:101), zb(1,89:3:101)+fence(1,89:3:101)-2.5, '-o k');
        A(4)=plot(x(1,89:3:101), zb(1,89:3:101)+fence(1,89:3:101)-2.75, '-o k');
        A(4)=plot(x(1,89:3:101), zb(1,89:3:101)+fence(1,89:3:101)-3, '-o k');
        
        A(4)=plot(x(1,107:3:136), zb(1,107:3:136)+fence(1,107:3:136)-1.25, '-o k');
        A(4)=plot(x(1,107:3:136), zb(1,107:3:136)+fence(1,107:3:136)-1.5, '-o k');
        A(4)=plot(x(1,107:3:136), zb(1,107:3:136)+fence(1,107:3:136)-1.75, '-o k');
        A(4)=plot(x(1,107:3:136), zb(1,107:3:136)+fence(1,107:3:136)-2, '-o k');
        A(4)=plot(x(1,107:3:136), zb(1,107:3:136)+fence(1,107:3:136)-2.25, '-o k');
        A(4)=plot(x(1,107:3:136), zb(1,107:3:136)+fence(1,107:3:136)-2.5, '-o k');
        A(4)=plot(x(1,107:3:136), zb(1,107:3:136)+fence(1,107:3:136)-2.75, '-o k');
        for j = 1:length(x1)
            A(5)=plot([x1(1,j),x1(1,j)],[zb(1,j)-1 zb(1,j)+fence(1,j)-1],'Color',[.5 .4 .3]);
            set(A(5), 'LineWidth', 13);
        end
	    A(6)=fill(x1(index([ 1 1:end end])),...
               [baseline waves(index) baseline],...
                blue2,'EdgeColor','none'); alpha(.5)                     
        A(7)=fill(x1(index([ 1 1:end end])),...
               [baseline zs1(index) baseline],...
               'b','EdgeColor','none'); alpha(.5)                
        A(8)=fill(x1(index([1 1:end end])),...
               [baseline zb1(index) baseline],...
               brown,'EdgeColor','none');
	set(gcf,'color','w');
        set(A(1),'LineWidth',1.75);
        set(A(2),'LineWidth',1.75);
        set(A(3),'LineWidth',1.75);
        legend([A(1),A(2),A(3)],'waves','tide','bathymetry','Location','NorthWest');
        axis([460 562 -4 7]);	
	set(gca,'LooseInset',get(gca,'TightInset'));
        xlabel('Cross-Section Length (m)','FontSize',12);
        ylabel('Elevation above NAVD88 (m)','FontSize',12);
        str=strcat({'Virginia Beach, Croatan' ; datestr(delta_t2(i),' mm-dd-yyyy HH:MM')});
        title(str,'FontSize',16,'Fontweight','bold');
        F(i)=getframe(i);
        writeVideo(video,F(i));
	set(gcf,'color','w');
         else
             break;
         end
end    
close(video);



close all;