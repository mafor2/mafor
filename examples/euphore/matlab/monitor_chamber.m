clear
load ../data/monitor_20090511.dat;
load ../data/ptrms_20090511.dat;
load ../data/smpsams_20090511.dat;
load ../output/concout.res;
load ../output/aerconc.res;


time=concout(:,1);
time=time/3600.;
time=time-2;  % convert to UT
timeaer=concout(:,1);
timeaer=timeaer/3600.;
timeaer=timeaer-2;
tmax=6*60*5;

TM=300.
R=83.145E5/6.0221E23
M=101325./(R*TM)
ppb=M/1.E9
ppt=M/1.e12

o3mod  =concout(:,911);
nomod  =concout(:,932);
no2mod =concout(:,927);
hono   =concout(:,810);
ohmod  =concout(:,915);
hno3mod=concout(:,877);
comod  =concout(:,903);
mea    =concout(:,401);
%LTMB   
tmb    =concout(:,312);  
%H2NCHO
formd  =concout(:,607);   
%MEANNO2
nitmd  =concout(:,291);   

o3mod=o3mod/ppb;
nomod=nomod/ppb;
no2mod=no2mod/ppb;
honomod=hono/ppb;
hno3mod=hno3mod/ppb;
meamod=mea/ppb;
tmbmod=tmb/ppb;
formod=formd/ppb;
nitmod=nitmd/ppb;
comod=comod/ppb;


% ORG and MEA condensation
pmorgct=aerconc(:,18)+aerconc(:,19)+aerconc(:,20)+aerconc(:,21);
pmammot=aerconc(:,22)+aerconc(:,23)+aerconc(:,24)+aerconc(:,25);
pmnitrt=aerconc(:,26)+aerconc(:,27)+aerconc(:,28)+aerconc(:,29);
pmecbct=aerconc(:,30)+aerconc(:,31)+aerconc(:,32)+aerconc(:,33);
pmorgt =pmorgct*1.e-3;
%pmxxxt =pmxxxxt*1.e-3; 
pmno3t =pmnitrt*1.e-3;
pmnamt =pmammot*1.e-3;
pmbcet =pmecbct*1.e-3; 
% TOTAL in ug/m3
pmtotal=(pmecbct+pmammot+pmnitrt+pmorgct)*1e-3;


%%%TIME        O3        NO        NO2        JNO2
tmeas   =monitor_20090511(:,1)/3600;
o3meas  =monitor_20090511(:,2);
nomeas  =monitor_20090511(:,3);
no2meas =monitor_20090511(:,4);
jno2meas=monitor_20090511(:,5);

%%%TIME     DMA   MeMA  MeFo   NDMA  DMNNO2
tptr    =ptrms_20090511(:,1)/3600.;
meaptr  =ptrms_20090511(:,2);
forptr  =ptrms_20090511(:,3);
meaerr   =meaptr*0.15;

%%% Total Aerosol Mass (AMS) in ug/m3
timea    = smpsams_20090511(:,1);
pmtotct  = smpsams_20090511(:,6);
timea=timea/3600.;

figure(1)
[AX,H1,H2] = plotyy(tptr,meaptr,timea,pmtotct,'plot');
hold
ax=gca;
set(ax,'linewidth',1.5,'fontsize',13)
plot(tptr,forptr, 'rs','linewidth',1.5,'MarkerSize',4.1)
plot(time,ohmod*5e-5,'Color',[.8 .8 .8],'LineWidth',2.1)
plot(time,meamod ,'-k','LineWidth',1.6)
plot(time,formod,':r','LineWidth',2.6)
%plot(time,nitmod*50,'k--','LineWidth',1.6)
plot(timeaer,pmtotal*0.5,'b--','LineWidth',1.6)
%errorbar(tptr,meaptr,meaerr,'Color',[.5 .5 .5],'Marker','none','LineStyle','none','linewidth',1.,'MarkerSize',4.1)
line([10.80 10.80],[0 1200],'LineStyle','--','linewidth',2.,'Color','k')
line([14.65 14.65],[0 1200],'LineStyle','--','linewidth',2.,'Color','k')
set(AX(1),'linewidth',1.5,'fontsize',15.,'YColor',[0 0 0])
set(AX(2),'linewidth',1.5,'fontsize',15.,'YColor',[0 0 0])
set(get(AX(1),'Ylabel'),'String','Mixing ratio (ppbv)','FontSize',12,'Color',[0 0 0],'FontName','Arial')
set(get(AX(2),'Ylabel'),'String','Tot. Aerosol Mass (\mug m^{-3})','FontSize',10,'Color',[0 0 0],'FontName','Arial') 
title('11.05.2009 MEA+OH Exp','FontSize',20)
legend('MEA PTR-MS','H2NCHO PTR-MS','OH model','MEA model','H2NCHO model','Tot. Aerosol model','Location','NorthEast')
xlabel('Time UTC (hrs)')
tx=10.0:1:15.0;
set(AX(1),'xtick',tx);
set(AX(2),'xtick',tx);
set(AX(1),'XTickLabel',{'10:00','11:00','12:00','13:00','14:00','15:00'});
set(AX(2),'XTickLabel',{'10:00','11:00','12:00','13:00','14:00','15:00'});
set(AX(1),'Xlim',[10 15]); 
set(AX(2),'Xlim',[10 15]);
set(AX(1),'Ylim',[0 140]); 
set(AX(1),'ytick',[0,40,80,120,160,200]);
set(AX(2),'Ylim',[0 280]); 
set(AX(2),'ytick',[0,50,100,150,200,250,300]);
set(H1,'Color',[0 0 0],'Marker','o','LineStyle','none','linewidth',2.5,'MarkerSize',6.1)
set(H2,'LineStyle','--' ,'Color',[0 0 1],'LineWidth',2.5,'Marker','o','MarkerSize',6.1)
%uncomment below line to save as jpg
print -djpg '../afigs_chamber/2009-05-11-mea1_mod.jpg'

figure(2)
[AX,H1,H2] = plotyy(tmeas,o3meas,tmeas,jno2meas*1e3,'plot');
hold
ax=gca;
set(ax,'linewidth',1.5,'fontsize',16)
plot(tmeas,nomeas, 'gs','linewidth',1.5,'MarkerSize',4.1)
plot(tmeas,no2meas,'ms','linewidth',1.5,'MarkerSize',4.1)
plot(time(1:tmax),o3mod(1:tmax), '-k','LineWidth',2.5)
plot(time(1:tmax),nomod(1:tmax), 'Color',[.1 .45 .1],'LineWidth',2.5)
plot(time(1:tmax),no2mod(1:tmax),'Color',[.45 .1 .45],'LineWidth',2.5)
%plot(time(1:tmax),honomod(1:tmax)*10,'-r','LineWidth',3)
%plot(time(1:tmax),ohmod(1:tmax)*1e-5, '--b','LineWidth',2.5)
line([10.80 10.80],[0 1200],'LineStyle','--','linewidth',2.,'Color','k')
line([14.65 14.65],[0 1200],'LineStyle','--','linewidth',2.,'Color','k')
set(AX(1),'linewidth',1.5,'fontsize',15.,'YColor',[0 0 0])
set(AX(2),'linewidth',1.5,'fontsize',15.,'YColor',[0 0 0])
set(get(AX(1),'Ylabel'),'String','O_3, NO, NO_2 (ppb)','FontSize',14,'Color',[0 0 0],'FontName','Arial')
set(get(AX(2),'Ylabel'),'String','Photol. Frequency j(NO_2)*1E3 (s^{-1})','FontSize',14,'Color',[0 0 0],'FontName','Arial') 
xlabel('Time UTC (hrs)','FontSize',14)
title('11.05.2009 MEA+OH Exp','FontSize',20)
legend('O3 meas.','NO meas.','NO2 meas.','O3 mod','NO mod','NO2 mod','Location','NorthWest')
tx=10.0:1:15.0;
set(AX(1),'xtick',tx);
set(AX(2),'xtick',tx);
set(AX(1),'XTickLabel',{'10:00','11:00','12:00','13:00','14:00','15:00'});
set(AX(2),'XTickLabel',{'10:00','11:00','12:00','13:00','14:00','15:00'});
set(AX(1),'Xlim',[10 15]); 
set(AX(2),'Xlim',[10 15]);
set(AX(1),'Ylim',[0 1050]); 
set(AX(1),'ytick',[0,200,400,600,800,1000]);
set(AX(2),'Ylim',[0 15]); 
set(AX(2),'ytick',[0,3,6,9,12,15]);
set(H1,'Color',[0 0 1],'Marker','s','LineStyle','none','linewidth',1.5,'MarkerSize',4.1)
set(H2,'LineStyle','-' ,'Color',[.6 .6 .6],'LineWidth',4.1)
%uncomment below line to save as jpg
print -djpg '../afigs_chamber/2009-05-11-photo1.jpg'

