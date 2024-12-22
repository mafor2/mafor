clear

load ../data/ptrms_20090511.dat;
load ../data/smpsams_20090511.dat;

load ../output/aerconc.res;
load ../output/soadis.res;

time=aerconc(:,1);
time=time/3600.;
time=time-2;  % convert to UT
timeaer=aerconc(:,1);
timeaer=timeaer/3600.;
timeaer=timeaer-2;
tmax=6*60*5;

TM=300.
R=83.145E5/6.0221E23
M=101325./(R*TM)
ppb=M/1.E9
ppt=M/1.e12

% soadis.res output
% all in ug/m3
% 2-10  C*   (SOA1 ... SOA9)
% 11-19 Cgas
% 20-28 SOA
soagas1=soadis(:,11);
soagas2=soadis(:,12);
soagas3=soadis(:,13);
soapar1=soadis(:,20);
soapar2=soadis(:,21);
soapar3=soadis(:,22);
soatot=soapar1+soapar2+soapar3;

% ORG and MEA condensation
pmorgct=aerconc(:,22)+aerconc(:,23)+aerconc(:,24)+aerconc(:,25);
pmammot=aerconc(:,27)+aerconc(:,28)+aerconc(:,29)+aerconc(:,30);
pmnitrt=aerconc(:,32)+aerconc(:,33)+aerconc(:,34)+aerconc(:,35);
pmecbct=aerconc(:,37)+aerconc(:,38)+aerconc(:,39)+aerconc(:,40);

pmorgt =pmorgct*1.e-3;
%pmxxxt =pmxxxxt*1.e-3;
pmno3t =pmnitrt*1.e-3;
pmnamt =pmammot*1.e-3;
pmbcet =pmecbct*1.e-3;
% TOTAL in ug/m3
pmtotal=(pmecbct+pmammot+pmnitrt+pmorgct)*1e-3;
% total OC in ug/m3
pmsoat = pmorgct*1e-3;
pmamin = pmammot*1e-3;

%%%TIME     DMA   MeMA  MeFo   NDMA  DMNNO2
tptr      =ptrms_20090511(:,1)/3600.;
meaptr    =ptrms_20090511(:,2);
meaerr    =meaptr*0.15;

%%% Total Aerosol Mass (AMS) in ug/m3
timea      =smpsams_20090511(:,1);
pmtotct    =smpsams_20090511(:,6);
%%% AMS
totmassams =smpsams_20090511(:,6);
nitrateams =smpsams_20090511(:,7);
organicams =smpsams_20090511(:,8);
ammoniaams =smpsams_20090511(:,9);

timea=timea/3600.;

fsize=12;

figure(1)
[AX,H1,H2] = plotyy(timea,organicams,timea,organicams,'plot');
hold
ax=gca;
set(ax,'linewidth',1.5,'fontsize',fsize);
plot(timeaer,soagas1*1.0 ,'--b','LineWidth',1.8)
plot(timeaer,soagas2*1.0 ,'--g','LineWidth',1.8)
plot(timeaer,soagas3*1.0 ,'--r','LineWidth',1.8)
plot(timeaer,soapar1*1.0 , '-b','LineWidth',1.8)
plot(timeaer,soapar2*1.0 , '-g','LineWidth',1.8)
plot(timeaer,soapar3*1.0 , '-r','LineWidth',1.8)
plot(timeaer,soatot *1.0 , 'Color',[.6 .6 .6],'LineWidth',2.6)
plot(timeaer,pmsoat *1.0, 'k--','LineWidth',2.0)
%errorbar(tptr,meaptr,meaerr,'Color',[.5 .5 .5],'Marker','none','LineStyle','none','linewidth',1.,'MarkerSize',4.1)
line([10.80 10.80],[0 1200],'LineStyle','--','linewidth',2.,'Color','k')
line([14.65 14.65],[0 1200],'LineStyle','--','linewidth',2.,'Color','k')
set(AX(1),'linewidth',1.5,'fontsize',fsize,'YColor',[0 0 0])
set(AX(2),'linewidth',1.5,'fontsize',fsize,'YColor',[0 0 0])
set(get(AX(1),'Ylabel'),'String','Mass concentration (\mug m^{-3})','FontSize',fsize,'Color',[0 0 0],'FontName','Arial')
set(get(AX(2),'Ylabel'),'String','Organic Aerosol Mass (\mug m^{-3})','FontSize',fsize,'Color',[0 0 0],'FontName','Arial')
title('11.05.2009 MEA+OH Exp','FontSize',fsize)
legend('OC AMS','SOA1(g) mod','SOA2(g) mod','SOA3(g) mod', ...
                  'SOA1(p) mod','SOA2(p) mod','SOA3(p) mod', ...
                  'SOA mod','OA mod','Location','NorthWest')
xlabel('Time UTC (hrs)')
tx=10.0:1:15.0;
set(AX(1),'tickdir','out');
set(AX(1),'xtick',tx);
set(AX(2),'xtick',tx);
set(AX(1),'XTickLabel',{'10:00','11:00','12:00','13:00','14:00','15:00'});
set(AX(2),'XTickLabel',{'10:00','11:00','12:00','13:00','14:00','15:00'});
set(AX(1),'Xlim',[10 15]);
set(AX(2),'Xlim',[10 15]);
set(AX(1),'Ylim',[0 82]);
set(AX(1),'ytick',[0,10,20,30,40,50,60,70,80,90,100]);
set(AX(2),'tickdir','out');
set(AX(2),'Ylim',[0 82]);
set(AX(2),'ytick',[0,10,20,30,40,50,60,70,80,90,100]);
set(H1,'Color',[0 0 0],'Marker','o','linewidth',1.5,'MarkerSize',4.1)
set(H2,'LineStyle','--' ,'Color',[0 0 0],'LineWidth',2.0,'MarkerSize',5.1)
%uncomment below line to save as jpg
print -djpg '../afigs_chamber/2009-05-11-soadis_mod.jpg'
