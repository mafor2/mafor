clear

%case4a
load ../output/concout.res

% SPECIES - LWC
smax=937;      %kpp_parameters NSPEC
smax=smax+1    %add 1 for time column
smax=smax+1;   %NSPEC+1 = LWC(1)

% TIME
time=concout(:,1);
tim=time/3600.;
tim=tim-5;
tmax=74*60*6;

% MIXING RATIO
tm=288.15;
NA=6.022E23;  %1/mol
R=83.145E5/NA;
M=101325./(R*tm);
ppb=M/1.E9;
ppt=M/1.E12;

% FROM AEROSOL MODE CS 
%LWCcs = concout(:,smax);
%LWCcs=3.53E-11;  %l/cm^3   (constant)
% FOR INCLOUD RUN
% average in fog:
LWCcs=7.0E-08;  
% prescribed LWC:
LWCcst = concout(:,smax);

% aq. species  ind_xxx_a01
ind_Hp_a04    = concout(:,925);
ind_H2SO4_a04 = concout(:,136);
ind_HSO4m_a04 = concout(:,685);
ind_SO4mm_a04 = concout(:,909);
ind_SO2_a04   = concout(:,651);
ind_OH_a04    = concout(:,934);
ind_HO2_a04   = concout(:,926);
ind_H2O2_a04  = concout(:,920);
ind_O3_a04    = concout(:,910);
ind_HONO_a04  = concout(:,919);
ind_DMA_a04   = concout(:,898);
ind_DMAp_a04  = concout(:,385);
ind_NDMA_a04  = concout(:,888);
%DMNNO2
ind_DMN_a04   = concout(:,880);  


%gas-phase
oh    = concout(:,915);
ho2   = concout(:,913);
o3    = concout(:,911);
h2o2  = concout(:,843);
hono  = concout(:,810);
dma   = concout(:,889);
ndma  = concout(:,501);
%DMNNO2
dnno2 = concout(:,453);


% unit conversion
% conversion molecules/cm^3(air) to mol/dm^3(water)
for i=1:tmax
 if (LWCcst(i)<1.e-10) 
   conv=0.;
 else
   %conv=1e3/(NA*LWCcst(i));
   conv=1e3/(NA*LWCcs);
 end
 pH_a04(i)        = -log10(ind_Hp_a04(i)*conv);
 ind_SVI_a04(i)   = (ind_H2SO4_a04(i)+ind_HSO4m_a04(i)+ind_SO4mm_a04(i))*conv;
 ind_H2SO4_a04(i) = ind_H2SO4_a04(i)*conv;
 ind_H2O2_a04(i)  = ind_H2O2_a04(i)*conv;
 ind_SO2_a04(i)   = ind_SO2_a04(i)*conv;
 ind_OH_a04(i)    = ind_OH_a04(i)*conv;
 ind_HO2_a04(i)   = ind_HO2_a04(i)*conv;
 ind_O3_a04(i)    = ind_O3_a04(i)*conv;
 ind_HONO_a04(i)  = ind_HONO_a04(i)*conv;
 ind_DMA_a04(i)   = ind_DMA_a04(i)*conv+ind_DMAp_a04(i)*conv;
 ind_NDMA_a04(i)  = ind_NDMA_a04(i)*conv;
 ind_DMN_a04(i)   = ind_DMN_a04(i)*conv;

 %gas-phase
 ind_O3(i)        = o3(i)/ppb;
 ind_OH(i)        = oh(i);
 ind_HO2(i)       = ho2(i);
 ind_H2O2(i)      = h2o2(i)/ppt;
 ind_HONO(i)      = hono(i)/ppt;
 ind_DMA(i)       = dma(i)/ppb;
 ind_NDMA(i)      = ndma(i)/ppt;
 ind_DMN(i)       = dnno2(i)/ppt;
end

%set for getting the transparency right in Octave
%comment below line if using MATLAB
graphics_toolkit ("gnuplot")

figure(1);clf
axes('linewidth',2,'fontsize',20)
hold
ax=gca;
set(ax,'linewidth',2,'fontsize',15)
XD=4:12:74;

h(1)=subplot(2,5,1); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,3.0,3.0,0,0,3.0,3.0,0,0,3.0,3.0,0,0,3.0,3.0,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:360:tmax),ind_HONO_a04(1:360:tmax)*1e8,'k-','linewidth',2.6);
title('HONOaq (1.e-8M)')
set(h(1),'xtick',XD);
set(h(1),'XTickLabel',{'08','20','08','20','08','20'})
set(h(1),'xlim', [2,74]) 
set(h(1),'ylim', [0.0,3.0]) 

h(2)=subplot(2,5,2); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,25.0,25.0,0,0,25.0,25.0,0,0,25.0,25.0,0,0,25.0,25.0,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:tmax),ind_O3_a04(1:tmax)*1e10,'k-','linewidth',2.6);
title('O3aq (1.e-10M)')
set(h(2),'xtick',XD);
set(h(2),'XTickLabel',{'08','20','08','20','08','20'})
set(h(2),'xlim', [2,74]) 

h(3)=subplot(2,5,3); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,4.0,4.0,0,0,4.0,4.0,0,0,4.0,4.0,0,0,4.0,4.0,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:360:tmax),ind_H2O2_a04(1:360:tmax)*1e4,'k-','linewidth',2.6);
title('H2O2aq (1.e-4M)')
set(h(3),'xtick',XD);
set(h(3),'XTickLabel',{'08','20','08','20','08','20'})
set(h(3),'xlim', [2,74]) 

h(4)=subplot(2,5,4); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,10,10,0,0,10,10,0,0,10,10,0,0,10,10,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:360:tmax),ind_OH_a04(1:360:tmax)*1e12,'k-','linewidth',2.6);
title('OHaq (1.e-12M)')
set(h(4),'xtick',XD);
set(h(4),'XTickLabel',{'08','20','08','20','08','20'})
set(h(4),'xlim', [2,74]) 

h(5)=subplot(2,5,5); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,10.0,10.0,0,0,10.0,10.0,0,0,10.0,10.0,0,0,10.0,10.0,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:360:tmax),ind_HO2_a04(1:360:tmax)*1e8,'k-','linewidth',2.6);
title('HO2aq (1.e-8M)')
set(h(5),'xtick',XD);
set(h(5),'XTickLabel',{'08','20','08','20','08','20'})
set(h(5),'xlim', [2,74]) 

h(6)=subplot(2,5,6); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,400,400,0,0,400,400,0,0,400,400,0,0,400,400,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:tmax),ind_HONO(1:tmax),'k-','linewidth',2.6);
title('HONOg (ppt)')
set(h(6),'xtick',XD);
set(h(6),'XTickLabel',{'08','20','08','20','08','20'})
set(h(6),'xlim', [2,74]) 

%gas-phase
h(7)=subplot(2,5,7); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,100,100,0,0,100,100,0,0,100,100,0,0,100,100,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:tmax),ind_O3(1:tmax),'k-','linewidth',2.6);
title('O3g (ppb)')
set(h(7),'xtick',XD);
set(h(7),'XTickLabel',{'08','20','08','20','08','20'})
set(h(7),'xlim', [2,74]) 

h(8)=subplot(2,5,8); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,2000,2000,0,0,2000,2000,0,0,2000,2000,0,0,2000,2000,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:360:tmax),ind_H2O2(1:360:tmax),'k-','linewidth',2.6); 
title('H2O2g (ppt)')
set(h(8),'xtick',XD);
set(h(8),'XTickLabel',{'08','20','08','20','08','20'})
set(h(8),'xlim', [2,74]) 

h(9)=subplot(2,5,9); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,20,20,0,0,20,20,0,0,20,20,0,0,20,20,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:tmax),ind_OH(1:tmax)*1e-6,'k-','linewidth',2.6);
title('OHg (1E6cm-3)')
set(h(9),'xtick',XD);
set(h(9),'XTickLabel',{'08','20','08','20','08','20'})
set(h(9),'xlim', [2,74]) 

h(10)=subplot(2,6,12); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,6,6,0,0,6,6,0,0,6,6,0,0,6,6,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:tmax),ind_HO2(1:tmax)*1e-8,'k-','linewidth',2.6);
title('HO2g (1E8cm^{-3})')
set(h(10),'xtick',XD);
set(h(10),'XTickLabel',{'08','20','08','20','08','20'})
set(h(10),'xlim', [2,74]) 


%FaceAlpha (transparency) is not supported in Octave
%because it is not supported in Postscript

%MATLAB
%print -dbmp '../afigs_aqchem/marine_aqchem_base.bmp'
%OCTAVE
print -djpg '../afigs_aqchem/marine_aqchem_base.jpg'
