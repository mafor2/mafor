clear

%case4a
load ../output/concout.res

% SPECIES - LWC
smax=1254;      %kpp_parameters NSPEC
smax=smax+1    %add 1 for time column
%%smax=smax+1;   %NSPEC+1 = LWC(1)
smax=smax+3;   %NSPEC+3 = LWC(3) CS mode

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
LWCcs=1.19E-07;
% prescribed LWC:
LWCcst = concout(:,smax);

% aq. species  ind_xxx_a03
ind_Hp_a03    = concout(:, 1237);
ind_H2SO4_a03 = concout(:,  186);
ind_HSO4m_a03 = concout(:,  646);
ind_SO4mm_a03 = concout(:, 1162);
ind_SO2_a03   = concout(:,  830);
ind_OH_a03    = concout(:, 1197);
ind_HO2_a03   = concout(:, 1227);
ind_H2O2_a03  = concout(:, 1233);
ind_O3_a03    = concout(:, 1240);
ind_HONO_a03  = concout(:, 1245);
ind_DMA_a03   = concout(:, 1167);
ind_DMAp_a03  = concout(:,  490);
ind_NDMA_a03  = concout(:, 1120);
%DMNNO2
ind_DMN_a03   = concout(:, 1110);


%gas-phase
oh    = concout(:, 1191);
ho2   = concout(:, 1199);
o3    = concout(:, 1238);
h2o2  = concout(:, 1061);
hono  = concout(:, 1094);
dma   = concout(:,  884);
ndma  = concout(:,  966);
%DMNNO2
dnno2 = concout(:,  942);


% unit conversion
% conversion molecules/cm^3(air) to mol/dm^3(water)
for i=1:tmax
 if (LWCcst(i)<1.e-10)
   conv=0.;
 else
   %conv=1e3/(NA*LWCcst(i));
   conv=1e3/(NA*LWCcs);
 end
 pH_a04(i)        = -log10(ind_Hp_a03(i)*conv);
 ind_SVI_a03(i)   = (ind_H2SO4_a03(i)+ind_HSO4m_a03(i)+ind_SO4mm_a03(i))*conv;
 ind_H2SO4_a03(i) = ind_H2SO4_a03(i)*conv;
 ind_H2O2_a03(i)  = ind_H2O2_a03(i)*conv;
 ind_SO2_a03(i)   = ind_SO2_a03(i)*conv;
 ind_OH_a03(i)    = ind_OH_a03(i)*conv;
 ind_HO2_a03(i)   = ind_HO2_a03(i)*conv;
 ind_O3_a03(i)    = ind_O3_a03(i)*conv;
 ind_HONO_a03(i)  = ind_HONO_a03(i)*conv;
 ind_DMA_a03(i)   = ind_DMA_a03(i)*conv+ind_DMAp_a03(i)*conv;
 ind_NDMA_a03(i)  = ind_NDMA_a03(i)*conv;
 ind_DMN_a03(i)   = ind_DMN_a03(i)*conv;

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
%graphics_toolkit ("gnuplot")
%commented in Octave v9.2.0:
%warning: using the gnuplot graphics toolkit is discouraged

fsize=7;

figure(1);clf
axes('linewidth',2,'fontsize',fsize)
hold
ax=gca;
set(ax,'linewidth',2,'fontsize',fsize,'tickdir','out')
XD=4:12:74;

h(1)=subplot(2,5,1); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,4.0,4.0,0,0,4.0,4.0,0,0,4.0,4.0,0,0,4.0,4.0,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:360:tmax),ind_HONO_a03(1:360:tmax)*1e8,'k-','linewidth',1.0);
title('HONOaq (1.e-8M)','fontsize',fsize)
set(h(1),'xtick',XD);
set(h(1),'XTickLabel',{'8','20','8','20','8','20'})
set(h(1),'xlim', [2,74])
set(h(1),'ylim', [0.0,3.0])
set(h(1),'linewidth',1,'fontsize',fsize)

h(2)=subplot(2,5,2); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,25.0,25.0,0,0,25.0,25.0,0,0,25.0,25.0,0,0,25.0,25.0,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:tmax),ind_O3_a03(1:tmax)*1e10,'k-','linewidth',0.9);
title('O3aq (1.e-10M)','fontsize',fsize)
set(h(2),'xtick',XD);
set(h(2),'XTickLabel',{'8','20','8','20','8','20'})
set(h(2),'xlim', [2,74])
set(h(2),'linewidth',1,'fontsize',fsize)

h(3)=subplot(2,5,3); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,4.0,4.0,0,0,4.0,4.0,0,0,4.0,4.0,0,0,4.0,4.0,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:360:tmax),ind_H2O2_a03(1:360:tmax)*1e4,'k-','linewidth',1.0);
title('H2O2aq (1.e-4M)','fontsize',fsize)
set(h(3),'xtick',XD);
set(h(3),'XTickLabel',{'8','20','8','20','8','20'})
set(h(3),'xlim', [2,74])
set(h(3),'linewidth',1,'fontsize',fsize)

h(4)=subplot(2,5,4); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,10.0,10.0,0,0,10.0,10.0,0,0,10.0,10.0,0,0,10.0,10.0,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:360:tmax),ind_OH_a03(1:360:tmax)*1e13,'k-','linewidth',1.0);
title('OHaq (1.e-13M)','fontsize',fsize)
set(h(4),'xtick',XD);
set(h(4),'XTickLabel',{'8','20','8','20','8','20'})
set(h(4),'xlim', [2,74])
set(h(4),'linewidth',1,'fontsize',fsize)

h(5)=subplot(2,5,5); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,10.0,10.0,0,0,10.0,10.0,0,0,10.0,10.0,0,0,10.0,10.0,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:360:tmax),ind_HO2_a03(1:360:tmax)*1e8,'k-','linewidth',1.0);
title('HO2aq (1.e-8M)','fontsize',fsize)
set(h(5),'xtick',XD);
set(h(5),'XTickLabel',{'8','20','8','20','8','20'})
set(h(5),'xlim', [2,74])
set(h(5),'linewidth',1,'fontsize',fsize)

%gas-phase

h(6)=subplot(2,5,6); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,400,400,0,0,400,400,0,0,400,400,0,0,400,400,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:tmax),ind_HONO(1:tmax),'k-','linewidth',1.0);
title('HONOg (ppt)','fontsize',fsize)
set(h(6),'xtick',XD);
set(h(6),'XTickLabel',{'8','20','8','20','8','20'})
set(h(6),'xlim', [2,74])
set(h(6),'linewidth',1,'fontsize',fsize)

h(7)=subplot(2,5,7); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,100,100,0,0,100,100,0,0,100,100,0,0,100,100,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:tmax),ind_O3(1:tmax),'k-','linewidth',1.0);
title('O3g (ppb)','fontsize',fsize)
set(h(7),'xtick',XD);
set(h(7),'XTickLabel',{'8','20','8','20','8','20'})
set(h(7),'xlim', [2,74])
set(h(7),'linewidth',1,'fontsize',fsize)

h(8)=subplot(2,5,8); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,2000,2000,0,0,2000,2000,0,0,2000,2000,0,0,2000,2000,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:360:tmax),ind_H2O2(1:360:tmax),'k-','linewidth',1.0);
title('H2O2g (ppt)','fontsize',fsize)
set(h(8),'xtick',XD);
set(h(8),'XTickLabel',{'8','20','8','20','8','20'})
set(h(8),'xlim', [2,74])
set(h(8),'linewidth',1,'fontsize',fsize)

h(9)=subplot(2,5,9); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,20,20,0,0,20,20,0,0,20,20,0,0,20,20,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:tmax),ind_OH(1:tmax)*1e-6,'k-','linewidth',1.0);
title('OHg (1E6cm^{-3})','fontsize',fsize)
set(h(9),'xtick',XD);
set(h(9),'XTickLabel',{'8','20','8','20','8','20'})
set(h(9),'xlim', [2,74])
set(h(9),'linewidth',1,'fontsize',fsize)

h(10)=subplot(2,5,10); fill([2,2,6,6,21,21,30,30,44,44,54,54,69,69,74,74], ...
    [0,7,7,0,0,7,7,0,0,7,7,0,0,7,7,0],'b','FaceAlpha',0.4,'LineStyle','none')
hold on;
plot(tim(1:tmax),ind_HO2(1:tmax)*1e-8,'k-','linewidth',1.0);
title('HO2g (1E8cm^{-3})','fontsize',fsize)
set(h(10),'xtick',XD);
set(h(10),'XTickLabel',{'8','20','8','20','8','20'})
set(h(10),'xlim', [2,74])
set(h(10),'linewidth',1,'fontsize',fsize)

%MATLAB
%print -dbmp '../afigs_aqchem/marine_aqchem_base.bmp'
%OCTAVE
print -djpg '../afigs_aqchem/marine_aqchem_base.jpg'

