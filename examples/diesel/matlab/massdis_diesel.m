clear

load ../output/wetdp.res;
load ../output/size_dism.res;
load ../output/size_dis.res;

infiledp='wetdp.res';
indp=strrep(infiledp,'.res','');
dp=eval(indp);
[rowdp,coldp]=size(dp);

infile='size_dis.res';
in=strrep(infile,'.res','');
ys=eval(in);
[rows,cols]=size(ys);  

% 0.0s
diameter_bin01=dp(1,2:coldp)*1e9; 
% 0.1s
diameter_bin02=dp(1+1,2:coldp)*1e9; 
% 0.9s
diameter_bin03=dp(1+9,2:coldp)*1e9;
% 2.7s
diameter_bin04=dp(1+27,2:coldp)*1e9; 

% Mass distribution
infileam='size_dism.res';
inam=strrep(infileam,'.res','');
yam=eval(inam);
[row,col]=size(yam);                    %row=xxx col=61

a=1;
b=22;   %  1 x 21 + 1
c=190;  %  9 x 21 + 1 
d=568;  % 27 x 21 + 1

%TOTAL MASS kg/m3-->ng/m3
dmdlogdp_bin01=yam(a,2:col)   *1e12  *2.303;
dmdlogdp_bin02=yam(b,2:col)   *1e12  *2.303;
dmdlogdp_bin03=yam(c,2:col)   *1e12  *2.303; 
dmdlogdp_bin04=yam(d,2:col)   *1e12  *2.303; 

%Sulfate ng/m3
dmdlogdp_bin01_sul=yam(a+1,2:col)    *2.303;
dmdlogdp_bin02_sul=yam(b+1,2:col)    *2.303; 
dmdlogdp_bin03_sul=yam(c+1,2:col)    *2.303; 
dmdlogdp_bin04_sul=yam(d+1,2:col)    *2.303; 

%SOA-2 ng/m3 (BLOV)
dmdlogdp_bin01_adi=yam(a+7,2:col)    *2.303; 
dmdlogdp_bin02_adi=yam(b+7,2:col)    *2.303; 
dmdlogdp_bin03_adi=yam(c+7,2:col)    *2.303; 
dmdlogdp_bin04_adi=yam(d+7,2:col)    *2.303; 

%SOA-9 ng/m3 (PELV)
dmdlogdp_bin01_elv=yam(a+14,2:col)   *2.303;
dmdlogdp_bin02_elv=yam(b+14,2:col)   *2.303; 
dmdlogdp_bin03_elv=yam(c+14,2:col)   *2.303; 
dmdlogdp_bin04_elv=yam(d+14,2:col)   *2.303; 

%SOOT
dmdlogdp_bin01_ebc=yam(a+17,2:col)   *2.303; 
dmdlogdp_bin02_ebc=yam(b+17,2:col)   *2.303; 
dmdlogdp_bin03_ebc=yam(c+17,2:col)   *2.303; 
dmdlogdp_bin04_ebc=yam(d+17,2:col)   *2.303; 

%POM
dmdlogdp_bin01_pom=yam(a+19,2:col)   *2.303; 
dmdlogdp_bin02_pom=yam(b+19,2:col)   *2.303; 
dmdlogdp_bin03_pom=yam(c+19,2:col)   *2.303; 
dmdlogdp_bin04_pom=yam(d+19,2:col)   *2.303; 

%H2O
dmdlogdp_bin01_wat=yam(a+20,2:col)   *2.303; 
dmdlogdp_bin02_wat=yam(b+20,2:col)   *2.303; 
dmdlogdp_bin03_wat=yam(c+20,2:col)   *2.303; 
dmdlogdp_bin04_wat=yam(d+20,2:col)   *2.303;

figure(1);clf
axes('linewidth',1.6,'fontsize',20)
hold
ax=gca;
set(ax,'linewidth',1.6,'fontsize',18)

%TIME 0.0 s
h(1)=subplot(2,2,1); loglog(diameter_bin01,dmdlogdp_bin01,      '-g' ,'LineWidth',2.4); hold on; ...
    loglog(diameter_bin01,dmdlogdp_bin01_ebc,  '--k' ,'LineWidth',2.1); ...
    loglog(diameter_bin01,dmdlogdp_bin01_elv,  'k:' ,'LineWidth',2.6); ...
    loglog(diameter_bin01,dmdlogdp_bin01_sul,  '-.k' ,'LineWidth',1.8); ...
    loglog(diameter_bin01,dmdlogdp_bin01_adi,  '--ko' ,'LineWidth',1.1,'MarkerSize',3.0);...
    loglog(diameter_bin01,dmdlogdp_bin01_pom,  '--ro' ,'LineWidth',1.1,'MarkerSize',3.0);...
    loglog(diameter_bin01,dmdlogdp_bin01_wat,  '--b'  ,'LineWidth',2.1,'MarkerSize',3.0);
title('Time 0.0s','FontSize',12)
%xlabel('D_p (nm)','FontSize',12)
[g]=legend('Total','Soot','OM_l','Sulfate','OM_s','OM nv','H2O','Location','southeastoutside');
set(g,'fontsize',9);
ylabel('dM/dlogDp (ng/m3)','FontSize',12)
set(h(1),'xtick',[1,10,100,1000]);
set(h(1), 'xticklabel', [1,10,100,1000]);
set(h(1),'XLim',[1. 600.],'Ylim',[1.e1 1.e7]);

%TIME 0.1 s
h(2)=subplot(2,2,2); loglog(diameter_bin02,dmdlogdp_bin02,      '-g' ,'LineWidth',2.4); hold on; ...
    loglog(diameter_bin02,dmdlogdp_bin02_ebc,  '--k' ,'LineWidth',2.1); ...
    loglog(diameter_bin02,dmdlogdp_bin02_elv,  'k:' ,'LineWidth',2.6); ...
    loglog(diameter_bin02,dmdlogdp_bin02_sul,  '-.k' ,'LineWidth',1.8); ...
    loglog(diameter_bin02,dmdlogdp_bin02_adi,  '--ko' ,'LineWidth',1.1,'MarkerSize',3.0);...
    loglog(diameter_bin02,dmdlogdp_bin02_pom,  '--ro' ,'LineWidth',1.1,'MarkerSize',3.0);...
    loglog(diameter_bin02,dmdlogdp_bin02_wat,  '--b'  ,'LineWidth',2.1,'MarkerSize',3.0);
title('Time 0.1s','FontSize',12)
%xlabel('D_p (nm)','FontSize',12)
%ylabel('dM/dlogDp (ng/m3)','FontSize',12)
set(h(2),'xtick',[1,10,100,1000]);
set(h(2), 'xticklabel', [1,10,100,1000]);
set(h(2),'XLim',[1. 600.],'Ylim',[1.e1 1.e7]); 

%TIME 0.9 s
h(3)=subplot(2,2,3); loglog(diameter_bin03,dmdlogdp_bin03,      '-g' ,'LineWidth',2.4); hold on; ...
    loglog(diameter_bin03,dmdlogdp_bin03_ebc,  '--k' ,'LineWidth',2.1); ...
    loglog(diameter_bin03,dmdlogdp_bin03_elv,  'k:' ,'LineWidth',2.6); ...
    loglog(diameter_bin03,dmdlogdp_bin03_sul,  '-.k' ,'LineWidth',1.8); ...
    loglog(diameter_bin03,dmdlogdp_bin03_adi,  '--ko' ,'LineWidth',1.1,'MarkerSize',3.0);...
    loglog(diameter_bin03,dmdlogdp_bin03_pom,  '--ro' ,'LineWidth',1.1,'MarkerSize',3.0);...
    loglog(diameter_bin03,dmdlogdp_bin03_wat,  '--b'  ,'LineWidth',2.1,'MarkerSize',3.0);
title('Time 0.9s','FontSize',12)
xlabel('D_p (nm)','FontSize',12)
ylabel('dM/dlogDp (ng/m3)','FontSize',12)
set(h(3),'xtick',[1,10,100,1000]);
set(h(3), 'xticklabel', [1,10,100,1000]);
set(h(3),'XLim',[1. 600.],'Ylim',[1.e1 1.e7]); 

%TIME 2.7 s
h(4)=subplot(2,2,4); loglog(diameter_bin04,dmdlogdp_bin04,      '-g' ,'LineWidth',2.4); hold on; ...
    loglog(diameter_bin04,dmdlogdp_bin04_ebc,  '--k' ,'LineWidth',2.1); ...
    loglog(diameter_bin04,dmdlogdp_bin04_elv,  'k:' ,'LineWidth',2.6); ...
    loglog(diameter_bin04,dmdlogdp_bin04_sul,  '-.k' ,'LineWidth',1.8); ...
    loglog(diameter_bin04,dmdlogdp_bin04_adi,  '--ko' ,'LineWidth',1.1,'MarkerSize',3.0);...
    loglog(diameter_bin04,dmdlogdp_bin04_pom,  '--ro' ,'LineWidth',1.1,'MarkerSize',3.0);...
    loglog(diameter_bin04,dmdlogdp_bin04_wat,  '--b'  ,'LineWidth',2.1,'MarkerSize',3.0);
title('Time 2.7s','FontSize',12)
xlabel('D_p (nm)','FontSize',12)
%ylabel('dM/dlogDp (ng/m3)','FontSize',12)
set(h(4),'xtick',[1,10,100,1000]);
set(h(4), 'xticklabel', [1,10,100,1000]);
set(h(4),'XLim',[1. 600.],'Ylim',[1.e1 1.e7]); 

%print -dbmp '../afigs_aging/diesel-sizedismass-spec.bmp'
%print -depsc '../afigs_aging/diesel-sizedismass-spec.eps'
print -djpg '../afigs_aging/diesel-sizedismass-spec.jpg'


figure(2)
axes('linewidth',2,'fontsize',16)
%modelled at 25m
loglog(diameter_bin01,dmdlogdp_bin01,      '-k' ,'LineWidth',2.4)
hold
ax=gca;
set(ax,'linewidth',2,'fontsize',13)
%modelled at xxx m
loglog(diameter_bin02,dmdlogdp_bin02,      '-r' ,'LineWidth',2.4)
loglog(diameter_bin03,dmdlogdp_bin03,      '-g' ,'LineWidth',2.4)
loglog(diameter_bin04,dmdlogdp_bin04,      '-b' ,'LineWidth',2.4)

legend('MAFOR 0.0s','MAFOR 0.1s','MAFOR 0.9s','MAFOR 2.7s','Location','NorthEastOutside')  

xlabel('D_p (nm)','FontSize',16)
ylabel('dM/dlogDp (ng/m3)','FontSize',16)
set(gca,'xtick',[1,10,100,1000]);
set(gca, 'xticklabel', [1,10,100,1000]);
set(gca,'XLim',[1. 600.],'Ylim',[1.e1 1.e7]); 

%print -dbmp '../afigs_aging/diesel-sizedismass.bmp'
%print -depsc '../afigs_aging/diesel-sizedismass.eps'
print -djpg '../afigs_aging/diesel-sizedismass.jpg'
