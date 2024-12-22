clear

load ../output/iconw1/wetdp.res;
load ../output/iconw1/size_dis.res;

%FIT dN/dlogDp

% SMPS at 2.7s (end of experiment)
bga16_0_nm=[9.37	9.93	10.5	11.2	11.8	12.5	13.3	14.1	14.9	15.8	16.8	17.8	18.8	19.9	21.1	22.4	23.7	25.2	26.7	28.3	29.9	31.7	33.6	35.6	37.8	40	42.4	45	47.7	50.5	53.5	56.7	60.1	63.7	67.5	71.6	75.8	80.4	85.2	90.3	95.7	101	107	114	121	128	136	144	152	161	171	181	192	204	216	229	242	257	272	288	306	324	343	364	386	409	433	459	487	516 ];

bga16_0_dndlogdp=[1.90E+00	7.84E+00	5.30E+00	2.65E+00	4.36E+00	3.37E+01	1.96E+02	6.56E+02	1.57E+03	2.83E+03	4.31E+03	5.47E+03	6.74E+03	7.77E+03	8.27E+03	8.73E+03	9.19E+03	9.72E+03	1.02E+04	1.04E+04	1.06E+04	1.07E+04	1.07E+04	1.07E+04	1.06E+04	1.05E+04	1.03E+04	1.01E+04	9.87E+03	9.53E+03	9.11E+03	8.60E+03	8.04E+03	7.44E+03	6.82E+03	6.22E+03	5.67E+03	5.17E+03	4.73E+03	4.34E+03	4.00E+03	3.11E+03	2.88E+03	2.67E+03	2.46E+03	2.27E+03	2.08E+03	1.90E+03	1.72E+03	1.55E+03	1.37E+03	1.18E+03	1.01E+03	8.52E+02	6.99E+02	5.60E+02	4.40E+02	3.42E+02	2.66E+02	2.06E+02	1.49E+02	9.92E+01	6.56E+01	4.43E+01	3.12E+01	2.23E+01	1.65E+01	1.22E+01	9.08E+00	7.93E+00 ];

% SMPS meas dN/dlogDp [# cm-3]
bg_number_dmps=bga16_0_dndlogdp(1:69);
bg_dndlogdp_dmps=bg_number_dmps;
bg_diameter_dmps(1:69)=0.;

for i=1:69
  bg_diameter_dmps(i)=bga16_0_nm(i);
end

a = size(bg_number_dmps)
b = size(bg_diameter_dmps)

%new structure of size_dis.res
% 1st line: dry diameter
% 2nd line: dlogDp
% 3rd line: dNdlogDp backgr
% 4th line: dNdlogDp(t=0)
% first value is model_time
%%%
infile='size_dis.res';
in=strrep(infile,'.res','');
y=eval(in);
[row,col]=size(y);                    %row=xxx col=61

infiledp='wetdp.res';
indp=strrep(infiledp,'.res','');
dp=eval(indp);
[rowdp,coldp]=size(dp);

%dndlogdp_binbg=y(3,2:col)   *1e-6  *2.303;
% 0.0s
dndlogdp_bin01=y(3,2:col)   *1e-6 *2.303;      % #/m3-->#/cm3
diameter_bin01=y(1,2:col)*1e9;
diameter_bin01=dp(1,2:coldp)*1e9;
% 1800s
dndlogdp_bin02=y(4+180,2:col)   *1e-6 *2.303;      % #/m3-->#/cm3
diameter_bin02=y(1+180,2:col)*1e9;
diameter_bin02=dp(1+180,2:coldp)*1e9;
% 3600s
dndlogdp_bin03=y(4+360,2:col)   *1e-6 *2.303;      % #/m3-->#/cm3
diameter_bin03=y(1+360,2:col)*1e9;
diameter_bin03=dp(1+360,2:coldp)*1e9;
% 7200s
dndlogdp_bin04=y(4+720,2:col)   *1e-6 *2.303;      % #/m3-->#/cm3
diameter_bin04=y(1+720,2:col)*1e9;
diameter_bin04=dp(1+720,2:coldp)*1e9;
% 10800s
dndlogdp_bin05=y(4+1080,2:col)   *1e-6 *2.303;      % #/m3-->#/cm3
diameter_bin05=y(1+1080,2:col)*1e9;
diameter_bin05=dp(1+1080,2:coldp)*1e9;
% 14400s
dndlogdp_bin06=y(4+1440,2:col)   *1e-6 *2.303;      % #/m3-->#/cm3
diameter_bin06=y(1+1440,2:col)*1e9;
diameter_bin06=dp(1+1440,2:coldp)*1e9;
% 18000s
dndlogdp_bin07=y(4+1800,2:col)   *1e-6 *2.303;      % #/m3-->#/cm3
diameter_bin07=y(1+1800,2:col)*1e9;
diameter_bin07=dp(1+1800,2:coldp)*1e9;

fsize=11;
fsizel=8;

%loglog
% first 8h: hourly. then every 8h
%1 hour= 1*60*6 = 360
% TIME UTC HERE AND ABOVE
figure(1)
axes('linewidth',1.5,'fontsize',fsize)
%SMPS
loglog(bg_diameter_dmps,bg_number_dmps,'bd','MarkerSize',5.0,'LineWidth',1.1)
hold
ax=gca;
set(ax,'linewidth',1.5,'fontsize',fsize)
%MAFOR   0.0s
loglog(diameter_bin01,dndlogdp_bin01,      '-k' ,'LineWidth',1.6)

%modelled at xxx m
loglog(diameter_bin02,dndlogdp_bin02,      '--y' ,'LineWidth',1.3)
loglog(diameter_bin03,dndlogdp_bin03,      '--m' ,'LineWidth',1.3)
loglog(diameter_bin04,dndlogdp_bin04,      '--g' ,'LineWidth',1.3)
loglog(diameter_bin05,dndlogdp_bin05,      '--c' ,'LineWidth',1.3)
loglog(diameter_bin06,dndlogdp_bin06,      '--b' ,'LineWidth',1.3)
loglog(diameter_bin07,dndlogdp_bin07,      '--r' ,'LineWidth',1.3)

%LEGEND
g=legend('SMPS bg','MAFOR bg','MAFOR 30min','MAFOR 60min','MAFOR 120min','MAFOR 180min','MAFOR 240min','MAFOR 300min','Location','NorthEastOutside')
set(g,'FontSize',fsizel)
xlabel('Wet diameter D_p (nm)','FontSize',fsize)
ylabel('dN/dlogDp (particles/cm^3)','FontSize',fsize)
title('ICONW=1 | Particle size distr. (wet Dp)','FontSize',fsize,'FontName','Arial')
%loglog
set(gca,'xtick',[1,10,50,100,500,1000]);
set(gca, 'xticklabel', [1,10,50,100,500,1000]);
set(gca,'XLim',[1. 700.],'Ylim',[5.e1 1.e5]);

print -djpg '../afigs_mesa/cabauw_sizedis_iconw1.jpg'


load ../output/iconw2/wetdp.res;
load ../output/iconw2/size_dis.res;

infile='size_dis.res';
in=strrep(infile,'.res','');
y=eval(in);
[row,col]=size(y);                    %row=xxx col=61

infiledp='wetdp.res';
indp=strrep(infiledp,'.res','');
dp=eval(indp);
[rowdp,coldp]=size(dp);

%dndlogdp_binbg=y(3,2:col)   *1e-6  *2.303;
% 0.0s
dndlogdp_bin01=y(3,2:col)   *1e-6 *2.303;      % #/m3-->#/cm3
diameter_bin01=y(1,2:col)*1e9;
diameter_bin01=dp(1,2:coldp)*1e9;
% 1800s
dndlogdp_bin02=y(4+180,2:col)   *1e-6 *2.303;      % #/m3-->#/cm3
diameter_bin02=y(1+180,2:col)*1e9;
diameter_bin02=dp(1+180,2:coldp)*1e9;
% 3600s
dndlogdp_bin03=y(4+360,2:col)   *1e-6 *2.303;      % #/m3-->#/cm3
diameter_bin03=y(1+360,2:col)*1e9;
diameter_bin03=dp(1+360,2:coldp)*1e9;
% 7200s
dndlogdp_bin04=y(4+720,2:col)   *1e-6 *2.303;      % #/m3-->#/cm3
diameter_bin04=y(1+720,2:col)*1e9;
diameter_bin04=dp(1+720,2:coldp)*1e9;
% 10800s
dndlogdp_bin05=y(4+1080,2:col)   *1e-6 *2.303;      % #/m3-->#/cm3
diameter_bin05=y(1+1080,2:col)*1e9;
diameter_bin05=dp(1+1080,2:coldp)*1e9;
% 14400s
dndlogdp_bin06=y(4+1440,2:col)   *1e-6 *2.303;      % #/m3-->#/cm3
diameter_bin06=y(1+1440,2:col)*1e9;
diameter_bin06=dp(1+1440,2:coldp)*1e9;
% 18000s
dndlogdp_bin07=y(4+1800,2:col)   *1e-6 *2.303;      % #/m3-->#/cm3
diameter_bin07=y(1+1800,2:col)*1e9;
diameter_bin07=dp(1+1800,2:coldp)*1e9;


%loglog
% first 8h: hourly. then every 8h
%1 hour= 1*60*6 = 360
% TIME UTC HERE AND ABOVE
figure(2)
axes('linewidth',1.5,'fontsize',fsize)
%SMPS
loglog(bg_diameter_dmps,bg_number_dmps,'bd','MarkerSize',5.0,'LineWidth',1.1)
hold
ax=gca;
set(ax,'linewidth',1.5,'fontsize',fsize)
%MAFOR   0.0s
loglog(diameter_bin01,dndlogdp_bin01,      '-k' ,'LineWidth',1.6)

%modelled at xxx m
loglog(diameter_bin02,dndlogdp_bin02,      '--y' ,'LineWidth',1.3)
loglog(diameter_bin03,dndlogdp_bin03,      '--m' ,'LineWidth',1.3)
loglog(diameter_bin04,dndlogdp_bin04,      '--g' ,'LineWidth',1.3)
loglog(diameter_bin05,dndlogdp_bin05,      '--c' ,'LineWidth',1.3)
loglog(diameter_bin06,dndlogdp_bin06,      '--b' ,'LineWidth',1.3)
loglog(diameter_bin07,dndlogdp_bin07,      '--r' ,'LineWidth',1.3)

%LEGEND
g=legend('SMPS bg','MAFOR bg','MAFOR 30min','MAFOR 60min','MAFOR 120min','MAFOR 180min','MAFOR 240min','MAFOR 300min','Location','NorthEastOutside')
set(g,'FontSize',fsizel)
xlabel('Wet diameter D_p (nm)','FontSize',fsize)
ylabel('dN/dlogDp (particles/cm^3)','FontSize',fsize)
title('ICONW=2 | Particle size distr. (wet Dp)','FontSize',fsize,'FontName','Arial')
%loglog
set(gca,'xtick',[1,10,50,100,500,1000]);
set(gca, 'xticklabel', [1,10,50,100,500,1000]);
set(gca,'XLim',[1. 700.],'Ylim',[5.e1 1.e5]);

print -djpg '../afigs_mesa/cabauw_sizedis_iconw2.jpg'

