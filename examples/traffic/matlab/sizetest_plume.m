clear


load ../output/wetdp.res
load ../output/size_dis.res;
load ../output/size_dism.res;

%SMPS dN/dlogDp
%MotorwayA16_MKA.xls
% At 25m downwind
ma16_0_nm=[10.9	11.3	11.8	12.2	12.6	13.1	13.6	14.1	14.6	15.1	15.7	16.3	16.8	17.5	18.1	18.8	19.5	20.2	20.9	21.7	22.5	23.3	24.1	25	25.9	26.9	27.9	28.9	30	31.1	32.2	33.4	34.6	35.9	37.2	38.5	40	41.4	42.9	44.5	46.1	47.8	49.6	51.4	53.3	55.2	57.3	59.4	61.5	63.8	66.1	68.5	71	73.7	76.4	79.1	82	85.1	88.2	91.4	94.7	98.2	101.8	105.5	109.4	113.4	117.6	121.9	126.3	131	135.8	140.7	145.9	151.2	156.8	162.5	168.5	174.7	181.1	187.7	194.6	201.7	209.1	216.7	224.7	232.9	241.4	250.3	259.5	269	278.8	289	299.6	310.6	322	333.8	346	358.7	371.8	385.4	399.5	414.2	429.4	445.1	461.4	478.3	495.8];
ma16_0_dndlogdp=[1.8E+04	3.8E+04	2.7E+04	2.2E+04	2.3E+04	2.5E+04	4.2E+04	3.0E+04	3.0E+04	4.0E+04	3.6E+04	4.0E+04	5.1E+04	5.8E+04	6.2E+04	6.0E+04	5.6E+04	5.3E+04	4.8E+04	4.6E+04	5.6E+04	4.8E+04	4.3E+04	4.3E+04	4.8E+04	5.1E+04	5.2E+04	6.2E+04	7.5E+04	8.0E+04	7.3E+04	6.5E+04	5.5E+04	4.8E+04	4.3E+04	3.8E+04	3.4E+04	3.3E+04	3.5E+04	3.2E+04	3.3E+04	3.2E+04	3.3E+04	3.3E+04	3.3E+04	3.5E+04	3.3E+04	3.3E+04	3.3E+04	3.2E+04	3.2E+04	3.2E+04	3.1E+04	3.1E+04	3.0E+04	2.9E+04	2.7E+04	2.6E+04	2.6E+04	2.4E+04	2.3E+04	2.2E+04	1.9E+04	1.9E+04	1.7E+04	1.6E+04	1.5E+04	1.6E+04	1.4E+04	1.2E+04	1.1E+04	1.1E+04	9.5E+03	8.2E+03	7.6E+03	7.3E+03	6.6E+03	6.1E+03	5.8E+03	5.4E+03	5.5E+03	5.2E+03	4.5E+03	4.1E+03	3.6E+03	3.0E+03	2.9E+03	2.5E+03	2.3E+03	2.1E+03	1.9E+03	3.1E+03	1.4E+03	1.2E+03	1.2E+03	9.3E+02	9.0E+02	7.4E+02	7.9E+02	1.5E+03	8.7E+02	6.4E+02	4.1E+02	3.6E+02	3.6E+02	2.1E+02	2.3E+02];
ma16_1_dndlogdp=[2.7E+04	2.9E+04	2.8E+04	2.6E+04	3.0E+04	3.6E+04	3.0E+04	3.0E+04	3.0E+04	2.9E+04	2.9E+04	3.5E+04	3.4E+04	3.0E+04	3.3E+04	2.8E+04	2.6E+04	2.5E+04	2.7E+04	2.5E+04	2.5E+04	2.6E+04	2.8E+04	2.7E+04	2.7E+04	2.9E+04	2.9E+04	2.9E+04	3.1E+04	2.8E+04	2.8E+04	2.8E+04	2.7E+04	2.7E+04	2.8E+04	2.7E+04	2.7E+04	2.7E+04	2.9E+04	2.8E+04	2.9E+04	2.9E+04	3.0E+04	2.9E+04	3.1E+04	3.1E+04	3.0E+04	3.1E+04	3.2E+04	3.1E+04	3.0E+04	3.0E+04	2.9E+04	2.9E+04	2.8E+04	2.7E+04	2.6E+04	2.5E+04	2.4E+04	2.3E+04	2.2E+04	2.1E+04	1.9E+04	1.7E+04	1.7E+04	1.5E+04	1.4E+04	1.2E+04	1.1E+04	1.1E+04	9.0E+03	7.8E+03	6.9E+03	6.5E+03	7.5E+03	6.0E+03	5.5E+03	4.9E+03	4.4E+03	3.7E+03	4.0E+03	3.3E+03	3.1E+03	2.9E+03	3.7E+03	2.6E+03	2.4E+03	2.6E+03	2.1E+03	1.8E+03	1.6E+03	1.5E+03	2.2E+03	2.0E+03	1.3E+03	1.7E+03	1.1E+03	1.0E+03	8.6E+02	2.5E+03	3.2E+03	1.3E+03	7.8E+02	9.1E+02	9.7E+02	1.3E+03	4.5E+02];
ma16_2_dndlogdp=[1.8E+04	2.6E+04	3.8E+04	2.7E+04	3.6E+04	2.5E+04	3.3E+04	2.7E+04	2.5E+04	2.1E+04	2.2E+04	3.8E+04	2.7E+04	2.4E+04	2.6E+04	2.8E+04	3.1E+04	3.4E+04	3.5E+04	3.7E+04	3.6E+04	3.7E+04	3.4E+04	3.1E+04	2.9E+04	3.0E+04	3.1E+04	3.1E+04	2.9E+04	2.8E+04	2.7E+04	2.7E+04	3.0E+04	2.7E+04	2.7E+04	2.8E+04	2.7E+04	2.8E+04	2.7E+04	2.8E+04	2.9E+04	2.9E+04	2.9E+04	2.8E+04	2.7E+04	2.9E+04	2.8E+04	2.9E+04	3.1E+04	2.9E+04	2.7E+04	2.7E+04	2.6E+04	2.5E+04	2.5E+04	2.4E+04	2.3E+04	2.3E+04	2.2E+04	2.1E+04	2.1E+04	2.0E+04	1.8E+04	1.7E+04	1.6E+04	1.4E+04	1.4E+04	1.3E+04	1.1E+04	9.7E+03	8.3E+03	7.6E+03	6.8E+03	6.3E+03	5.5E+03	5.2E+03	4.7E+03	3.6E+03	4.8E+03	3.7E+03	2.8E+03	3.6E+03	3.0E+03	3.6E+03	4.3E+03	3.1E+03	2.2E+03	1.9E+03	1.9E+03	2.4E+03	2.9E+03	2.2E+03	1.7E+03	3.0E+03	2.1E+03	1.6E+03	1.3E+03	8.2E+02	8.4E+02	2.9E+03	2.2E+03	1.5E+03	7.5E+02	4.9E+02	3.3E+02	3.0E+03	1.4E+03];
ma16_3_dndlogdp=[2.5E+04	3.1E+04	2.2E+04	3.0E+04	2.0E+04	1.9E+04	1.9E+04	2.0E+04	2.0E+04	2.6E+04	2.2E+04	1.8E+04	1.8E+04	1.8E+04	1.8E+04	2.2E+04	2.2E+04	2.0E+04	1.8E+04	2.7E+04	2.2E+04	2.1E+04	1.9E+04	2.0E+04	2.1E+04	2.0E+04	2.0E+04	2.0E+04	2.1E+04	2.1E+04	2.0E+04	2.0E+04	2.2E+04	2.2E+04	2.2E+04	2.2E+04	2.3E+04	2.3E+04	2.5E+04	2.5E+04	2.6E+04	2.7E+04	2.8E+04	2.8E+04	2.9E+04	3.0E+04	3.0E+04	3.1E+04	3.0E+04	3.0E+04	3.0E+04	3.0E+04	3.0E+04	2.9E+04	2.9E+04	2.9E+04	2.7E+04	2.9E+04	2.7E+04	2.4E+04	2.3E+04	2.2E+04	2.0E+04	2.0E+04	1.7E+04	1.7E+04	1.5E+04	1.3E+04	1.2E+04	1.1E+04	9.5E+03	8.5E+03	7.4E+03	6.6E+03	5.8E+03	5.5E+03	5.0E+03	4.9E+03	4.9E+03	4.5E+03	3.9E+03	3.6E+03	3.3E+03	3.7E+03	3.5E+03	3.3E+03	3.6E+03	3.1E+03	2.6E+03	2.4E+03	2.2E+03	1.9E+03	2.2E+03	1.8E+03	1.8E+03	1.6E+03	1.6E+03	1.6E+03	9.6E+02	6.7E+02	6.2E+02	3.8E+02	3.0E+02	2.7E+02	1.9E+02	1.6E+02	1.3E+02];

% Background at 25m
bga16_0_nm=[10.4	11.1	12	12.9	13.8	14.9	16	17.2	18.4	19.8	21.3	22.9	24.6	26.4	28.4	30.5	32.8	35.2	37.9	40.7	43.7	47	50.5	54.2	58.3	62.6	67.3	72.3	77.7	83.5	89.8	96.5	103.7	111.4	119.7	128.6	138.2	148.6	159.6	171.5	184.3	198.1	212.9	228.8	245.8	264.2	283.9	305.1	327.8	352.3	378.6	406.8	437.1	469.8];
bga16_0_dndlogdp=[1.7E+04	2.1E+04	2.3E+04	2.3E+04	2.1E+04	2.0E+04	1.9E+04	2.0E+04	2.6E+04	2.5E+04	2.3E+04	2.1E+04	2.1E+04	1.9E+04	1.7E+04	1.7E+04	1.7E+04	1.8E+04	1.8E+04	1.9E+04	2.0E+04	2.1E+04	2.2E+04	2.2E+04	2.2E+04	2.1E+04	2.1E+04	2.1E+04	2.0E+04	1.9E+04	1.7E+04	1.5E+04	1.2E+04	1.0E+04	9.0E+03	7.4E+03	6.0E+03	5.1E+03	4.2E+03	3.5E+03	2.8E+03	2.4E+03	2.0E+03	1.7E+03	1.4E+03	1.1E+03	8.3E+02	6.6E+02	4.5E+02	3.3E+02	2.5E+02	1.8E+02	1.4E+02	1.1E+02];

% SMPS meas dN/dlogDp [# cm-3]
bg_number_dmps=bga16_0_dndlogdp(1:54);
bg_dndlogdp_dmps=bg_number_dmps;
bg_diameter_dmps(1:54)=0.;
ma0_number_dmps=ma16_0_dndlogdp(1:107);
ma0_dndlogdp_dmps=ma0_number_dmps;
ma0_diameter_dmps(1:107)=0.;
ma1_number_dmps=ma16_1_dndlogdp(1:107);
ma1_dndlogdp_dmps=ma1_number_dmps;
ma2_number_dmps=ma16_2_dndlogdp(1:107);
ma2_dndlogdp_dmps=ma2_number_dmps;
ma3_number_dmps=ma16_3_dndlogdp(1:107);
ma3_dndlogdp_dmps=ma3_number_dmps;

for i=1:54
  bg_diameter_dmps(i)=bga16_0_nm(i);
end
% motorway A16 diamter_smps %[nm]
for i=1:107
  ma0_diameter_dmps(i)=ma16_0_nm(i);
end

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

diameter_bin=y(1,2:col)*1e9;          % m-->nm
diameter_bin=dp(4,2:coldp)*1e9;

% multiply by 2.303 to convert ln to log10
dndlogdp_binbg=y(3,2:col)   *1e-6 *2.303;
% 25m distance (t0=25m, u=5m/s)
dndlogdp_bin01=y(4,2:col)   *1e-6 *2.303;    % #/m3-->#/cm3
% 500m distance
dndlogdp_bin02=y(4+10,2:col) *1e-6 *2.303;
% 1000m distance
dndlogdp_bin03=y(4+20,2:col) *1e-6 *2.303;
% 2000m distance
dndlogdp_bin04=y(4+40,2:col) *1e-6 *2.303;
% 3000m distance
dndlogdp_bin05=y(4+60,2:col) *1e-6 *2.303;
% 4000m distance
dndlogdp_bin06=y(4+80,2:col) *1e-6 *2.303;
% 5000m distance
%dndlogdp_bin06=y(4+100,2:col) *1e-6 *2.303;
% 6000m distance
dndlogdp_bin07=y(4+120,2:col) *1e-6 *2.303;

% Mass distribution
% multiply by 2.303 to convert ln to log10
infileam='size_dism.res';
inam=strrep(infileam,'.res','');
yam=eval(inam);
[row,col]=size(yam);                    %row=xxx col=61
% 25m distance (t0=25m, u=5m/s)
dmdlogdp_bin01=yam(2,2:col)   *1e12 *2.303;    % kg/m3-->ng/m3
% 500m distance
dmdlogdp_bin02=yam(2+10,2:col) *1e12 *2.303;
% 1000m distance
dmdlogdp_bin03=yam(2+20,2:col) *1e12 *2.303;
% 2000m distance
dmdlogdp_bin04=yam(2+40,2:col) *1e12 *2.303;
% 3000m distance
dmdlogdp_bin05=yam(2+60,2:col) *1e12 *2.303;
% 4000m distance
dmdlogdp_bin06=yam(2+80,2:col) *1e12 *2.303;
% 5000m distance
%dmdlogdp_bin06=yam(2+100,2:col) *1e12 *2.303;
% 6000m distance
dmdlogdp_bin07=yam(2+120,2:col) *1e12 *2.303;

fsize=11;

% first 8h: hourly. then every 8h
%1 hour= 1*60*6 = 360
% TIME UTC HERE AND ABOVE
figure(1)
axes('linewidth',1.5,'fontsize',fsize)
%SMPS
loglog(ma0_diameter_dmps,ma16_0_dndlogdp,'ks-','MarkerSize',3.5,'LineWidth',0.9)
hold
ax=gca;
set(ax,'linewidth',1.5,'fontsize',fsize,'tickdir','out')
%loglog(ma0_diameter_dmps,ma16_2_dndlogdp,'rd-','MarkerSize',3.5,'LineWidth',0.9)
loglog(bg_diameter_dmps,bga16_0_dndlogdp,'kd-','MarkerSize',3.5,'LineWidth',0.9)
%modelled background size distribution
%loglog(diameter_bin,dndlogdp_binbg,      '-k' ,'LineWidth',2.1)
%modelled at 25m
loglog(diameter_bin,dndlogdp_bin01,'Color',[.6 .6 .6],'LineStyle','-','LineWidth',2.5)
%modelled at xxx m
loglog(diameter_bin,dndlogdp_bin02,      '--r' ,'LineWidth',1.8)
loglog(diameter_bin,dndlogdp_bin03,      '--g' ,'LineWidth',1.8)
loglog(diameter_bin,dndlogdp_bin04,      '--y' ,'LineWidth',1.8)
loglog(diameter_bin,dndlogdp_bin05,      '--c' ,'LineWidth',1.8)
loglog(diameter_bin,dndlogdp_bin06,      '--m' ,'LineWidth',1.8)
loglog(diameter_bin,dndlogdp_bin07,      '--b' ,'LineWidth',1.8)

grid

legend('A16dw,25m','backgr,25m','25m model','500m model','1000m model','2000m model','3000m model','4000m model','6000m model','Location','NorthEastOutside')
xlabel('Diameter D_p (nm)','FontSize',fsize)
ylabel('dN/dlog_{10}Dp (particles/cm^3)','FontSize',fsize)
%loglog
set(gca,'xtick',[1,10,100,1000]);
set(gca, 'xticklabel', [1,10,100,1000]);
set(gca,'XLim',[5. 1000.],'Ylim',[1.e1 1.e6]);
%uncomment below line to save as jpg
print -djpg '../afigs_plume/ma16-sizedisnum-all-dist.jpg'


% first 8h: hourly. then every 8h
%1 hour= 1*60*6 = 360
% TIME UTC HERE AND ABOVE
figure(2)
axes('linewidth',1.5,'fontsize',fsize)
%modelled at 25m
loglog(diameter_bin,dmdlogdp_bin01,'Color',[.6 .6 .6],'LineStyle','-','LineWidth',2.5)
hold
ax=gca;
set(ax,'linewidth',1.5,'fontsize',fsize,'tickdir','out')
%modelled at xxx m
loglog(diameter_bin,dmdlogdp_bin02,      '--r' ,'LineWidth',1.8)
loglog(diameter_bin,dmdlogdp_bin03,      '--g' ,'LineWidth',1.8)
loglog(diameter_bin,dmdlogdp_bin04,      '--y' ,'LineWidth',1.8)
loglog(diameter_bin,dmdlogdp_bin05,      '--c' ,'LineWidth',1.8)
loglog(diameter_bin,dmdlogdp_bin06,      '--m' ,'LineWidth',1.8)
loglog(diameter_bin,dmdlogdp_bin07,      '--b' ,'LineWidth',1.8)
grid
legend('25m model','500m model','1000m model','2000m model','3000m model','4000m model','6000m model','Location','NorthWest')
xlabel('Diameter D_p (nm)','FontSize',fsize)
ylabel('dM/dlog_{10}Dp (ng/m3)','FontSize',fsize)
set(gca,'xtick',[1,10,100,1000]);
set(gca, 'xticklabel', [1,10,100,1000]);
set(gca,'XLim',[5. 1000.],'Ylim',[1.e1 1.e5]);

%uncomment below line to save as jpg
print -djpg '../afigs_plume/ma16-sizedismas-all-dist.jpg'
