clear
load ../output/size_dis.res;
load ../output/wetdp.res;

%DMPS dN/dlogDp
load ../observ/stat3_all.txt

% DMPS meas N [cm-3]
%number_dmps=stat3_all(2,2:31);
number_dmps_00=stat3_all(17,2:31);     %209.406
number_dmps_10=stat3_all(31,2:31);     %209.504
number_dmps_20=stat3_all(45,2:31);     %209.601
number_dmps_70=stat3_all(70,2:31);     %209.802
dlogdp_dmps=0.078;
dndlogdp_dmps_00=number_dmps_00/dlogdp_dmps;
dndlogdp_dmps_10=number_dmps_10/dlogdp_dmps;
dndlogdp_dmps_20=number_dmps_20/dlogdp_dmps;
dndlogdp_dmps_70=number_dmps_70/dlogdp_dmps;

%%diameter_dmps=stat3_all(1,2:31);   %[nm]
diameter_dmps( 1)=2.8;
diameter_dmps( 2)=3.7 ;
diameter_dmps( 3)=4.3;
diameter_dmps( 4)=5.094;
diameter_dmps( 5)=6.089;
diameter_dmps( 6)=7.279;
diameter_dmps( 7)=8.701;
diameter_dmps( 8)=10.401;
diameter_dmps( 9)=12.433;
diameter_dmps(10)=14.862;
diameter_dmps(11)=17.766;
diameter_dmps(12)=21.237;
diameter_dmps(13)=29.122;
diameter_dmps(14)=34.812;
diameter_dmps(15)=41.613;
diameter_dmps(16)=49.743;
diameter_dmps(17)=59.872;
diameter_dmps(18)=71.569;
diameter_dmps(19)=85.552;
diameter_dmps(20)=102.27;
diameter_dmps(21)=122.25;
diameter_dmps(22)=146.13;
diameter_dmps(23)=174.68;
diameter_dmps(24)=208.81;
diameter_dmps(25)=249.61;
diameter_dmps(26)=300.44;
diameter_dmps(27)=359.13;
diameter_dmps(28)=429.3;
diameter_dmps(29)=513.18;
diameter_dmps(30)=613.44;

%new structure of size_dis.res
% 1st line: dry diameter
% 2nd line: dlogDp
% 3rd line: dNdlogDp(t=0)
% 4th line: wet diameter
% first value is model_time
%%%
infile='size_dis.res';
in=strrep(infile,'.res','');
y=eval(in);
[row,col]=size(y);

infiledp='wetdp.res';
indp=strrep(infiledp,'.res','');
dp=eval(indp);
[rowdp,coldp]=size(dp);
diameter_bin=dp(3,2:coldp)*1e9;     % m-->nm

% multiply by 2.303 to convert ln to log10
dndlogdp_bin01=y(5,2:col)    *1.e-6 *2.303;
dndlogdp_bin02=y(854,2:col)  *1.e-6 *2.303;
dndlogdp_bin03=y(1684,2:col) *1.e-6 *2.303;
dndlogdp_bin04=y(3424,2:col) *1.e-6 *2.303;

diameter_bin02=dp(852,2:coldp)*1e9;

imax=15;
nmodes=4;
bmax=nmodes*imax;
%1 hour= 1*60*6 = 360

figure(1)
semilogx(diameter_dmps(1:30),dndlogdp_dmps_00(1:30),'ks','LineWidth',1.5)
hold
ax=gca;
set(ax,'linewidth',1.5,'fontsize',15.5)
semilogx(diameter_bin,dndlogdp_bin01,      '-k' ,'LineWidth',2.5)
%209.504  (= 2hr 20min)
semilogx(diameter_dmps(1:30),dndlogdp_dmps_10(1:30),'gs','LineWidth',1.5)
semilogx(diameter_bin,dndlogdp_bin02,      '-g' ,'LineWidth',2.5)
%209.601  (= 4hr 40min)
semilogx(diameter_dmps(1:30),dndlogdp_dmps_20(1:30),'bs','LineWidth',1.5)
semilogx(diameter_bin,dndlogdp_bin03,      '-b' ,'LineWidth',2.5)
%209.802  (= 9hr 30min)
semilogx(diameter_dmps(1:30),dndlogdp_dmps_70(1:30),'rs','LineWidth',1.5)
semilogx(diameter_bin,dndlogdp_bin04,      '-r' ,'LineWidth',2.5)

legend('209.406 DMPS','209.406 model','209.504 DMPS','209.504 model','209.601 DMPS','209.601 model','209.802 DMPS','209.802 model','Location','NorthEast')
xlabel('D_p (nm)','FontSize',20,'FontName','Arial')
ylabel('dN/dlog_{10}Dp (particles/cm^{3})','FontSize',20,'FontName','Arial')

set(gca,'xtick',[1,10,100,1000]);
set(gca, 'xticklabel', [1,10,100,1000]);
set(gca,'XLim',[1. 1000.],'Ylim',[0 2800]); 

%uncomment below line to save as jpg
print -djpg '../afigs_cov/sizedis_209.jpg'

