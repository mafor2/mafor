clear

load ../output/iconw1/aerconc.res;
load ../output/iconw1/total_n.res;

time=aerconc(:,1);
time=time/3600.;
tmax=6*60*7;
%time column
infile='aerconc.res';
in=strrep(infile,'.res','');
y=eval(in);
[row,col]=size(y);
t=y(3:row,1);
t=t/3600.;
%jd=234;
%t=t/24;
%time=t+jd;
time=t-9;

%AEROSOL MASS CONCENTRATION
pmsulct=aerconc(:,2)+aerconc(:,3)+aerconc(:,4)+aerconc(:,5);
pmmsact=aerconc(:,6)+aerconc(:,7)+aerconc(:,8)+aerconc(:,9);
%no msa_p
%xxxx=seasalt (no BPBA)
pmxxxct=aerconc(:,10)+aerconc(:,11)+aerconc(:,12)+aerconc(:,13);
pmorgct=aerconc(:,14)+aerconc(:,15)+aerconc(:,16)+aerconc(:,17);
pmamoct=aerconc(:,18)+aerconc(:,19)+aerconc(:,20)+aerconc(:,21);
pmnitct=aerconc(:,22)+aerconc(:,23)+aerconc(:,24)+aerconc(:,25);
pmecbct=aerconc(:,26)+aerconc(:,27)+aerconc(:,28)+aerconc(:,29);
%pmdusct=aerconc(:,30)+aerconc(:,31)+aerconc(:,32)+aerconc(:,33);
pmseact=aerconc(:,34)+aerconc(:,35)+aerconc(:,36)+aerconc(:,37);
pmwatct=aerconc(:,38)+aerconc(:,39)+aerconc(:,40)+aerconc(:,41);

pmsulf =pmsulct*1.e-3;
pmmsap =pmmsact*1.e-3;
pmammo =pmamoct*1.e-3;
pmanit =pmnitct*1.e-3;
pmwatt =pmwatct*1.e-3;
pmdust =pmxxxct*1.e-3;
pmorgt =pmorgct*1.e-3;
pmecbt =pmecbct*1.e-3; 
pmseat =pmseact*1.e-3;

% TOTAL in ug/m3 (dry mass)
pmtotal=(pmorgct+pmecbct+pmsulct+pmxxxct+pmamoct+pmnitct+pmseact)*1e-3;

% TOTAL number conc in #/cm3
pnctotal =1.e-6*(total_n(:,2) +total_n(:,4)+total_n(:,5)+total_n(:,6));  % >3nm
% 10-25nm
pncnucl =1.e-6*(total_n(:,7)-total_n(:,3));
% 25-100nm
pncaitk =1.e-6*total_n(:,8);
% >100nm
pncaccu =1.e-6*total_n(:,6);


figure(1)
plot(time(1:tmax),pmtotal(1:tmax),'-b','LineWidth',2.1)
hold
ax=gca;
set(ax,'linewidth',1.5,'fontsize',14)
plot(time(1:tmax),pmorgt(1:tmax),'--m','LineWidth',2.1)
plot(time(1:tmax),pmsulf(1:tmax),'--y','LineWidth',2.1)
plot(time(1:tmax),pmanit(1:tmax),'--g','LineWidth',2.1)
%plot(time(1:tmax),pmmsap(1:tmax),'--c','LineWidth',2.1)
plot(time(1:tmax),pmseat(1:tmax),'--c','LineWidth',2.1)
plot(time(1:tmax),pmammo(1:tmax),'--r','LineWidth',2.1)
xlabel('Time (hours)','FontSize',14,'FontName','Arial')
ylabel('Particle mass conc. (\mug m^{-3})','FontSize',14,'FontName','Arial')
title('ICONW=1 | Particle mass conc.','FontSize',17,'FontName','Arial')
legend('total dry','OC','SO4','NO3','NaCl','NH4','Location','NorthWest')

axis([ -0.1, 7,   0.0, 20.0]) 
JD=0:1:8;
set(gca,'xtick',JD);
JDtext=['0';'1';'2';'3';'4';'5';'6';'7';'8']; 
set(gca, 'xticklabel', JDtext);
grid on

print -djpg '../afigs_mesa/cabauw_partmass_iconw1.jpg'


load ../output/iconw2/aerconc.res;
load ../output/iconw2/total_n.res;

time=aerconc(:,1);
time=time/3600.;
tmax=6*60*7;
%time column
infile='aerconc.res';
in=strrep(infile,'.res','');
y=eval(in);
[row,col]=size(y);
t=y(3:row,1);
t=t/3600.;
%jd=234;
%t=t/24;
%time=t+jd;
time=t-9;

%AEROSOL MASS CONCENTRATION
pmsulct=aerconc(:,2)+aerconc(:,3)+aerconc(:,4)+aerconc(:,5);
pmmsact=aerconc(:,6)+aerconc(:,7)+aerconc(:,8)+aerconc(:,9);
%no msa_p
%xxxx=seasalt (no BPBA)
pmxxxct=aerconc(:,10)+aerconc(:,11)+aerconc(:,12)+aerconc(:,13);
pmorgct=aerconc(:,14)+aerconc(:,15)+aerconc(:,16)+aerconc(:,17);
pmamoct=aerconc(:,18)+aerconc(:,19)+aerconc(:,20)+aerconc(:,21);
pmnitct=aerconc(:,22)+aerconc(:,23)+aerconc(:,24)+aerconc(:,25);
pmecbct=aerconc(:,26)+aerconc(:,27)+aerconc(:,28)+aerconc(:,29);
%pmdusct=aerconc(:,30)+aerconc(:,31)+aerconc(:,32)+aerconc(:,33);
pmseact=aerconc(:,34)+aerconc(:,35)+aerconc(:,36)+aerconc(:,37);
pmwatct=aerconc(:,38)+aerconc(:,39)+aerconc(:,40)+aerconc(:,41);
pmsulf =pmsulct*1.e-3;
pmmsap =pmmsact*1.e-3;
pmammo =pmamoct*1.e-3;
pmanit =pmnitct*1.e-3;
pmwatt =pmwatct*1.e-3;
pmdust =pmxxxct*1.e-3;
pmorgt =pmorgct*1.e-3;
pmecbt =pmecbct*1.e-3; 
pmseat =pmseact*1.e-3;

% TOTAL in ug/m3 (dry mass)
pmtotal=(pmorgct+pmecbct+pmsulct+pmxxxct+pmamoct+pmnitct+pmseact)*1e-3;

% TOTAL number conc in #/cm3
pnctotal =1.e-6*(total_n(:,2) +total_n(:,4)+total_n(:,5)+total_n(:,6));  % >3nm
% 10-25nm
pncnucl =1.e-6*(total_n(:,7)-total_n(:,3));
% 25-100nm
pncaitk =1.e-6*total_n(:,8);
% >100nm
pncaccu =1.e-6*total_n(:,6);


figure(2)
plot(time(1:tmax),pmtotal(1:tmax),'-b','LineWidth',2.1)
hold
ax=gca;
set(ax,'linewidth',1.5,'fontsize',15)
plot(time(1:tmax),pmorgt(1:tmax),'--m','LineWidth',2.1)
plot(time(1:tmax),pmsulf(1:tmax),'--y','LineWidth',2.1)
plot(time(1:tmax),pmanit(1:tmax),'--g','LineWidth',2.1)
%plot(time(1:tmax),pmmsap(1:tmax),'--c','LineWidth',2.1)
plot(time(1:tmax),pmseat(1:tmax),'--c','LineWidth',2.1)
plot(time(1:tmax),pmammo(1:tmax),'--r','LineWidth',2.1)
xlabel('Time from start (h)','FontSize',14,'FontName','Arial')
ylabel('Particle mass conc. (\mug m^{-3})','FontSize',14,'FontName','Arial')
title('ICONW=2 | Particle mass conc.','FontSize',17,'FontName','Arial')
legend('total dry','OC','SO4','NO3','NaCl','NH4','Location','NorthWest')

axis([ -0.1, 7,   0.0, 20.0]) 
JD=0:1:8;
set(gca,'xtick',JD);
JDtext=['0';'1';'2';'3';'4';'5';'6';'7';'8']; 
set(gca, 'xticklabel', JDtext);
grid on

print -djpg '../afigs_mesa/cabauw_partmass_iconw2.jpg'
