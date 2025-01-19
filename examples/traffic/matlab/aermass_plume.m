clear

load ../output/aerconc.res;
load ../output/total_n.res;
load ../output/plume.res;

bgntot=21362.1+28367.7+22087.5+1167.5;

time=aerconc(:,1);
time=time/3600.;
tmax=6*60*1;

for i=1:tmax
 pncbg(i)=bgntot;
end

%EC from MAAP
ecmaap=[1. 0.749047493 0.608520109 0.535246452];
tmaap=[39600 39610 39620 39630];
tmaap=tmaap/3600.;

%PNC from Menno
pncsmps       =[1 0.726037573 0.64795775  0.630396336];
pncl25smps    =[1 0.69797206  0.628781587 0.523065688];
pnc25_100smps =[1 0.715098648 0.621513036 0.651665691];
pncgt100smps  =[1 0.896830392 0.873825402 0.920414649];


%AEROSOL MASS CONCENTRATION
pmsulct=aerconc(:,2)+aerconc(:,3)+aerconc(:,4)+aerconc(:,5)+aerconc(:,6);
%no msa_p
%no io3_p
pmxxxct=aerconc(:,17)+aerconc(:,18)+aerconc(:,19)+aerconc(:,20)+aerconc(:,21); % bpba
pmorgct=aerconc(:,22)+aerconc(:,23)+aerconc(:,24)+aerconc(:,25)+aerconc(:,26);
pmamoct=aerconc(:,27)+aerconc(:,28)+aerconc(:,29)+aerconc(:,30)+aerconc(:,31);
pmnitct=aerconc(:,32)+aerconc(:,33)+aerconc(:,34)+aerconc(:,35)+aerconc(:,36);
pmecbct=aerconc(:,37)+aerconc(:,38)+aerconc(:,39)+aerconc(:,40)+aerconc(:,41);
pmdusct=aerconc(:,42)+aerconc(:,43)+aerconc(:,44)+aerconc(:,45)+aerconc(:,46);
pmsalct=aerconc(:,47)+aerconc(:,48)+aerconc(:,49)+aerconc(:,50)+aerconc(:,51);
pmwatct=aerconc(:,52)+aerconc(:,53)+aerconc(:,54)+aerconc(:,55)+aerconc(:,56);
pmsulf =pmsulct*1.e-3;
pmbpba =pmxxxct*1.e-3;
pmorgt =pmorgct*1.e-3;
pmammo =pmamoct*1.e-3;
pmanit =pmnitct*1.e-3;
pmecbt =pmecbct*1.e-3;
pmdust =pmdusct*1.e-3;
pmsalt =pmsalct*1.e-3;
pmwatt =pmwatct*1.e-3;

% TOTAL in ug/m3 (dry mass)
pmtotal=(pmorgct+pmecbct+pmsulct+pmxxxct+pmamoct+pmnitct)*1e-3;
%pmtotal=(pmorgct+pmecbct+pmwatct+pmsulct+pmxxxct+pmamoct+pmnitct)*1e-3;

% TOTAL number conc in #/cm3
pnctotal=1.e-6*(total_n(:,3)+total_n(:,4)+total_n(:,5)+total_n(:,6)+total_n(:,7));
% 10-25nm
pncnucl =1.e-6*(total_n(:,8) -total_n(:,3));
% 25-100nm
pncaitk =1.e-6*total_n(:,9);
% >100nm (AS+CS)
pncaccu =1.e-6*(total_n(:,6)+total_n(:,7))

% GAS PHASE ORGANICS in ug/m^3
covtot=plume(:,9)+plume(:,10);

% Coagulation Sink
coags=total_n(:,11);
conds=total_n(:,12);

fsize=12;

figure(1)
plot(time(1:tmax),pmtotal(1:tmax),'-b','LineWidth',2.1)
hold
ax=gca;
set(ax,'linewidth',1.5,'fontsize',fsize)
plot(time(1:tmax),pmecbt(1:tmax),'--r','LineWidth',2.1)
plot(time(1:tmax),pmorgt(1:tmax),'--m','LineWidth',2.1)
plot(time(1:tmax),pmsulf(1:tmax)+pmanit(1:tmax),'--g','LineWidth',2.1)
xlabel('Time (hours)','FontSize',fsize,'FontName','Arial')
ylabel('Particle mass conc. (\mug m^{-3})','FontSize',fsize,'FontName','Arial')
title('Particle mass concentrations','FontSize',fsize)
legend('total','EC','C22+C28','SO4+NO3','Location','NorthEast')
tx=11.0:0.01666:11.084;
set(gca,'xtick',tx);
txtext=['00:00';'00:01';'00:02';'00:03';'00:04';'00:05'];
set(gca, 'xticklabel', txtext)
axis([11.0,11.0833,0,25])
grid on
%uncomment below line to save as jpg
print -djpg '../afigs_plume/partmass-ma16.jpg'


figure(2)
plot(time(1:tmax),pnctotal(1:tmax)/pnctotal(1),'-k','LineWidth',1.1)
hold
ax=gca;
set(ax,'linewidth',1.5,'fontsize',fsize)
plot(time(1:tmax),pncnucl(1:tmax)/pncnucl(1),'-r','LineWidth',1.1)
plot(time(1:tmax),pncaitk(1:tmax)/pncaitk(1),'-m','LineWidth',1.1)
plot(time(1:tmax),pncaccu(1:tmax)/pncaccu(1),'-b','LineWidth',1.1)
plot(tmaap,pncsmps,'ks','MarkerSize',7.,'LineWidth',1.4)
plot(tmaap,pncl25smps,'rs','MarkerSize',7.,'LineWidth',1.4)
plot(tmaap,pnc25_100smps,'ms','MarkerSize',7.,'LineWidth',1.4)
plot(tmaap,pncgt100smps,'bs','MarkerSize',7.,'LineWidth',1.4)
%plot(time(1:tmax),pncbg(1:tmax)/pnctotal(1),'-g','LineWidth',2.1)
%title('MAFOR: only dilution','FontSize',16)
title('MAFOR: all processes','FontSize',fsize)
xlabel('Distance to road (m)','FontSize',fsize)
ylabel('Relative number concentration','FontSize',fsize)
legend('PNC model','PNC 10-25nm model','PNC 25-100nm model','PNC>100nm model','PNC SMPS','PNC<25nm SMPS','PNC 25-100nm SMPS','PNC>100nm SMPS','Location','SouthWest')
tx=10.998612:0.002776:11.015268;
set(gca,'xtick',tx);
txtext=[{'0','50','100','150','200','250','300'}];
set(gca, 'xticklabel', txtext)
axis([10.99861,11.015268,0,1.05])
grid on
%uncomment below line to save as jpg
print -djpg '../afigs_plume/pncrel-ma16-all.jpg'


figure(3)
semilogy(time(1:tmax),conds(1:tmax),'-b','LineWidth',3.1)
hold
ax=gca;
set(ax,'linewidth',1.5,'fontsize',fsize)
semilogy(time(1:tmax),coags(1:tmax),'-r','LineWidth',3.1)
xlabel('Time (hours)','FontSize',14,'FontName','Arial')
ylabel('Coag and Cond Sink (s^{-1})','FontSize',fsize,'FontName','Arial')
legend('CondS','CoagS','Location','NorthEast')
tx=11.0:0.01666:11.084;
set(gca,'xtick',tx);
txtext=['00:00';'00:01';'00:02';'00:03';'00:04';'00:05'];
set(gca, 'xticklabel', txtext)
axis([11.0,11.0833,1e-4,1.])
grid on
%uncomment below line to save as jpg
print -djpg '../afigs_plume/coagsink-ma16.jpg'


figure(4)
plot(time(2:tmax),pmorgt(2:tmax)+covtot(2:tmax),'-b','LineWidth',2.1)
hold
ax=gca;
set(ax,'linewidth',1.5,'fontsize',fsize)
plot(time(2:tmax),covtot(2:tmax),'--r','LineWidth',2.1)
plot(time(2:tmax),pmorgt(2:tmax),'--m','LineWidth',2.1)
xlabel('Time (hours)','FontSize',fsize,'FontName','Arial')
ylabel('Concentration (\mug m^{-3})','FontSize',fsize,'FontName','Arial')
title('Budget COV+OCp','FontSize',fsize)
legend('total','C22+C28(gas)','C22+C28 (p)','Location','NorthEast')
tx=11.0:0.01666:11.084;
set(gca,'xtick',tx);
txtext=['00:00';'00:01';'00:02';'00:03';'00:04';'00:05'];
set(gca, 'xticklabel', txtext)
axis([11.0,11.0833,0,10])
grid on
%uncomment below line to save as jpg
print -djpg '../afigs_plume/ocbudget-ma16.jpg'

