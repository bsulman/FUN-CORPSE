% Figure 1: Conceptual diagram (not Matlab)

%%
% Figure 2: Validation figure
clear
close all

model_output_dir='../../Model-output';

d=load([model_output_dir '/' 'FUN-CORPSE-outputs-NPP-1.0.mat']);

f2=figure(2);clf
f2.Position = [440    33   560   765];

% (A) Root exudation and fungal biomass at endmembers
subplot(211);cla;hold on

per_ECM=d.per_ECM;

load('../Obs-data/Daily_FUN_C_allocation_output.mat','AM_rhizo_percent','ECM_rhizo_percent');

if size(per_ECM,2)>1
    per_ECM=d.per_ECM(:,2);
end

exud_obs=readtable('../Obs-data/exudation.csv');

fs=15;
t=1:length(d.totalrhizoCflux);
ty=t/365;

s=1;e=365*2;
xx=s:e;

[~,ii]=sort(per_ECM);

model=plot(per_ECM(ii),mean(d.totalrhizoCflux(xx,ii)+d.totalfungalprod(xx,ii))*365*1e3,'k-','LineWidth',2.0);
obs=plot(per_ECM,exud_obs.TotalExudation_gC_m2_yr_,'s','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k');

title('(a): Rhizosphere C allocation','FontSize',fs)
xlabel('Percent ECM','FontSize',fs)
ylabel('C flux (gC m^{-2} year^{-1})','FontSize',fs)

h=legend([model,obs],{'Mod','Obs'});
h.set('FontSize',fs);
h.EdgeColor='w';
h.Location='northwest';
set(gca,'FontSize',fs);

%
% (B) Rhizosphere stimulation at endmembers
plots=1:45;

total_bulk=add_outputs(d.ECM_bulk_outputs_e,d.AM_bulk_outputs_e);
total_rhiz=add_outputs(d.ECM_rhiz_outputs_e,d.AM_rhiz_outputs_e);
total_soil=add_outputs(total_bulk,total_rhiz);

total_bulk_noe=add_outputs(d.ECM_bulk_outputs_noe,d.AM_bulk_outputs_noe);
total_rhiz_noe=add_outputs(d.ECM_rhiz_outputs_noe,d.AM_rhiz_outputs_noe);
total_soil_noe=add_outputs(total_bulk_noe,total_rhiz_noe);

total_bulk=total_rhiz_noe;

AMplots=find(per_ECM<=20);ECMplots=find(per_ECM>=80);
rhiz_frac=repmat((AM_rhizo_percent+ECM_rhizo_percent)/100,[20,1]);

rhiz_norm=total_rhiz.protectedC.deadmic;
bulk_norm=total_bulk.protectedC.deadmic;

AMplots=find(per_ECM<=20);ECMplots=find(per_ECM>=80);
mfrac_e_AM=nanmean(total_rhiz.livingMicrobeC(xx,AMplots),2)./nanmean(rhiz_norm(xx,AMplots),2);
mfrac_noe_AM=nanmean(total_bulk.livingMicrobeC(xx,AMplots),2)./nanmean(bulk_norm(xx,AMplots),2);

fast_noe_AM=nanmean(total_bulk.unprotectedC.fast(xx,AMplots),2)./nanmean(bulk_norm(xx,AMplots),2);
fast_e_AM=nanmean(total_rhiz.unprotectedC.fast(xx,AMplots),2)./nanmean(rhiz_norm(xx,AMplots),2);
decomp_e_AM=nanmean(total_rhiz.CO2prod(xx,AMplots),2)./nanmean(rhiz_norm(xx,AMplots),2);
decomp_noe_AM=nanmean(total_bulk.CO2prod(xx,AMplots),2)./nanmean(bulk_norm(xx,AMplots),2);
Nmin_e_AM=nanmean(total_rhiz.N_mineralization(xx,AMplots),2)./nanmean(rhiz_norm(xx,AMplots),2);
Nmin_noe_AM=nanmean(total_bulk.N_mineralization(xx,AMplots),2)./nanmean(bulk_norm(xx,AMplots),2);
Nimm_e_AM=nanmean(total_rhiz.N_immobilization(xx,AMplots),2)./nanmean(rhiz_norm(xx,AMplots),2);
Nimm_noe_AM=nanmean(total_bulk.N_immobilization(xx,AMplots),2)./nanmean(bulk_norm(xx,AMplots),2);
Ndecomp_e_AM=nanmean(total_rhiz.N_decomposed.slow(xx,AMplots),2)./nanmean(rhiz_norm(xx,AMplots),2);
Ndecomp_noe_AM=nanmean(total_bulk.N_decomposed.slow(xx,AMplots),2)./nanmean(bulk_norm(xx,AMplots),2);

mfrac_e_ECM=nanmean(total_rhiz.livingMicrobeC(xx,ECMplots),2)./nanmean(rhiz_norm(xx,ECMplots),2);
mfrac_noe_ECM=nanmean(total_bulk.livingMicrobeC(xx,ECMplots),2)./nanmean(bulk_norm(xx,ECMplots),2);

fast_noe_ECM=nanmean(total_bulk.unprotectedC.fast(xx,ECMplots),2)./nanmean(bulk_norm(xx,ECMplots),2);
fast_e_ECM=nanmean(total_rhiz.unprotectedC.fast(xx,ECMplots),2)./nanmean(rhiz_norm(xx,ECMplots),2);
decomp_e_ECM=nanmean(total_rhiz.CO2prod(xx,ECMplots),2)./nanmean(rhiz_norm(xx,ECMplots),2);
decomp_noe_ECM=nanmean(total_bulk.CO2prod(xx,ECMplots),2)./nanmean(bulk_norm(xx,ECMplots),2);
Nmin_e_ECM=nanmean(total_rhiz.N_mineralization(xx,ECMplots),2)./nanmean(rhiz_norm(xx,ECMplots),2);
Nmin_noe_ECM=nanmean(total_bulk.N_mineralization(xx,ECMplots),2)./nanmean(bulk_norm(xx,ECMplots),2);
Nimm_e_ECM=nanmean(total_rhiz.N_immobilization(xx,ECMplots),2)./nanmean(rhiz_norm(xx,ECMplots),2);
Nimm_noe_ECM=nanmean(total_bulk.N_immobilization(xx,ECMplots),2)./nanmean(bulk_norm(xx,ECMplots),2);
Ndecomp_e_ECM=nanmean(total_rhiz.N_decomposed.slow(xx,ECMplots),2)./nanmean(rhiz_norm(xx,ECMplots),2);
Ndecomp_noe_ECM=nanmean(total_bulk.N_decomposed.slow(xx,ECMplots),2)./nanmean(bulk_norm(xx,ECMplots),2);

ax=subplot(212);cla;hold on
fast_ratio_AM=mean(fast_e_AM)./mean(fast_noe_AM)*100-100;fast_ratio_ECM=mean(fast_e_ECM)./mean(fast_noe_ECM)*100-100;
mfrac_ratio_AM=mean(mfrac_e_AM)./mean(mfrac_noe_AM)*100-100;mfrac_ratio_ECM=mean(mfrac_e_ECM)./mean(mfrac_noe_ECM)*100-100;
decomp_ratio_AM=mean(decomp_e_AM)./mean(decomp_noe_AM)*100-100;decomp_ratio_ECM=mean(decomp_e_ECM)./mean(decomp_noe_ECM)*100-100;
Nmin_ratio_AM=mean(Nmin_e_AM)./mean(Nmin_noe_AM)*100-100;Nmin_ratio_ECM=mean(Nmin_e_ECM)./mean(Nmin_noe_ECM)*100-100;

AM_Cmin_obs=readtable('../Obs-data/AM_Cmin.csv','ReadVariableNames',false);
ECM_Cmin_obs=readtable('../Obs-data/ECM_Cmin.csv','ReadVariableNames',false);
AM_Nmin_obs=readtable('../Obs-data/AM_Nmin.csv','ReadVariableNames',false);
ECM_Nmin_obs=readtable('../Obs-data/ECM_Nmin.csv','ReadVariableNames',false);

plot(decomp_ratio_AM,Nmin_ratio_AM,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8)
plot(decomp_ratio_ECM,Nmin_ratio_ECM,'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',8)
errorbar(mean(AM_Cmin_obs.Var2),mean(AM_Nmin_obs.Var2),std(AM_Nmin_obs.Var2)/2,'Color','k','Marker','s','MarkerSize',8,'MarkerFaceColor','k')
errorbar(mean(ECM_Cmin_obs.Var2),mean(ECM_Nmin_obs.Var2),std(ECM_Nmin_obs.Var2)/2,'Color','k','Marker','s','MarkerSize',8,'MarkerFaceColor','w')

% NOTE: This uses herrorbar, which is available on Matlab file exchange:
% https://www.mathworks.com/matlabcentral/fileexchange/3963-herrorbar
herrorbar(mean(AM_Cmin_obs.Var2),mean(AM_Nmin_obs.Var2),std(AM_Cmin_obs.Var2)/2,'k')
herrorbar(mean(ECM_Cmin_obs.Var2),mean(ECM_Nmin_obs.Var2),std(ECM_Cmin_obs.Var2)/2,'k')

xlim([0,80]);ylim([0,80])

xlabel('C mineralization stimulation (%)')
ylabel('N mineralization stimulation (%)')

title('(b): Rhizosphere stimulation','FontSize',fs)
h=legend('Mod AM','Mod ECM','Obs AM','Obs ECM');
h.EdgeColor='w';
h.Location='northwest';
set(gca,'FontSize',fs);


%% Figure 3: Ambient CO2 soil C, N mineralization, and C mineralization with and without exudation

f3=figure(3);
f3.Position=[843,385,1000,1000];
clf
% d=load('FUN-CORPSE-outputs-NPP-1.0.mat');

if size(per_ECM,2)>1
    per_ECM=d.per_ECM(:,2);
end

d_reversed=load([model_output_dir '/' 'FUN-CORPSE-outputs-NPP-1.0-myc-reversed.mat']);

% Colors
c3=[1.0,0.0,0.0];
c2=[0.5,0.0,0.5];
c1=[0.0,0.0,1.0];
c4=[0.0,1.0,0.0];

total_bulk=add_outputs(d.ECM_bulk_outputs_e,d.AM_bulk_outputs_e);
total_rhiz=add_outputs(d.ECM_rhiz_outputs_e,d.AM_rhiz_outputs_e);
total_soil=add_outputs(total_bulk,total_rhiz);

[totalsoilC,totalsoilN]=total_pool(total_soil,true);
unprotectedC=total_soil.unprotectedC.fast+total_soil.unprotectedC.slow+total_soil.unprotectedC.deadmic;
protectedC=total_soil.protectedC.fast+total_soil.protectedC.slow+total_soil.protectedC.deadmic;
unprotectedN=total_soil.unprotectedN.fast+total_soil.unprotectedN.slow+total_soil.unprotectedN.deadmic;
protectedN=total_soil.protectedN.fast+total_soil.protectedN.slow+total_soil.protectedN.deadmic;

total_bulk_noe=add_outputs(d.ECM_bulk_outputs_noe,d.AM_bulk_outputs_noe);
total_rhiz_noe=add_outputs(d.ECM_rhiz_outputs_noe,d.AM_rhiz_outputs_noe);
total_soil_noe=add_outputs(total_bulk_noe,total_rhiz_noe);

[totalsoilC_noe,totalsoilN_noe]=total_pool(total_soil_noe,true);
unprotectedC_noe=total_soil_noe.unprotectedC.fast+total_soil_noe.unprotectedC.slow+total_soil_noe.unprotectedC.deadmic;
protectedC_noe=total_soil_noe.protectedC.fast+total_soil_noe.protectedC.slow+total_soil_noe.protectedC.deadmic;
unprotectedN_noe=total_soil_noe.unprotectedN.fast+total_soil_noe.unprotectedN.slow+total_soil_noe.unprotectedN.deadmic;
protectedN_noe=total_soil_noe.protectedN.fast+total_soil_noe.protectedN.slow+total_soil_noe.protectedN.deadmic;


total_bulk_rev=add_outputs(d_reversed.ECM_bulk_outputs_e,d_reversed.AM_bulk_outputs_e);
total_rhiz_rev=add_outputs(d_reversed.ECM_rhiz_outputs_e,d_reversed.AM_rhiz_outputs_e);
total_soil_rev=add_outputs(total_bulk_rev,total_rhiz_rev);

[totalsoilC_rev,totalsoilN_rev]=total_pool(total_soil_rev,true);
unprotectedC_rev=total_soil_rev.unprotectedC.fast+total_soil_rev.unprotectedC.slow+total_soil_rev.unprotectedC.deadmic;
protectedC_rev=total_soil_rev.protectedC.fast+total_soil_rev.protectedC.slow+total_soil_rev.protectedC.deadmic;
unprotectedN_rev=total_soil_rev.unprotectedN.fast+total_soil_rev.unprotectedN.slow+total_soil_rev.unprotectedN.deadmic;
protectedN_rev=total_soil_rev.protectedN.fast+total_soil_rev.protectedN.slow+total_soil_rev.protectedN.deadmic;

total_bulk_rev_noe=add_outputs(d_reversed.ECM_bulk_outputs_noe,d_reversed.AM_bulk_outputs_noe);
total_rhiz_rev_noe=add_outputs(d_reversed.ECM_rhiz_outputs_noe,d_reversed.AM_rhiz_outputs_noe);
total_soil_rev_noe=add_outputs(total_bulk_rev_noe,total_rhiz_rev_noe);

[totalsoilC_rev_noe,totalsoilN_rev_noe]=total_pool(total_soil_rev_noe,true);
unprotectedC_rev_noe=total_soil_rev_noe.unprotectedC.fast+total_soil_rev_noe.unprotectedC.slow+total_soil_rev_noe.unprotectedC.deadmic;
protectedC_rev_noe=total_soil_rev_noe.protectedC.fast+total_soil_rev_noe.protectedC.slow+total_soil_rev_noe.protectedC.deadmic;
unprotectedN_rev_noe=total_soil_rev_noe.unprotectedN.fast+total_soil_rev_noe.unprotectedN.slow+total_soil_rev_noe.unprotectedN.deadmic;
protectedN_rev_noe=total_soil_rev_noe.protectedN.fast+total_soil_rev_noe.protectedN.slow+total_soil_rev_noe.protectedN.deadmic;

s2=365*14+1;e2=365*19;
s3=365*44+1;e3=365*49;
s=365*4+1;e=365*5;

xx=s:e;
xx2=s2:e2;
xx3=s3:e3;
fs=15;

[~,ii]=sort(per_ECM);


subplot(331);cla;hold on
soilC=plot(per_ECM(ii),mean(totalsoilC(xx,ii),1),'ko','MarkerFaceColor',[0.2,0.2,0.2]);
soilC_noe=plot(per_ECM(ii),mean(totalsoilC_noe(xx,ii),1),'ks','MarkerFaceColor',[0.9,0.9,0.9]);
title('(a): Total C','FontSize',fs);ylabel('Total C (kg m^{-2})','FontSize',fs)
set(gca,'FontSize',fs)

subplot(334);cla;hold on
h_unprotectedC=plot(per_ECM(ii),mean(unprotectedC(xx,ii),1),'ko','MarkerFaceColor',[0.2,0.2,0.2]);
h_unprotectedC_noe=plot(per_ECM(ii),mean(unprotectedC_noe(xx,ii),1),'ks','MarkerFaceColor',[0.9,0.9,0.9]);
title('(b): Unprotected C','FontSize',fs);ylabel('Unprotected C (kg m^{-2})','FontSize',fs)
set(gca,'FontSize',fs)

subplot(337);cla;hold on
h_protectedC=plot(per_ECM(ii),mean(protectedC(xx,ii),1),'ko','MarkerFaceColor',[0.2,0.2,0.2]);
h_protectedC_noe=plot(per_ECM(ii),mean(protectedC_noe(xx,ii),1),'ks','MarkerFaceColor',[0.9,0.9,0.9]);
title('(c): Protected C','FontSize',fs);ylabel('Protected C (kg m^{-2})','FontSize',fs)
set(gca,'FontSize',fs)
xlabel('Percent ECM','FontSize',fs)


subplot(332);cla;hold on
soilN=plot(per_ECM(ii),mean(totalsoilN(xx,ii),1),'ko','MarkerFaceColor',[0.2,0.2,0.2]);
soilN_noe=plot(per_ECM(ii),mean(totalsoilN_noe(xx,ii),1),'ks','MarkerFaceColor',[0.9,0.9,0.9]);
title('(d): Total N','FontSize',fs);ylabel('Total N (kg m^{-2})','FontSize',fs)
set(gca,'FontSize',fs)

subplot(335);cla;hold on
h_unprotectedN=plot(per_ECM(ii),mean(unprotectedN(xx,ii),1),'ko','MarkerFaceColor',[0.2,0.2,0.2]);
h_unprotectedN_noe=plot(per_ECM(ii),mean(unprotectedN_noe(xx,ii),1),'ks','MarkerFaceColor',[0.9,0.9,0.9]);
title('(e): Unprotected N','FontSize',fs);ylabel('Unprotected N (kg m^{-2})','FontSize',fs)
set(gca,'FontSize',fs)

subplot(338);cla;hold on
h_protectedN=plot(per_ECM(ii),mean(protectedN(xx,ii),1),'ko','MarkerFaceColor',[0.2,0.2,0.2]);
h_protectedN_noe=plot(per_ECM(ii),mean(protectedN_noe(xx,ii),1),'ks','MarkerFaceColor',[0.9,0.9,0.9]);
title('(f): Protected N','FontSize',fs);ylabel('Protected N (kg m^{-2})','FontSize',fs)
set(gca,'FontSize',fs)
xlabel('Percent ECM')


subplot(333);cla;hold on

yy=mod(t,365)>150&mod(t,365)<=350&t>s&t<e;

t2=1:length(d_reversed.totalrhizoCflux);
yy2=mod(t2,365)>150&mod(t2,365)<=350&t2>s2&t2<e2;
yy3=mod(t2,365)>150&mod(t2,365)<=350&t2>s3&t2<e3;

yy=xx;yy2=xx2;yy3=xx3;

Nmin_e=mean(total_soil.N_mineralization(yy,:),1)*1e3;
Nmin_noe=mean(total_soil_noe.N_mineralization(yy,:),1)*1e3;

withe=plot(per_ECM(ii),Nmin_e(ii)*365,'ko','MarkerFaceColor',[0.2,0.2,0.2]);
noe=plot(per_ECM(ii),Nmin_noe(ii)*365,'ks','MarkerFaceColor',[0.9,0.9,0.9]);

h=legend([withe,noe],{'Coupled rhizosphere','De-coupled rhizosphere'});

h.EdgeColor='w';
title('(g): N mineralization','FontSize',fs);
ylabel('N min (g m^{-2} year^{-1})','FontSize',fs)

set(h,'FontSize',fs*0.8);
set(gca,'FontSize',fs);

turnover_e=mean(d.total_outputs_e.decomp.slow(yy,plots)./d.total_outputs_e.unprotectedC.slow(yy,plots));
turnover_noe=mean(d.total_outputs_noe.decomp.slow(yy,plots)./d.total_outputs_noe.unprotectedC.slow(yy,plots));

subplot(336);cla;hold on
withe=plot(per_ECM(ii),turnover_e(ii),'ko','MarkerFaceColor',[0.2,0.2,0.2]);
noe=plot(per_ECM(ii),turnover_noe(ii),'ks','MarkerFaceColor',[0.9,0.9,0.9]);

ylabel('Slow C turnover rate (yr^{-1})','FontSize',fs)

title('(h): Soil C turnover rate','FontSize',fs)

set(gca,'FontSize',fs)
ylim([0,0.2])


subplot(339);cla;hold on
unprotectedC=total_soil.unprotectedC.fast+total_soil.unprotectedC.slow+total_soil.unprotectedC.deadmic;
microbefrac_e=mean(total_soil.livingMicrobeC(yy,plots));

unprotectedC_noe=total_soil_noe.unprotectedC.fast+total_soil_noe.unprotectedC.slow+total_soil_noe.unprotectedC.deadmic;
microbefrac_noe=mean(total_soil_noe.livingMicrobeC(yy,plots));



withe=plot(per_ECM(ii),microbefrac_e(ii)*1000,'ko','MarkerFaceColor',[0.2,0.2,0.2]);
noe=plot(per_ECM(ii),microbefrac_noe(ii)*1000,'ks','MarkerFaceColor',[0.9,0.9,0.9]);

ylabel('Microbial biomass (gC m^{-2})','FontSize',fs)
xlabel('Percent ECM','FontSize',fs)
title('(i): Microbial biomass','FontSize',fs)
set(gca,'FontSize',fs)

% Figure 5: Change over time at end members with reversed mycorrhizae
f5=figure(5);
f5.Position=[843,385,1000,1000];
t=(1:length(totalsoilC))/365;
t_rev=(1:length(totalsoilC_rev))/365;

subplot(331);cla;hold on
plot(t,smooth(mean(totalsoilC(:,AMplots),2),365),'b-','LineWidth',2.0);
plot(t,smooth(mean(totalsoilC(:,ECMplots),2),365),'r-','LineWidth',2.0);
plot(t_rev,smooth(mean(totalsoilC_rev(:,AMplots),2),365),'b--','LineWidth',2.0);
plot(t_rev,smooth(mean(totalsoilC_rev(:,ECMplots),2),365),'r--','LineWidth',2.0);
title('(a): Total C','FontSize',fs);ylabel('Total C (kg m^{-2})','FontSize',fs)
set(gca,'FontSize',fs)
xlim([0,50])

subplot(334);cla;hold on
plot(t,smooth(mean(unprotectedC(:,AMplots),2),365),'b-','LineWidth',2.0);
plot(t,smooth(mean(unprotectedC(:,ECMplots),2),365),'r-','LineWidth',2.0);
plot(t_rev,smooth(mean(unprotectedC_rev(:,AMplots),2),365),'b--','LineWidth',2.0);
plot(t_rev,smooth(mean(unprotectedC_rev(:,ECMplots),2),365),'r--','LineWidth',2.0);
title('(b): Unprotected C','FontSize',fs);ylabel('Unprotected C (kg m^{-2})','FontSize',fs)
set(gca,'FontSize',fs)
xlim([0,50])

subplot(337);cla;hold on
plot(t,mean(protectedC(:,AMplots),2),'b-','LineWidth',2.0);
plot(t_rev,mean(protectedC_rev(:,AMplots),2),'b--','LineWidth',2.0);
plot(t,mean(protectedC(:,ECMplots),2),'r-','LineWidth',2.0);
plot(t_rev,mean(protectedC_rev(:,ECMplots),2),'r--','LineWidth',2.0);
title('(c): Protected C','FontSize',fs);ylabel('Protected C (kg m^{-2})','FontSize',fs)
set(gca,'FontSize',fs)
xlabel('Year','FontSize',fs)
xlim([0,50])


subplot(332);cla;hold on
plot(t,smooth(mean(totalsoilN(:,AMplots),2),365),'b-','LineWidth',2.0);
plot(t,smooth(mean(totalsoilN(:,ECMplots),2),365),'r-','LineWidth',2.0);
plot(t_rev,smooth(mean(totalsoilN_rev(:,AMplots),2),365),'b--','LineWidth',2.0);
plot(t_rev,smooth(mean(totalsoilN_rev(:,ECMplots),2),365),'r--','LineWidth',2.0);
title('(d): Total N','FontSize',fs);ylabel('Total N (kg m^{-2})','FontSize',fs)
set(gca,'FontSize',fs)
h=legend('AM Veg, AM Soil','ECM Veg, ECM Soil','ECM Veg, AM Soil','AM Veg, ECM Soil');
set(h,'Orientation','horizontal');
h.EdgeColor='w';
h.Position=[.46,.964,.447,.027];
xlim([0,50])

subplot(335);cla;hold on
plot(t,smooth(mean(unprotectedN(:,AMplots),2),365),'b-','LineWidth',2.0);
plot(t,smooth(mean(unprotectedN(:,ECMplots),2),365),'r-','LineWidth',2.0);
plot(t_rev,smooth(mean(unprotectedN_rev(:,AMplots),2),365),'b--','LineWidth',2.0);
plot(t_rev,smooth(mean(unprotectedN_rev(:,ECMplots),2),365),'r--','LineWidth',2.0);
title('(e): Unprotected N','FontSize',fs);ylabel('Unprotected N (kg m^{-2})','FontSize',fs)
set(gca,'FontSize',fs)
xlim([0,50])

subplot(338);cla;hold on
plot(t,mean(protectedN(:,AMplots),2),'b-','LineWidth',2.0);
plot(t,mean(protectedN(:,ECMplots),2),'r-','LineWidth',2.0);
plot(t_rev,mean(protectedN_rev(:,AMplots),2),'b--','LineWidth',2.0);
plot(t_rev,mean(protectedN_rev(:,ECMplots),2),'r--','LineWidth',2.0);
title('(f): Protected N','FontSize',fs);ylabel('Protected N (kg m^{-2})','FontSize',fs)
set(gca,'FontSize',fs)
xlabel('Year')
xlim([0,50])

subplot(333);cla;hold on

Nmin_AM=smooth(mean(total_soil.N_mineralization(:,AMplots),2)*1e3,365);
Nmin_ECM=smooth(mean(total_soil.N_mineralization(:,ECMplots),2)*1e3,365);
Nmin_AM_rev=smooth(mean(total_soil_rev.N_mineralization(:,AMplots),2)*1e3,365);
Nmin_ECM_rev=smooth(mean(total_soil_rev.N_mineralization(:,ECMplots),2)*1e3,365);

x=183:365*50-183;

plot(t(x),Nmin_AM(x)*365,'b-','LineWidth',2.0);
plot(t(x),Nmin_ECM(x)*365,'r-','LineWidth',2.0);
plot(t(x),Nmin_AM_rev(x)*365,'b--','LineWidth',2.0);
plot(t(x),Nmin_ECM_rev(x)*365,'r--','LineWidth',2.0);
title('(g): N mineralization','FontSize',fs)
ylabel('N min (g m^{-2} year^{-1})','FontSize',fs)
set(gca,'FontSize',fs);
xlim([0,50])


turnover_AM=smooth(mean(d.total_outputs_e.decomp.slow(:,AMplots)./d.total_outputs_e.unprotectedC.slow(:,AMplots),2),365);
turnover_ECM=smooth(mean(d.total_outputs_e.decomp.slow(:,ECMplots)./d.total_outputs_e.unprotectedC.slow(:,ECMplots),2),365);

turnover_AM_rev=smooth(mean(d_reversed.total_outputs_e.decomp.slow(:,AMplots)./d_reversed.total_outputs_e.unprotectedC.slow(:,AMplots),2),365);
turnover_ECM_rev=smooth(mean(d_reversed.total_outputs_e.decomp.slow(:,ECMplots)./d_reversed.total_outputs_e.unprotectedC.slow(:,ECMplots),2),365);

subplot(336);cla;hold on
plot(t(x),turnover_AM(x),'b-','LineWidth',2.0);
plot(t(x),turnover_ECM(x),'r-','LineWidth',2.0);
plot(t(x),turnover_AM_rev(x),'b--','LineWidth',2.0);
plot(t(x),turnover_ECM_rev(x),'r--','LineWidth',2.0);
ylabel('Slow C turnover rate (yr^{-1})','FontSize',fs)
title('(h): Soil C turnover rate','FontSize',fs)
set(gca,'FontSize',fs)
xlim([0,50])


subplot(339);cla;hold on
unprotectedC=total_soil.unprotectedC.fast+total_soil.unprotectedC.slow+total_soil.unprotectedC.deadmic;
microbefrac_AM=smooth(mean(total_soil.livingMicrobeC(:,AMplots),2),365);
microbefrac_ECM=smooth(mean(total_soil.livingMicrobeC(:,ECMplots),2),365);

unprotectedC_rev=total_soil_rev.unprotectedC.fast+total_soil_rev.unprotectedC.slow+total_soil_rev.unprotectedC.deadmic;
microbefrac_AM_rev=smooth(mean(total_soil_rev.livingMicrobeC(:,AMplots),2),365);
microbefrac_ECM_rev=smooth(mean(total_soil_rev.livingMicrobeC(:,ECMplots),2),365);

plot(t(x),microbefrac_AM(x)*1000,'b-','LineWidth',2.0);
plot(t(x),microbefrac_ECM(x)*1000,'r-','LineWidth',2.0);
plot(t(x),microbefrac_AM_rev(x)*1000,'b--','LineWidth',2.0);
plot(t(x),microbefrac_ECM_rev(x)*1000,'r--','LineWidth',2.0);
ylabel('Microbial biomass (gC m^{-2})','FontSize',fs)
xlabel('Year','FontSize',fs)
title('(i): Microbial biomass','FontSize',fs)
set(gca,'FontSize',fs)
xlim([0,50])


% Figure 4: Soil decomposition per unit N mineralization

d1=d;

total_bulk1=add_outputs(d1.ECM_bulk_outputs_e,d1.AM_bulk_outputs_e);
total_rhiz1=add_outputs(d1.ECM_rhiz_outputs_e,d1.AM_rhiz_outputs_e);
total_soil1=add_outputs(total_bulk1,total_rhiz1);
total_bulk1_noe=add_outputs(d1.ECM_bulk_outputs_noe,d1.AM_bulk_outputs_noe);
total_rhiz1_noe=add_outputs(d1.ECM_rhiz_outputs_noe,d1.AM_rhiz_outputs_noe);
total_soil1_noe=add_outputs(total_bulk1_noe,total_rhiz1_noe);

Nmin_e=mean(total_soil1.N_mineralization(1:365,:),1)*1e3*365;%gN/m2/year
Nmin_noe=mean(total_soil1_noe.N_mineralization(1:365,:),1)*1e3*365;

decomp_e=mean(d1.total_outputs_e.decomp.slow(1:365,plots))*1e3;
decomp_noe=mean(d1.total_outputs_noe.decomp.slow(1:365,plots))*1e3;

figure(4);clf

subplot(111);cla;hold on
plot(per_ECM(ii),(decomp_e(ii)-decomp_noe(ii))./(Nmin_e(ii)-Nmin_noe(ii)),'ko','MarkerFaceColor','k')
title('SOC loss per unit extra N mineralization','FontSize',fs)
xlabel('Percent ECM','FontSize',fs)
ylabel('SOC loss per unit N min (gC gN^{-1})','FontSize',fs)
set(gca,'FontSize',fs)


%% Figure 6: Differences with elevated NPP (actual difference from control, not percent diff)

d_control=load([model_output_dir '/' 'FUN-CORPSE-outputs-NPP-1.0.mat']);
d_12=load([model_output_dir '/' 'FUN-CORPSE-outputs-NPP-1.2.mat']);
d_highlitter=load([model_output_dir '/' 'FUN-CORPSE-outputs-NPP-1.2-same-exud.mat']);


if size(per_ECM,2)>1
    per_ECM=d.per_ECM(:,2);
end


t=1:length(d_control.totalrhizoCflux);


[~,ii]=sort(per_ECM);

s=365*7+1;e=365*9;
xx=s:e;
t=1:length(d_control.totalrhizoCflux);

yy=mod(t,365)>250&mod(t,365)<=350&t>365*7&t<365*10;
yy=find(yy);


f6=figure(6);clf
f6.Position=[606    21   1000   784];
[~,ii]=sort(per_ECM);

s=365*7+1;e=365*9;
xx=s:e;
t=1:length(d_control.totalrhizoCflux);

d2=d_12;
d=d_highlitter;

% Root exudation, ambient and elevated
subplot(321);cla;hold on

exud_ratio=plot(per_ECM(ii),(mean(d2.totalrhizoCflux(xx,ii)+d2.totalfungalprod(xx,ii))-(mean(d.totalrhizoCflux(xx,ii)+d.totalfungalprod(xx,ii))))*365*1e3,'ko','MarkerFaceColor','k');

title('(a): Rhizosphere C allocation (Absolute change)','FontSize',fs)
xlabel('Percent ECM','FontSize',fs)
ylabel('C flux difference (gC m^{-2} year^{-1})','FontSize',fs)

set(gca,'FontSize',fs)

total_bulk=add_outputs(d.ECM_bulk_outputs_e,d.AM_bulk_outputs_e);
total_rhiz=add_outputs(d.ECM_rhiz_outputs_e,d.AM_rhiz_outputs_e);
total_soil=add_outputs(total_bulk,total_rhiz);

total_bulk2=add_outputs(d2.ECM_bulk_outputs_e,d2.AM_bulk_outputs_e);
total_rhiz2=add_outputs(d2.ECM_rhiz_outputs_e,d2.AM_rhiz_outputs_e);
total_soil2=add_outputs(total_bulk2,total_rhiz2);


% N mineralization, ambient and elevated
subplot(323);cla;hold on

Nmin_e=mean(d.total_outputs_e.N_mineralization(yy,ii),1)*1e3;
Nmin_e2=mean(d2.total_outputs_e.N_mineralization(yy,ii),1)*1e3;

elev_ratio=plot(per_ECM(ii),(Nmin_e2-Nmin_e)*365,'ko','MarkerFaceColor','k');

title('(b): N mineralization','FontSize',fs)
xlabel('Percent ECM','FontSize',fs)
ylabel('N min difference (gN m^{-2} year^{-1})','FontSize',fs)
set(gca,'FontSize',fs)


% Soil C pools
subplot(325);cla;hold on

[totalsoilC,totalsoilN]=total_pool(total_soil,true);
[totalsoilC2,totalsoilN2]=total_pool(total_soil2,true);

elev=plot(per_ECM(ii),(mean(totalsoilC2(xx,ii),1)-mean(totalsoilC(xx,ii),1))*1e3,'ko','MarkerFaceColor','k');

title('(c): Soil C','FontSize',fs)
xlabel('Percent ECM','FontSize',fs)
ylabel('Soil C stock difference (gC m^{-2})','FontSize',fs)
set(gca,'FontSize',fs)


%% Percent differences


% Root exudation, ambient and elevated
subplot(322);cla;hold on
exud_ratio=plot(per_ECM(ii),(mean(d2.totalrhizoCflux(xx,ii)+d2.totalfungalprod(xx,ii))./(mean(d.totalrhizoCflux(xx,ii)+d.totalfungalprod(xx,ii))))*100-100,'ko','MarkerFaceColor','k');

title('(d): Rhizosphere C allocation (Percent change)','FontSize',fs)
xlabel('Percent ECM','FontSize',fs)
ylabel('C flux change (%)','FontSize',fs)

set(gca,'FontSize',fs)

total_bulk=add_outputs(d.ECM_bulk_outputs_e,d.AM_bulk_outputs_e);
total_rhiz=add_outputs(d.ECM_rhiz_outputs_e,d.AM_rhiz_outputs_e);
total_soil=add_outputs(total_bulk,total_rhiz);

total_bulk2=add_outputs(d2.ECM_bulk_outputs_e,d2.AM_bulk_outputs_e);
total_rhiz2=add_outputs(d2.ECM_rhiz_outputs_e,d2.AM_rhiz_outputs_e);
total_soil2=add_outputs(total_bulk2,total_rhiz2);

% N mineralization, ambient and elevated
subplot(324);cla;hold on
Nmin_e=mean(total_soil.N_mineralization(yy,ii),1)*1e3;
Nmin_e2=mean(total_soil2.N_mineralization(yy,ii),1)*1e3;

Nmin_e=mean(d.total_outputs_e.N_mineralization(yy,ii),1)*1e3;
Nmin_e2=mean(d2.total_outputs_e.N_mineralization(yy,ii),1)*1e3;

elev_ratio=plot(per_ECM(ii),(Nmin_e2./Nmin_e*100-100),'ko','MarkerFaceColor','k');

title('(e): N mineralization','FontSize',fs)
xlabel('Percent ECM','FontSize',fs)
ylabel('N min difference (%)','FontSize',fs)
set(gca,'FontSize',fs)

% Soil C pools
subplot(326);cla;hold on

[totalsoilC,totalsoilN]=total_pool(total_soil,true);
[totalsoilC2,totalsoilN2]=total_pool(total_soil2,true);

[totalsoilC_unprotected,totalsoilN_unprotected]=total_pool(total_soil,false);
[totalsoilC2_unprotected,totalsoilN2_unprotected]=total_pool(total_soil2,false);


elev=plot(per_ECM(ii),mean(totalsoilC2(xx,ii),1)./mean(totalsoilC(xx,ii),1)*100-100,'ko','MarkerFaceColor','k');

title('(f): Soil C','FontSize',fs)
xlabel('Percent ECM','FontSize',fs)
ylabel('Soil C stock change (%)','FontSize',fs)
set(gca,'FontSize',fs)
