%This code takes the data from the mycorrhizal gradient plots in Southern
%Indiana to generate daily NPP and daily N demand.

%Data Structure for input and output variables is 365 rows and the number
%of columns corresponds to the number of gradient plots (here 45).

addpath ../Plotting-code

load ../Obs-data/Daily_FUN_C_allocation_output.mat
load ../Obs-data/froot_dynamics.mat


tic

if ~exist('reverse_plant_myco','var')
  reverse_plant_myco=false;
end
if ~exist('reverse_myco_costs','var')
  reverse_myco_costs=false;
end

%(2) FIXED CONSTANTS
woodCN=337.6; %Do not have data on a mycorrhizal basis right now.
if ~reverse_plant_myco
  AMrootCN=32.7;%The AM and ECM CN ratios are the mean of plots that are more than 70% of each assoication and are significantly different.
  ECMrootCN=48.2;
  AMleafCN=18.5; %These are the means of leaf samples shot down in 2012.
  ECMleafCN=18.6;
  AMstorageCN=1; %This number is meant to be the C to N ratio of glutamic acid
  ECMstorageCN=1;
else
  ECMrootCN=32.7;%The AM and ECM CN ratios are the mean of plots that are more than 70% of each assoication and are significantly different.
  AMrootCN=48.2;
  ECMleafCN=18.5; %These are the means of leaf samples shot down in 2012.
  AMleafCN=18.6;
  ECMstorageCN=1; %This number is meant to be the C to N ratio of glutamic acid
  AMstorageCN=1;
end

if ~exist('exudate_CN','var')
  exudate_CN=0;
end

if ~exist('cost_param_mult','var')
    cost_param_mult=1.0;
end

%(3) COST PARAMETERS
%Fixation cost parameters (cost_fixN) based off of Houlton et al. 2008
s_FIX = -30;
a_FIX = -3.62;
b_FIX = 0.27;
c_FIX = 25.15;
%Resorb cost parameters (cost_resorb)
kR=0.1/20;

%Mycorrhizal cost parameters
if ~reverse_plant_myco && ~reverse_myco_costs
  AkN=.1/365*cost_param_mult/2; %AM cost
  AkC=.05/365*cost_param_mult*2e6; %AM cost
  EkN=.05/365*cost_param_mult/2; %ECM cost
  EkC=.3/365*cost_param_mult*2.5e6; %ECM cost
else
  AkN=.05/365*cost_param_mult/2; %AM cost
  AkC=.3/365*cost_param_mult*2.5e6; %AM cost
  EkN=.1/365*cost_param_mult/2; %ECM cost
  EkC=.05/365*cost_param_mult*2e6; %ECM cost
end

%Non_Mycorrhizal cost parameters
NkN=0.6/365*cost_param_mult;
NkC=0.01/365*cost_param_mult*2e6;


% 100% and 0% ECM cause issues, so make sure a small percentage exists
for ii=1:45
    if per_ECM(ii,2)==0
        per_ECM(ii,2)=0.1;
    end
    if per_ECM(ii,2)==100
        per_ECM(ii,2)=99.9;
    end
end

%(4)TIMESTEP AND NUMBER OF PLOTS.
if ~exist('nyears','var')
    nyears=1;
end

if exist('NPP_mult','var')
   daily_froot_production=daily_froot_production*NPP_mult;
   fine_roots_npp=fine_roots_npp*NPP_mult;

   daily_leaf_production=daily_leaf_production*NPP_mult;
   daily_wood_production=daily_wood_production*NPP_mult;

end

if exist('litter_mult','var')
    litter_production=litter_production*litter_mult;
    daily_froot_turnover=daily_froot_turnover*litter_mult;
end

if ~exist('rhiz_fungi_frac','var')
    rhiz_fungi_frac=1.0;
end

if ~exist('litter_fungi_frac','var')
  litter_fungi_frac=0.0;
end

if exist('constant_NPP','var')
    if constant_NPP
       daily_froot_production=repmat(mean(daily_froot_production,2),[1,45]);
       fine_roots_npp=repmat(mean(fine_roots_npp,2),[1,45]);
       daily_froot_turnover=repmat(mean(daily_froot_turnover,2),[1,45]);
       daily_leaf_production=repmat(mean(daily_leaf_production,2),[1,45]);
       daily_wood_production=repmat(mean(daily_wood_production,2),[1,45]);
       rootbiomass=repmat(mean(rootbiomass,2),[1,45]);
       ECMrb=repmat(mean(ECMrb,2),[1,45]);
       AMrb=repmat(mean(AMrb,2),[1,45]);
       litter_production=repmat(mean(litter_production,2),[1,45]);
    end
end

if exist('woodNPP_mult','var')
    daily_wood_production=daily_wood_production*woodNPP_mult;
end

total_rhiz=AM_rhizo_percent+ECM_rhizo_percent;
AM_rhizo_percent=total_rhiz;ECM_rhizo_percent=total_rhiz;

if exist('rhiz_mult','var')
    AM_rhizo_percent=AM_rhizo_percent*rhiz_mult;
    ECM_rhizo_percent=ECM_rhizo_percent*rhiz_mult;
    AM_rhizo_percent(AM_rhizo_percent>100)=100;
    ECM_rhizo_percent(AM_rhizo_percent>100)=100;
end

if ~exist('root_inputs_only_rhiz','var')
    root_inputs_only_rhiz=false;
end

if ~exist('atmo_N_dep')
    atmo_N_dep=0.0;
end

if ~reverse_plant_myco
  ECM_litter_composition=[0.1,0.90,0.0];
  ECM_rootlitter_composition=[0.1,0.90,0.0];

  AM_litter_composition=[0.3,0.7,0.0];
  AM_rootlitter_composition=[0.3,0.7,0.0];
else
  AM_litter_composition=[0.1,0.90,0.0];
  AM_rootlitter_composition=[0.1,0.90,0.0];

  ECM_litter_composition=[0.3,0.7,0.0];
  ECM_rootlitter_composition=[0.3,0.7,0.0];
end

number=45; %Number of plots or individual units
timestep=365*nyears; %Daily here


%(4) INGEST DAILY C ALLOCATION DATA and CALCULATE INTEGRATED NPP AND
%NDEMAND THAT INCLUDES RETRANSLOCATION DEMAND DURING SENESCENCE

% Here we are splitting NPP into its AM and ECM fractions and then also
%calculating the N needed to support that NPP

%First calculate amount of AM and ECM root production and N needed
AMrootNdemand=zeros(365,number);
ECMrootNdemand=zeros(365,number);

for i=1:number;
    for k=1:365
    AMrootNdemand(k,i)=((daily_froot_production(k,i)*(1-(per_ECM(i,2)/100)))/AMrootCN)/1000;
    ECMrootNdemand(k,i)=((daily_froot_production(k,i)*((per_ECM(i,2)/100)))/ECMrootCN)/1000;
    end
end

%Do the same for wood
AMwoodNdemand=zeros(365,number);
ECMwoodNdemand=zeros(365,number);

for i=1:number;
    for k=1:365
    AMwoodNdemand(k,i)=((daily_wood_production(k,i)*(1-(per_ECM(i,2)/100)))/woodCN)/1000;
    ECMwoodNdemand(k,i)=((daily_wood_production(k,i)*((per_ECM(i,2)/100)))/woodCN)/1000;
    end
end

%Do the same leaves
AMleafNdemand=zeros(365,number);
ECMleafNdemand=zeros(365,number);

for i=1:number;
    for k=1:365
    AMleafNdemand(k,i)=((daily_leaf_production(k,i)*(1-(per_ECM(i,2)/100)))/AMleafCN)/1000;
    ECMleafNdemand(k,i)=((daily_leaf_production(k,i)*((per_ECM(i,2)/100)))/ECMleafCN)/1000;
    end
end

%Finally do the same to replenish N stores during senescence
AMstorageNdemand=zeros(365,number);
ECMstorageNdemand=zeros(365,number);

AMstorageNdemand(303,:)=(sum(AMleafNdemand,1)*.7);
ECMstorageNdemand(303,:)=(sum(ECMleafNdemand,1)*.7);

%Calculate total N demand for each mycorrhizal type
AMtotalNdemand=zeros(365,number);
ECMtotalNdemand=zeros(365,number);

for i=1:number;
    for k=1:365
    AMtotalNdemand(k,i)=AMwoodNdemand(k,i)+AMleafNdemand(k,i)+AMrootNdemand(k,i);
    ECMtotalNdemand(k,i)=ECMwoodNdemand(k,i)+ECMleafNdemand(k,i)+ECMrootNdemand(k,i);
    end
end


%Calculate NPP for each type
AMNPP=zeros(365,number);
ECMNPP=zeros(365,number);

for i=1:number;
    for k=1:365
    AMNPP(k,i)=((daily_leaf_production(k,i)+daily_wood_production(k,i)+daily_froot_production(k,i))*(1-(per_ECM(i,2)/100)))/1000;
    ECMNPP(k,i)=((daily_leaf_production(k,i)+daily_wood_production(k,i)+daily_froot_production(k,i))*((per_ECM(i,2)/100)))/1000;
    end
end

%Calculate plantCN for each mycorrhizal type
AMplantCN=zeros(365,number);
ECMplantCN=zeros(365,number);

for i=1:number;
    for k=1:365
    AMplantCN(k,i)=AMNPP(k,i)/AMtotalNdemand(k,i);
    if isnan(AMplantCN(k,i))==1;
        AMplantCN(k,i)=0;
    end

    ECMplantCN(k,i)=ECMNPP(k,i)/ECMtotalNdemand(k,i);
    if isnan(ECMplantCN(k,i))==1;
        ECMplantCN(k,i)=0;
    end
    end
end

%Calculate leafN available for retranslocation for each mycorrhizal type
%and the amount of C in litter
AMleafN=zeros(365,number);
ECMleafN=zeros(365,number);

for i=1:number;
    for k=1:365
    AMleafN(k,i)=((litter_production(k,i)*(1-(per_ECM(i,2)/100)))/AMleafCN)/1000;
    ECMleafN(k,i)=((litter_production(k,i)*((per_ECM(i,2)/100)))/ECMleafCN)/1000;
    end
end



AMlitter_production=zeros(365,number);
ECMlitter_production=zeros(365,number);

litter_production=litter_production/1000;

for i=1:number;
    for k=1:365
    AMlitter_production(k,i)=litter_production(k,i)*(1-(per_ECM(i,2)/100));
    ECMlitter_production(k,i)=litter_production(k,i)*((per_ECM(i,2)/100));
    end
end

meanlitter_production=mean([AMlitter_production,ECMlitter_production],2);

%Initialize the N and C storage pools.  The N is used to grow leaves in spring which we
%also determine demand in the fall and the C is used for retranslocation.

AMstorageN=zeros(timestep,number);
ECMstorageN=zeros(timestep,number);
AMstorageC=zeros(365,number);
ECMstorageC=zeros(365,number);

AMstorageN(1,:)=sum(AMleafNdemand,1)*.7;
ECMstorageN(1,:)=sum(ECMleafNdemand,1)*.7;


AMstorageC=sum(AMlitter_production);  % Assuming that the storage pool of C is equal to the leaf litter C pool
ECMstorageC=sum(ECMlitter_production);


%% ECM portion of code
%(5) INITIALIZING MATRICES
ECMpassive=zeros(timestep,number);
ECMstorageNmob=zeros(timestep,number);
ECMfree=zeros(timestep,number);
ECMrsoilN=zeros(timestep,number);
ECMsoilN=zeros(timestep,number);
ECMcost_fixN=zeros(timestep,number);
ECMcost_resorb=zeros(timestep,number);
ECMcost_non_myc=zeros(timestep,number);
ECMcost_active=zeros(timestep,number);
ECMrec_cost_acq=zeros(timestep,number);
ECMcost_acq=zeros(timestep,number);
ECMNdeficit=zeros(timestep,number);
ECMNacq=zeros(timestep,number);
ECMCacq=zeros(timestep,number);
ECMCgrowth=zeros(timestep,number);
ECMCavailable=zeros(timestep,number);
ECMtotal_Nacq_active=zeros(timestep,number);
ECMtotal_Nacq_non_myc=zeros(timestep,number);
ECMtotal_Nacq_resorb=zeros(timestep,number);
ECMtotal_Nacq_fix=zeros(timestep,number);
ECMtotal_Nacq=zeros(timestep,number);
ECMtotal_Nacq_storage=zeros(timestep,number);
ECMtotal_Cgrowth=zeros(timestep,number);
ECMlitter_productionN=zeros(timestep,number);
ECMlitter_productionCN=zeros(timestep,number);
ECMrhizoCflux=zeros(timestep,number);
ECMfungalprod=zeros(timestep,number);

total_soil_N_uptake=zeros(timestep,number);
total_mineralN=zeros(timestep,number);

%Initialize CORPSE stuff
%Model parameters
% global params
if ~exist('params','var')
    params.vmaxref=[1500,28,1000]; %Relative maximum enzymatic decomp rates
    params.Ea=[37e3,54e3,50e3];    % Activation energy
    params.kC=[0.01,0.01,0.01];    % Michaelis-Menton parameter
    params.max_immobilization_rate=1.0; %Fraction of mineral N per day
    params.gas_diffusion_exp=2.5;  % Determines suppression of decomp at high soil moisture
    params.minMicrobeC=5e-4;       %Minimum microbial biomass (fraction of total C)
    params.Tmic=0.25;       % Microbial turnover rate
    params.et=0.6;          % Fraction of turnover not converted to CO2
    params.eup=[0.6,0.1,0.6]; % Carbon uptake efficiency
    params.nup=[0.3,0.3,0.3];  % Nitrogen uptake efficiency
    params.CN_microb=10.0;      % Microbe C:N ratio
    params.tProtected=100.0;    % Protected C turnover time (years)
    params.protection_rate=[0.6,0.0,4.0]; % Protected carbon formation rate (year-1)
    params.frac_turnover_min=0.2;  % Fraction of microbe turnover N that is mineralized
    params.frac_turnover_slow=0.3; % Fraction of microbe turnover that goes to slow pool

    params_AM=params;
    params_ECM=params;
end

params_litter_AM=params_AM;
params_litter_AM.protection_rate=[0,0,0];
params_litter_ECM=params_ECM;
params_litter_ECM.protection_rate=[0,0,0];

if ~exist('do_exudation','var')
    do_exudation=true;
end
if ~exist('do_litter','var')
    do_litter=true;
end

if ~exist('mineral_N_lost_fraction','var')
    mineral_N_lost_fraction=0.5;
end


if ~exist('plots_to_run','var')
    plots_to_run=1:number;
end


% Use the same litter inputs for all plots
if ~exist('constant_litter_inputs','var')
    constant_litter_inputs=false;
end

if ~exist('litter_transfer_to_soil','var')
    litter_transfer_to_soil=0.0;
end

if ~exist('initial_litter_fraction','var')
    initial_litter_fraction=0.3;
end

% This will work better if rhizo percent is greater than zero
ECM_rhizo_percent(ECM_rhizo_percent<0.1)=0.1;
AM_rhizo_percent(AM_rhizo_percent<0.1)=0.1;


%Set up soil data container and initial values for each plot
clear rhizosphere ECM_rhizospheres
clear bulksoil ECM_bulksoils
clear litterlayer ECM_litterlayers

if exist('use_restart_file','var')
    fprintf('Reading from restart file %s\n',use_restart_file)
    restart_data=load(use_restart_file);
    ECM_rhizospheres=restart_data.ECM_rhizospheres;
    ECM_bulksoils=restart_data.ECM_bulksoil;
    ECM_litterlayers=restart_data.ECM_litterlayers;
    mineralN=restart_data.ECM_mineral_N;
    ECMstorageN(1,:)=restart_data.ECMstorageN;
else

    rhiz_frac=ECM_rhizo_percent(1,:)/100;

    rhizosphere.litterC=repmat([5.6,1020+500,3.1],[number,1])*1e-3; %kgC of [fast,slow,dead microbe]
    rhizosphere.litterN=repmat([0.1,25.8+65,0.35],[number,1])*1e-3; %kgN of [fast,slow,dead microbe]
    rhizosphere.protectedC=repmat([140,0,600]*2,[number,1])*1e-3;
    rhizosphere.protectedN=repmat([13.5,0,65]*2,[number,1])*1e-3;
    rhizosphere.livingMicrobeC=zeros(number,1)+14.2e-3;
    rhizosphere.livingMicrobeN=rhizosphere.livingMicrobeC/params_ECM.CN_microb;
    rhizosphere.CO2=zeros(number,1);
    rhizosphere.Rtot=zeros(number,1);
    %These are for conservation checks
    rhizosphere.originalLitterC=(sum(rhizosphere.litterC+rhizosphere.protectedC,2)+rhizosphere.livingMicrobeC);
    rhizosphere.originalLitterN=(sum(rhizosphere.litterN+rhizosphere.protectedN,2)+rhizosphere.livingMicrobeN);

    %Initial fix: Just multiply soil pools by per_ECM so total pools can be
    %added together at the end.

    % So the 0% ECM plots don't start out completely empty
    pECM=per_ECM(:,2);
    pECM(pECM<0.1)=0.1;
    bulksoil=multiply_cohorts(rhizosphere,(1.0-rhiz_frac').*pECM/100*(1-initial_litter_fraction));
    rhizosphere=multiply_cohorts(rhizosphere,rhiz_frac'.*pECM/100*(1-initial_litter_fraction));
    litterlayer=multiply_cohort(rhizosphere,initial_litter_fraction);

    litterlayer.litterC=litterlayer.litterC+litterlayer.protectedC;
    litterlayer.litterN=litterlayer.litterN+litterlayer.protectedN;
    litterlayer.protectedC=0.0;
    litterlayer.protectedN=0.0;

    ECM_rhizospheres=rhizosphere;
    ECM_bulksoils=bulksoil;
    ECM_litterlayers=litterlayer;

    mineralN=zeros(1,number);
end

ECM_rhiz_outputs=setup_CORPSE_outputs(timestep,number);
ECM_bulk_outputs=setup_CORPSE_outputs(timestep,number);
ECM_litter_outputs=setup_CORPSE_outputs(timestep,number);




%(5) INITIALIZING MATRICES
AMpassive=zeros(timestep,number);
AMstorageNmob=zeros(timestep,number);
AMfree=zeros(timestep,number);
AMrsoilN=zeros(timestep,number);
AMsoilN=zeros(timestep,number);
AMcost_fixN=zeros(timestep,number);
AMcost_resorb=zeros(timestep,number);
AMcost_non_myc=zeros(timestep,number);
AMcost_active=zeros(timestep,number);
AMrec_cost_acq=zeros(timestep,number);
AMcost_acq=zeros(timestep,number);
AMNdeficit=zeros(timestep,number);
AMNacq=zeros(timestep,number);
AMCacq=zeros(timestep,number);
AMCgrowth=zeros(timestep,number);
AMCavailable=zeros(timestep,number);
AMtotal_Nacq_active=zeros(timestep,number);
AMtotal_Nacq_non_myc=zeros(timestep,number);
AMtotal_Nacq_resorb=zeros(timestep,number);
AMtotal_Nacq_fix=zeros(timestep,number);
AMtotal_Nacq=zeros(timestep,number);
AMtotal_Nacq_storage=zeros(timestep,number);
AMtotal_Cgrowth=zeros(timestep,number);
AMlitter_productionN=zeros(timestep,number);
AMlitter_productionCN=zeros(timestep,number);
AMrhizoCflux=zeros(timestep,number);
AMfungalprod=zeros(timestep,number);

%Set up soil data container and initial values for each plot
clear rhizosphere AM_rhizospheres
clear bulksoil AM_bulksoils

if exist('use_restart_file','var')
    fprintf('Reading from restart file %s\n',use_restart_file)
    restart_data=load(use_restart_file);
    AM_rhizospheres=restart_data.AM_rhizospheres;
    AM_bulksoils=restart_data.AM_bulksoil;
    AM_litterlayers=restart_data.AM_litterlayers;
    mineralN=restart_data.AM_mineral_N;
    AMstorageN(1,:)=restart_data.AMstorageN;
else

    rhiz_frac=AM_rhizo_percent(1,:)/100;

    rhizosphere.litterC=repmat([5.6,1020+500,3.1],[number,1])*1e-3; %kgC of [fast,slow,dead microbe]
    rhizosphere.litterN=repmat([0.1,25.8+65,0.35],[number,1])*1e-3; %kgN of [fast,slow,dead microbe]
    rhizosphere.protectedC=repmat([140,0,600]*2,[number,1])*1e-3;
    rhizosphere.protectedN=repmat([13.5,0,65]*2,[number,1])*1e-3;
    rhizosphere.livingMicrobeC=zeros(number,1)+14.2e-3;
    rhizosphere.livingMicrobeN=rhizosphere.livingMicrobeC/params_AM.CN_microb;
    rhizosphere.CO2=zeros(number,1);
    rhizosphere.Rtot=zeros(number,1);
    %These are for conservation checks
    rhizosphere.originalLitterC=(sum(rhizosphere.litterC+rhizosphere.protectedC,2)+rhizosphere.livingMicrobeC);
    rhizosphere.originalLitterN=(sum(rhizosphere.litterN+rhizosphere.protectedN,2)+rhizosphere.livingMicrobeN);

    %So the 0% AM plots don't start out completely empty
    pAM=100-per_ECM(:,2);
    pAM(pAM<0.1)=0.1;
    bulksoil=multiply_cohorts(rhizosphere,(1.0-rhiz_frac').*(pAM/100)*(1-initial_litter_fraction));
    rhizosphere=multiply_cohorts(rhizosphere,rhiz_frac'.*(pAM/100)*(1-initial_litter_fraction));
    litterlayer=multiply_cohort(rhizosphere,initial_litter_fraction);

    litterlayer.litterC=litterlayer.litterC+litterlayer.protectedC;
    litterlayer.litterN=litterlayer.litterN+litterlayer.protectedN;
    litterlayer.protectedC=0.0;
    litterlayer.protectedN=0.0;

    AM_rhizospheres=rhizosphere;
    AM_bulksoils=bulksoil;
    AM_litterlayers=litterlayer;

    mineralN=zeros(1,number);

end

AM_rhiz_outputs=setup_CORPSE_outputs(timestep,number);
AM_bulk_outputs=setup_CORPSE_outputs(timestep,number);
AM_litter_outputs=setup_CORPSE_outputs(timestep,number);


if exist('theta_override','var')
    soilVWC=theta_override;
end

soilVWC_orig=soilVWC;
soilVWC_factor=1.0;

%%

%(6)  MAIN LOOP

for i2=1:timestep
    if mod(i2,365)==0
        fprintf('Starting year %d\n',floor(i2/365))
    end
    i=mod(i2-1,365)+1;



   %Update CORPSE soil state first
    [rhizosphere,orhiz]=update_cohorts(ECM_rhizospheres,soilT(i)+273.15,soilVWC(i),mineralN,1.0/365.0,params_ECM);
    mineralN(:)=orhiz.mineralN;
    check_cohorts(rhizosphere);

    %ECM_rhiz_outputs=record_CORPSE_outputs(ECM_rhiz_outputs,rhizosphere,orhiz,i2,j);
    ECM_rhiz_outputs.unprotectedC.fast(i2,:)=rhizosphere.litterC(:,1);
    ECM_rhiz_outputs.unprotectedC.slow(i2,:)=rhizosphere.litterC(:,2);
    ECM_rhiz_outputs.unprotectedC.deadmic(i2,:)=rhizosphere.litterC(:,3);
    ECM_rhiz_outputs.protectedC.fast(i2,:)=rhizosphere.protectedC(:,1);
    ECM_rhiz_outputs.protectedC.slow(i2,:)=rhizosphere.protectedC(:,2);
    ECM_rhiz_outputs.protectedC.deadmic(i2,:)=rhizosphere.protectedC(:,3);
    ECM_rhiz_outputs.unprotectedN.fast(i2,:)=rhizosphere.litterN(:,1);
    ECM_rhiz_outputs.unprotectedN.slow(i2,:)=rhizosphere.litterN(:,2);
    ECM_rhiz_outputs.unprotectedN.deadmic(i2,:)=rhizosphere.litterN(:,3);
    ECM_rhiz_outputs.protectedN.fast(i2,:)=rhizosphere.protectedN(:,1);
    ECM_rhiz_outputs.protectedN.slow(i2,:)=rhizosphere.protectedN(:,2);
    ECM_rhiz_outputs.protectedN.deadmic(i2,:)=rhizosphere.protectedN(:,3);
    ECM_rhiz_outputs.N_mineralization(i2,:)=orhiz.N_mineralization;
    ECM_rhiz_outputs.N_immobilization(i2,:)=orhiz.N_immobilization;
    ECM_rhiz_outputs.livingMicrobeC(i2,:)=rhizosphere.livingMicrobeC;

    ECM_rhiz_outputs.decomp.fast(i2,:)=orhiz.decomp(:,1);
    ECM_rhiz_outputs.decomp.slow(i2,:)=orhiz.decomp(:,2);
    ECM_rhiz_outputs.decomp.deadmic(i2,:)=orhiz.decomp(:,3);
    ECM_rhiz_outputs.CO2prod(i2,:)=orhiz.CO2prod;

    ECM_rhiz_outputs.N_decomposed.fast(i2,:)=orhiz.totalN_decomposed(:,1);
    ECM_rhiz_outputs.N_decomposed.slow(i2,:)=orhiz.totalN_decomposed(:,2);
    ECM_rhiz_outputs.N_decomposed.deadmic(i2,:)=orhiz.totalN_decomposed(:,3);



    [bulksoil,obulk]=update_cohorts(ECM_bulksoils,soilT(i)+273.15,soilVWC(i),mineralN,1.0/365.0,params_ECM);
    mineralN=obulk.mineralN;
    check_cohorts(bulksoil);

    ECM_bulk_outputs.unprotectedC.fast(i2,:)=bulksoil.litterC(:,1);
    ECM_bulk_outputs.unprotectedC.slow(i2,:)=bulksoil.litterC(:,2);
    ECM_bulk_outputs.unprotectedC.deadmic(i2,:)=bulksoil.litterC(:,3);
    ECM_bulk_outputs.protectedC.fast(i2,:)=bulksoil.protectedC(:,1);
    ECM_bulk_outputs.protectedC.slow(i2,:)=bulksoil.protectedC(:,2);
    ECM_bulk_outputs.protectedC.deadmic(i2,:)=bulksoil.protectedC(:,3);
    ECM_bulk_outputs.unprotectedN.fast(i2,:)=bulksoil.litterN(:,1);
    ECM_bulk_outputs.unprotectedN.slow(i2,:)=bulksoil.litterN(:,2);
    ECM_bulk_outputs.unprotectedN.deadmic(i2,:)=bulksoil.litterN(:,3);
    ECM_bulk_outputs.protectedN.fast(i2,:)=bulksoil.protectedN(:,1);
    ECM_bulk_outputs.protectedN.slow(i2,:)=bulksoil.protectedN(:,2);
    ECM_bulk_outputs.protectedN.deadmic(i2,:)=bulksoil.protectedN(:,3);
    ECM_bulk_outputs.N_mineralization(i2,:)=obulk.N_mineralization;
    ECM_bulk_outputs.N_immobilization(i2,:)=obulk.N_immobilization;
    ECM_bulk_outputs.livingMicrobeC(i2,:)=bulksoil.livingMicrobeC;

    ECM_bulk_outputs.decomp.fast(i2,:)=obulk.decomp(:,1);
    ECM_bulk_outputs.decomp.slow(i2,:)=obulk.decomp(:,2);
    ECM_bulk_outputs.decomp.deadmic(i2,:)=obulk.decomp(:,3);
    ECM_bulk_outputs.CO2prod(i2,:)=obulk.CO2prod;

    ECM_bulk_outputs.N_decomposed.fast(i2,:)=obulk.totalN_decomposed(:,1);
    ECM_bulk_outputs.N_decomposed.slow(i2,:)=obulk.totalN_decomposed(:,2);
    ECM_bulk_outputs.N_decomposed.deadmic(i2,:)=obulk.totalN_decomposed(:,3);


    [litterlayer,olitter]=update_cohorts(ECM_litterlayers,soilT(i)+273.15,soilVWC(i),mineralN,1.0/365.0,params_litter_ECM);
    mineralN=olitter.mineralN;
    check_cohorts(litterlayer);

    %ECM_litter_outputs=record_CORPSE_outputs(ECM_litter_outputs,litterlayer,obulk,i2,j);
    ECM_litter_outputs.unprotectedC.fast(i2,:)=litterlayer.litterC(:,1);
    ECM_litter_outputs.unprotectedC.slow(i2,:)=litterlayer.litterC(:,2);
    ECM_litter_outputs.unprotectedC.deadmic(i2,:)=litterlayer.litterC(:,3);
    ECM_litter_outputs.protectedC.fast(i2,:)=litterlayer.protectedC(:,1);
    ECM_litter_outputs.protectedC.slow(i2,:)=litterlayer.protectedC(:,2);
    ECM_litter_outputs.protectedC.deadmic(i2,:)=litterlayer.protectedC(:,3);
    ECM_litter_outputs.unprotectedN.fast(i2,:)=litterlayer.litterN(:,1);
    ECM_litter_outputs.unprotectedN.slow(i2,:)=litterlayer.litterN(:,2);
    ECM_litter_outputs.unprotectedN.deadmic(i2,:)=litterlayer.litterN(:,3);
    ECM_litter_outputs.protectedN.fast(i2,:)=litterlayer.protectedN(:,1);
    ECM_litter_outputs.protectedN.slow(i2,:)=litterlayer.protectedN(:,2);
    ECM_litter_outputs.protectedN.deadmic(i2,:)=litterlayer.protectedN(:,3);
    ECM_litter_outputs.N_mineralization(i2,:)=olitter.N_mineralization;
    ECM_litter_outputs.N_immobilization(i2,:)=olitter.N_immobilization;
    ECM_litter_outputs.livingMicrobeC(i2,:)=litterlayer.livingMicrobeC;

    ECM_litter_outputs.decomp.fast(i2,:)=olitter.decomp(:,1);
    ECM_litter_outputs.decomp.slow(i2,:)=olitter.decomp(:,2);
    ECM_litter_outputs.decomp.deadmic(i2,:)=olitter.decomp(:,3);
    ECM_litter_outputs.CO2prod(i2,:)=olitter.CO2prod;

    ECM_litter_outputs.N_decomposed.fast(i2,:)=olitter.totalN_decomposed(:,1);
    ECM_litter_outputs.N_decomposed.slow(i2,:)=olitter.totalN_decomposed(:,2);
    ECM_litter_outputs.N_decomposed.deadmic(i2,:)=olitter.totalN_decomposed(:,3);


    %Change in rhizosphere size
    rhiz_frac=ECM_rhizo_percent(i,:)/100;
    current_rhiz=rhiz_frac;
    if i+1<=365
        next_rhiz=ECM_rhizo_percent(i+1,:)/100;
        rhiz_diff=next_rhiz-current_rhiz;


        rhiz_diff=rhiz_diff';

        newrhiz=multiply_cohorts(bulksoil,max(rhiz_diff,0));
        rhizosphere=add_cohorts(rhizosphere,newrhiz);
        bulksoil=add_cohorts(bulksoil,multiply_cohort(newrhiz,-1.0));

        newbulk=multiply_cohorts(rhizosphere,-min(rhiz_diff,0));
        bulksoil=add_cohorts(bulksoil,newbulk);
        rhizosphere=add_cohorts(rhizosphere,multiply_cohort(newbulk,-1.0));

    else
        % Mix soil fractions at end of year
        totalsoil=add_cohorts(rhizosphere,bulksoil);
        rhizosphere=multiply_cohorts(totalsoil,rhiz_frac');
        bulksoil=multiply_cohorts(totalsoil,1.0-rhiz_frac');
    end



    % Transfer fraction of litter to soil
    newsoil=multiply_cohort(litterlayer,litter_transfer_to_soil);
    rhizosphere=add_cohorts(rhizosphere,multiply_cohorts(newsoil,rhiz_frac'));
    bulksoil=add_cohorts(bulksoil,multiply_cohorts(newsoil,1-rhiz_frac'));
    litterlayer=multiply_cohort(litterlayer,1.0-litter_transfer_to_soil);


    check_cohorts(rhizosphere);check_cohorts(bulksoil);check_cohorts(litterlayer);

    ECM_rhizospheres=rhizosphere;
    ECM_bulksoils=bulksoil;
    ECM_litterlayers=litterlayer;


    %%%%%%%%%%%%

    % %(6)  RUNNING AM PORTION OF THE FUN CODE

    %Update CORPSE soil state first
    [rhizosphere,orhiz]=update_cohorts(AM_rhizospheres,soilT(i)+273.15,soilVWC(i),mineralN(:),1.0/365.0,params_AM);
    mineralN(:)=orhiz.mineralN;
    check_cohorts(rhizosphere);

    %AM_rhiz_outputs=record_CORPSE_outputs(AM_rhiz_outputs,rhizosphere,orhiz,i2,j);
    AM_rhiz_outputs.unprotectedC.fast(i2,:)=rhizosphere.litterC(:,1);
    AM_rhiz_outputs.unprotectedC.slow(i2,:)=rhizosphere.litterC(:,2);
    AM_rhiz_outputs.unprotectedC.deadmic(i2,:)=rhizosphere.litterC(:,3);
    AM_rhiz_outputs.protectedC.fast(i2,:)=rhizosphere.protectedC(:,1);
    AM_rhiz_outputs.protectedC.slow(i2,:)=rhizosphere.protectedC(:,2);
    AM_rhiz_outputs.protectedC.deadmic(i2,:)=rhizosphere.protectedC(:,3);
    AM_rhiz_outputs.unprotectedN.fast(i2,:)=rhizosphere.litterN(:,1);
    AM_rhiz_outputs.unprotectedN.slow(i2,:)=rhizosphere.litterN(:,2);
    AM_rhiz_outputs.unprotectedN.deadmic(i2,:)=rhizosphere.litterN(:,3);
    AM_rhiz_outputs.protectedN.fast(i2,:)=rhizosphere.protectedN(:,1);
    AM_rhiz_outputs.protectedN.slow(i2,:)=rhizosphere.protectedN(:,2);
    AM_rhiz_outputs.protectedN.deadmic(i2,:)=rhizosphere.protectedN(:,3);
    AM_rhiz_outputs.N_mineralization(i2,:)=orhiz.N_mineralization;
    AM_rhiz_outputs.N_immobilization(i2,:)=orhiz.N_immobilization;
    AM_rhiz_outputs.livingMicrobeC(i2,:)=rhizosphere.livingMicrobeC;

    AM_rhiz_outputs.decomp.fast(i2,:)=orhiz.decomp(:,1);
    AM_rhiz_outputs.decomp.slow(i2,:)=orhiz.decomp(:,2);
    AM_rhiz_outputs.decomp.deadmic(i2,:)=orhiz.decomp(:,3);
    AM_rhiz_outputs.CO2prod(i2,:)=orhiz.CO2prod;

    AM_rhiz_outputs.N_decomposed.fast(i2,:)=orhiz.totalN_decomposed(:,1);
    AM_rhiz_outputs.N_decomposed.slow(i2,:)=orhiz.totalN_decomposed(:,2);
    AM_rhiz_outputs.N_decomposed.deadmic(i2,:)=orhiz.totalN_decomposed(:,3);

    [bulksoil,obulk]=update_cohorts(AM_bulksoils,soilT(i)+273.15,soilVWC(i),mineralN,1.0/365.0,params_AM);
    mineralN(:)=obulk.mineralN;
    check_cohorts(bulksoil);

    %AM_bulk_outputs=record_CORPSE_outputs(AM_bulk_outputs,bulksoil,obulk,i2,j);
    AM_bulk_outputs.unprotectedC.fast(i2,:)=bulksoil.litterC(:,1);
    AM_bulk_outputs.unprotectedC.slow(i2,:)=bulksoil.litterC(:,2);
    AM_bulk_outputs.unprotectedC.deadmic(i2,:)=bulksoil.litterC(:,3);
    AM_bulk_outputs.protectedC.fast(i2,:)=bulksoil.protectedC(:,1);
    AM_bulk_outputs.protectedC.slow(i2,:)=bulksoil.protectedC(:,2);
    AM_bulk_outputs.protectedC.deadmic(i2,:)=bulksoil.protectedC(:,3);
    AM_bulk_outputs.unprotectedN.fast(i2,:)=bulksoil.litterN(:,1);
    AM_bulk_outputs.unprotectedN.slow(i2,:)=bulksoil.litterN(:,2);
    AM_bulk_outputs.unprotectedN.deadmic(i2,:)=bulksoil.litterN(:,3);
    AM_bulk_outputs.protectedN.fast(i2,:)=bulksoil.protectedN(:,1);
    AM_bulk_outputs.protectedN.slow(i2,:)=bulksoil.protectedN(:,2);
    AM_bulk_outputs.protectedN.deadmic(i2,:)=bulksoil.protectedN(:,3);
    AM_bulk_outputs.N_mineralization(i2,:)=obulk.N_mineralization;
    AM_bulk_outputs.N_immobilization(i2,:)=obulk.N_immobilization;
    AM_bulk_outputs.livingMicrobeC(i2,:)=bulksoil.livingMicrobeC;

    AM_bulk_outputs.decomp.fast(i2,:)=obulk.decomp(:,1);
    AM_bulk_outputs.decomp.slow(i2,:)=obulk.decomp(:,2);
    AM_bulk_outputs.decomp.deadmic(i2,:)=obulk.decomp(:,3);
    AM_bulk_outputs.CO2prod(i2,:)=obulk.CO2prod;

    AM_bulk_outputs.N_decomposed.fast(i2,:)=obulk.totalN_decomposed(:,1);
    AM_bulk_outputs.N_decomposed.slow(i2,:)=obulk.totalN_decomposed(:,2);
    AM_bulk_outputs.N_decomposed.deadmic(i2,:)=obulk.totalN_decomposed(:,3);


    [litterlayer,olitter]=update_cohorts(AM_litterlayers,soilT(i)+273.15,soilVWC(i),mineralN,1.0/365.0,params_litter_AM);
    mineralN(:)=olitter.mineralN;
    check_cohorts(litterlayer);

    %AM_litter_outputs=record_CORPSE_outputs(AM_litter_outputs,litterlayer,obulk,i2,j);
    AM_litter_outputs.unprotectedC.fast(i2,:)=litterlayer.litterC(:,1);
    AM_litter_outputs.unprotectedC.slow(i2,:)=litterlayer.litterC(:,2);
    AM_litter_outputs.unprotectedC.deadmic(i2,:)=litterlayer.litterC(:,3);
    AM_litter_outputs.protectedC.fast(i2,:)=litterlayer.protectedC(:,1);
    AM_litter_outputs.protectedC.slow(i2,:)=litterlayer.protectedC(:,2);
    AM_litter_outputs.protectedC.deadmic(i2,:)=litterlayer.protectedC(:,3);
    AM_litter_outputs.unprotectedN.fast(i2,:)=litterlayer.litterN(:,1);
    AM_litter_outputs.unprotectedN.slow(i2,:)=litterlayer.litterN(:,2);
    AM_litter_outputs.unprotectedN.deadmic(i2,:)=litterlayer.litterN(:,3);
    AM_litter_outputs.protectedN.fast(i2,:)=litterlayer.protectedN(:,1);
    AM_litter_outputs.protectedN.slow(i2,:)=litterlayer.protectedN(:,2);
    AM_litter_outputs.protectedN.deadmic(i2,:)=litterlayer.protectedN(:,3);
    AM_litter_outputs.N_mineralization(i2,:)=olitter.N_mineralization;
    AM_litter_outputs.N_immobilization(i2,:)=olitter.N_immobilization;
    AM_litter_outputs.livingMicrobeC(i2,:)=litterlayer.livingMicrobeC;

    AM_litter_outputs.decomp.fast(i2,:)=olitter.decomp(:,1);
    AM_litter_outputs.decomp.slow(i2,:)=olitter.decomp(:,2);
    AM_litter_outputs.decomp.deadmic(i2,:)=olitter.decomp(:,3);
    AM_litter_outputs.CO2prod(i2,:)=olitter.CO2prod;

    AM_litter_outputs.N_decomposed.fast(i2,:)=olitter.totalN_decomposed(:,1);
    AM_litter_outputs.N_decomposed.slow(i2,:)=olitter.totalN_decomposed(:,2);
    AM_litter_outputs.N_decomposed.deadmic(i2,:)=olitter.totalN_decomposed(:,3);



    %Change in rhizosphere size
    rhiz_frac=AM_rhizo_percent(i,:)/100;
    current_rhiz=rhiz_frac;
    if i+1<=365
        next_rhiz=AM_rhizo_percent(i+1,:)/100;
        rhiz_diff=next_rhiz-current_rhiz;

        rhiz_diff=rhiz_diff';
        newrhiz=multiply_cohorts(bulksoil,max(rhiz_diff,0));
        rhizosphere=add_cohorts(rhizosphere,newrhiz);
        bulksoil=add_cohorts(bulksoil,multiply_cohort(newrhiz,-1.0));

        newbulk=multiply_cohorts(rhizosphere,-min(rhiz_diff,0));
        bulksoil=add_cohorts(bulksoil,newbulk);
        rhizosphere=add_cohorts(rhizosphere,multiply_cohort(newbulk,-1.0));

    else
        % Mix soil fractions at end of year
        totalsoil=add_cohorts(rhizosphere,bulksoil);
        rhizosphere=multiply_cohorts(totalsoil,rhiz_frac');
        bulksoil=multiply_cohorts(totalsoil,1.0-rhiz_frac');
    end



     % Transfer fraction of litter to soil
    newsoil=multiply_cohort(litterlayer,litter_transfer_to_soil);
    rhizosphere=add_cohorts(rhizosphere,multiply_cohorts(newsoil,rhiz_frac'));
    bulksoil=add_cohorts(bulksoil,multiply_cohorts(newsoil,1-rhiz_frac'));
    litterlayer=multiply_cohort(litterlayer,1.0-litter_transfer_to_soil);


    check_cohorts(rhizosphere);check_cohorts(bulksoil);check_cohorts(litterlayer);

    AM_rhizospheres=rhizosphere;
    AM_bulksoils=bulksoil;
    AM_litterlayers=litterlayer;


    % Update inorganic N pool
    ECMsoilN(i2,:)=mineralN(:);
    AMsoilN(i2,:)=mineralN(:);

    % Start calculating FUN costs
    % Assume passive uptake (with transpiration) is zero
    ECMpassive(i2,:)=0;
    %Making this zero for now because of ET and SD values
    ECMpassive(i2,:)=min(ECMpassive(i2,:),ECMsoilN(i2,:));

    ECMrsoilN(i2,:)=ECMsoilN(i2,:)-ECMpassive(i2,:); %Soil N remaining

    %Sum Passive and Storage N allocation
    ECMfree(i2,:)=ECMpassive(i2,:);%+ECMstorageNmob(i2,j);

    %%%%Calculates the cost for each uptake strategy%%%%

    ECMcost_fixN(i2,:)=s_FIX*exp(a_FIX+b_FIX*soilT(i)*(1-0.5*soilT(i)/c_FIX))-2*s_FIX; %Cost fixation
    ECMcost_resorb(i2,:)=inf;%Cost retranslocation fixed to be high during green season
    ECMsoilN(i2,:)=ECMrsoilN(i2,:);
    ECMcost_non_myc(i2,:)=(NkN./ECMsoilN(i2,:))+(NkC./(AMrb(i,:)+ECMrb(i,:))); %Cost of non-mycorrhizal uptake
    ECMcost_active(i2,:)=(EkN./ECMsoilN(i2,:))+(EkC./(AMrb(i,:)+ECMrb(i,:))); %Cost of mycorrhizal uptake]

    ECMcost_storage=zeros(1,number)+2.0e-3;

    ECMcost_storage(ECMstorageN(i2,:)==0)=inf;

    %%%%%Adding in a newer retranslocation function that calculate at a
    %%%%%annual scale and the spreads it out over the senescence period.
    if i==303
        ECMleafN2=sum(ECMleafN);
        ECMcost_resorb(i2,:)=kR./ECMleafN2(1,:);
        clear ECMleafN2
    end


    %Use resistance network to get integrated cost
    ECMrec_cost_acq(i2,:)=(1./ECMcost_active(i2,:))+(1./ECMcost_fixN(i2,:))+(1./ECMcost_resorb(i2,:))+(1./ECMcost_non_myc(i2,:))+(1./ECMcost_storage);
    ECMcost_acq(i2,:)=1./ECMrec_cost_acq(i2,:);
    ECMleafN2=sum(ECMleafN);
    if i==303
        ECMcost_resorb(i2,:)=kR.*per_ECM(:,2)'/100./ECMleafN2(1,:);
        %ECMcost_fixN(i2,j)=inf;  %Turn these off at retranslocation
        %ECMcost_non_myc(i2,j)=inf;
        %ECMcost_active(i2,j)=inf;
        ECMrec_cost_acq(i2,:)=(1./ECMcost_active(i2,:))+(1./ECMcost_fixN(i2,:))+(1./ECMcost_resorb(i2,:))+(1./ECMcost_non_myc(i2,:)+(1./ECMcost_storage));
        ECMcost_acq(i2,:)=1./ECMrec_cost_acq(i2,:);
    end

    %This code below simultaneously solves eq 6a-6d from the Fisher et al.
    %(2010) to optimize C allocated to growth and N uptake. N uptake is
    %truncated if Ndemand is met or leafN pool or soilN pool is
    %exhausted.

    if i==303
      Ndemand_plus_storage=ECMstorageNdemand(i,:)+ECMtotalNdemand(i,:);
      %Add in storage C into NPP or C available pool
      npp_plus_storage=ECMNPP(i,:)+ECMstorageC(1,:);
    else
       npp_plus_storage=ECMNPP(i,:);
       Ndemand_plus_storage=ECMtotalNdemand(i,:);
    end



    ECMplantCN(i2,:)=npp_plus_storage./Ndemand_plus_storage;
    ECMplantCN(i2,isnan(ECMplantCN(i2,:)))=0;

    ECMCavailable(i2,:)=npp_plus_storage;
    ECMNdeficit(i2,:)=Ndemand_plus_storage-ECMfree(i2,:);

    xx=(ECMNdeficit(i2,:)<=0);

    ECMNacq(i2,xx)=0;
    ECMCacq(i2,xx)=0;
    ECMCgrowth(i2,xx)=ECMCgrowth(i2,xx);

    xx=(ECMNdeficit(i2,:)>0);
    ECMCacq(i2,xx)=(ECMCavailable(i2,xx)-ECMfree(i2,xx).*ECMplantCN(i,xx))./(1+ECMplantCN(i,xx)./ECMcost_acq(i2,xx));
    ECMNacq(i2,xx)=min(ECMNdeficit(i2,xx),ECMCacq(i2,xx)./ECMcost_acq(i2,xx));
    ECMCgrowth(i2,xx)=ECMCavailable(i2,xx)-(ECMNacq(i2,xx).*ECMcost_acq(i2,xx));

    ECMtotal_Nacq_active(i2,xx)=ECMCacq(i2,xx)./ECMcost_active(i2,xx);
    ECMtotal_Nacq_non_myc(i2,xx)=ECMCacq(i2,xx)./ECMcost_non_myc(i2,xx);
    ECMtotal_Nacq_resorb(i2,xx)=ECMCacq(i2,xx)./ECMcost_resorb(i2,xx);
    ECMtotal_Nacq_fix(i2,xx)=ECMCacq(i2,xx)./ECMcost_fixN(i2,xx);
    ECMtotal_Nacq_storage(i2,xx)=ECMCacq(i2,xx)./ECMcost_storage(xx);


    if exist('prescribed_ECM_rhizoCflux','var')
        ECMrhizoCflux(i2,:)=prescribed_ECM_rhizoCflux(i2,:);
    else
        ECMrhizoCflux(i2,:)=(ECMtotal_Nacq_non_myc(i2,:)+ECMtotal_Nacq_fix(i2,:)).*ECMcost_acq(i2,:);
    end

    if exist('prescribed_ECM_fungal_prod','var')
        ECMfungalprod(i2,:)=prescribed_ECM_fungal_prod(i2,:);
    else
        ECMfungalprod(i2,:)=(ECMtotal_Nacq_active(i2,:)).*ECMcost_acq(i2,:);
    end


    %Add rhizo C flux to rhizosphere as fast C
    %Should different sources use their different costs?
    if do_exudation

        ECM_rhizospheres.litterC(:,1)=ECM_rhizospheres.litterC(:,1)+ECMrhizoCflux(i2,:)';

        rhiz_fungi_frac=ECM_rhizo_percent(i,:)'/100;
        ECM_rhizospheres.litterC(:,1)=ECM_rhizospheres.litterC(:,1)+ECMfungalprod(i2,:)'.*rhiz_fungi_frac.*(1.0-litter_fungi_frac);
        ECM_bulksoils.litterC(:,1)=ECM_bulksoils.litterC(:,1)+ECMfungalprod(i2,:)'.*(1.0-rhiz_fungi_frac).*(1.0-litter_fungi_frac);
        ECM_litterlayers.litterC(:,1)=ECM_litterlayers.litterC(:,1)+ECMfungalprod(i2,:)'*litter_fungi_frac;

        if exudate_CN>0
          ECM_rhizospheres.litterN(:,1)=ECM_rhizospheres.litterN(:,1)+ECMfungalprod(i2,:)'.*rhiz_fungi_frac.*(1.0-litter_fungi_frac)/exudate_CN;
          ECM_bulksoils.litterN(:,1)=ECM_bulksoils.litterN(:,1)+ECMfungalprod(i2,:)'.*(1.0-rhiz_fungi_frac).*(1.0-litter_fungi_frac)/exudate_CN;
          ECM_litterlayers.litterN(:,1)=ECM_litterlayers.litterN(:,1)+ECMfungalprod(i2,:)'*litter_fungi_frac/exudate_CN;
        end
    end


    if any(ECMsoilN(i2,:)<0)
        error('soilN<0')
    end

    ECMtotal_Nacq_storage(i2,:)=min(ECMtotal_Nacq_storage(i2,:),ECMstorageN(i2,:));

    xx=ECMtotal_Nacq_resorb(i2,:)>ECMleafN2(1,:); %Truncate resorb if exceeds leaf N

    ECMtotal_Nacq_resorb(i2,xx)=ECMleafN2(1,xx);
    ECMCgrowth(i2,xx)=ECMCavailable(i2,xx)-(ECMNacq(i2,xx).*ECMcost_acq(i2,xx));


    %Update storage N pool to include retranslocation
    if i<303
        ECMstorageN(i2+1,:)=ECMstorageN(i2,:)-ECMtotal_Nacq_storage(i2,:);
    else
        ECMstorageN(i2+1,:)=ECMstorageN(i2,:)+ECMtotal_Nacq_resorb(i2,:)+ECMtotal_Nacq_active(i2,:)+ECMtotal_Nacq_non_myc(i2,:)+ECMtotal_Nacq_fix(i2,:);

        xx=ECMleafN2(:)>0;

        ECMresorbper(1,xx)=1-(ECMtotal_Nacq_resorb(303,xx)./ECMleafN2(1,xx));
        ECMlitter_productionN(i2,xx)=(ECMleafN(i,xx).*ECMresorbper(1,xx));
        ECMlitter_productionCN(i2,xx)=ECMlitter_production(i,xx)./ECMlitter_productionN(i2,xx);

        xx=ECMleafN2(:)<=0;
        ECMresorbper(1,xx)=0.0;
        ECMlitter_productionN(i2,xx)=0.0;
        ECMlitter_productionCN(i2,xx)=0.0;

    end

    %Add in ability to get N from storage in leaf_out period
    if i<303
      ECMstorageNmob(i2,:)=min(ECMstorageN(i2,:),ECMtotal_Nacq_storage(i2,:));
      ECMstorageN(i2+1,:)=ECMstorageN(i2,:)-ECMstorageNmob(i2,:);%Update storage N ECMount
    end

    if do_litter
        %Add litter to soil pools.
        if constant_litter_inputs
            litter_N_input=meanlitter_production(i)./repmat(ECMlitter_productionCN(i2,:),[3,1]).*repmat(ECM_litter_composition,[number,1])';
            litter_N_input(~isfinite(litter_N_input))=0.0;

        else
            litter_N_input=repmat(ECMlitter_productionN(i2,:),[3,1]).*repmat(ECM_litter_composition,[number,1])';
        end

        ECM_litterlayers.litterN=ECM_litterlayers.litterN+litter_N_input';

        if constant_litter_inputs
            litter_C_input=meanlitter_production(i).*repmat(ECM_litter_composition,[number,1])';
        else
            litter_C_input=repmat(ECMlitter_production(i,:),[3,1]).*repmat(ECM_litter_composition,[number,1])';
        end
        % ECM_rhizospheres(j).litterC=ECM_rhizospheres(j).litterC+litter_C_input*rhiz_frac;
        % ECM_bulksoils(j).litterC=ECM_bulksoils(j).litterC+litter_C_input*(1.0-rhiz_frac);
        ECM_litterlayers.litterC=ECM_litterlayers.litterC+litter_C_input';

        %Root litter
        if constant_litter_inputs
            rootlitter_N_input=repmat(mean(daily_froot_turnover(i,:))/1000.*(per_ECM(:,2)'/100)./ECMrootCN,[3,1]).*repmat(ECM_rootlitter_composition,[number,1])';
        else
            rootlitter_N_input=repmat(daily_froot_turnover(i,:)/1000.*(per_ECM(:,2)'/100)./ECMrootCN,[3,1]).*repmat(ECM_rootlitter_composition,[number,1])';
        end


        if constant_litter_inputs
            rootlitter_C_input=repmat(mean(daily_froot_turnover(i,:))/1000.*(per_ECM(:,2)'/100),[3,1]).*repmat(ECM_rootlitter_composition,[number,1])';
        else
            rootlitter_C_input=repmat(daily_froot_turnover(i,:)/1000.*(per_ECM(:,2)'/100),[3,1]).*repmat(ECM_rootlitter_composition,[number,1])';
        end

        if root_inputs_only_rhiz
            ECM_rhizospheres.litterC=ECM_rhizospheres.litterC+rootlitter_C_input';
            ECM_rhizospheres.litterN=ECM_rhizospheres.litterN+rootlitter_N_input';
        else
            rf=repmat(ECM_rhizo_percent(i,:)/100,[3,1]);
            ECM_rhizospheres.litterC=ECM_rhizospheres.litterC+rootlitter_C_input'.*rf';
            ECM_bulksoils.litterC=ECM_bulksoils.litterC+rootlitter_C_input'.*(1.0-rf)';
            ECM_rhizospheres.litterN=ECM_rhizospheres.litterN+rootlitter_N_input'.*rf';
            ECM_bulksoils.litterN=ECM_bulksoils.litterN+rootlitter_N_input'.*(1.0-rf)';
        end

      end

      check_cohorts(ECM_rhizospheres);

      ECMtotal_Nacq(i2,:)=ECMNacq(i2,:)+ECMfree(i2,:);
      ECMtotal_Cgrowth(i2,:)=sum(ECMCgrowth(i2,:));


    %%%%%%%%%%%%

    % AM plots
    AMpassive(i2,:)=0;
    AMpassive(i2,:)=min(AMpassive(i2,:),AMsoilN(i2,:));
    AMrsoilN(i2,:)=AMsoilN(i2,:)-AMpassive(i2,:); %Soil N remaining

    %Sum Passive and Storage N allocation
    AMfree(i2,:)=AMpassive(i2,:);%+AMstorageNmob(i2,j);

    %%%%Calculates the cost for each uptake strategy%%%%

    AMcost_fixN(i2,:)=s_FIX*exp(a_FIX+b_FIX*soilT(i)*(1-0.5*soilT(i)/c_FIX))-2*s_FIX; %Cost fixation
    AMcost_resorb(i2,:)=inf;%Cost retranslocation fixed to be high during green season
    AMsoilN(i2,:)=AMrsoilN(i2,:);
    AMcost_non_myc(i2,:)=(NkN./AMsoilN(i2,:))+(NkC./(ECMrb(i,:)+AMrb(i,:))); %Cost of non-mycorrhizal uptake
    AMcost_active(i2,:)=(AkN./AMsoilN(i2,:))+(AkC./(ECMrb(i,:)+AMrb(i,:))); %Cost of mycorrhizal uptake]

    AMcost_storage=zeros(1,number)+2.0e-3;

    AMcost_storage(AMstorageN(i2,:)==0)=inf;

    %%%%%Adding in a newer retranslocation function that calculate at a
    %%%%%annual scale and the spreads it out over the senescence period.
    if i==303
        AMleafN2=sum(AMleafN);
        AMcost_resorb(i2,:)=kR./AMleafN2(1,:);
        clear AMleafN2
    end

    %Use resistance network to get integrated cost
    AMrec_cost_acq(i2,:)=(1./AMcost_active(i2,:))+(1./AMcost_fixN(i2,:))+(1./AMcost_resorb(i2,:))+(1./AMcost_non_myc(i2,:))+(1./AMcost_storage);
    AMcost_acq(i2,:)=1./AMrec_cost_acq(i2,:);
    AMleafN2=sum(AMleafN);
    if i==303
        AMcost_resorb(i2,:)=kR.*(1-per_ECM(:,2)'/100)./AMleafN2(1,:);
        AMrec_cost_acq(i2,:)=(1./AMcost_active(i2,:))+(1./AMcost_fixN(i2,:))+(1./AMcost_resorb(i2,:))+(1./AMcost_non_myc(i2,:)+(1./AMcost_storage));
        AMcost_acq(i2,:)=1./AMrec_cost_acq(i2,:);
    end

   %This code below simultaneously solves eq 6a-6d from the Fisher et al.
   %(2010) to optimize C allocated to growth and N uptake. N uptake is
   %truncated if Ndemand is met or leafN pool or soilN pool is
   %exhausted.

   %Add in storage N demand to pool

    if i==303
      Ndemand_plus_storage=AMstorageNdemand(i,:)+AMtotalNdemand(i,:);
      %Add in storage C into NPP or C available pool
      npp_plus_storage=AMNPP(i,:)+AMstorageC(1,:);
    else
       npp_plus_storage=AMNPP(i,:);
       Ndemand_plus_storage=AMtotalNdemand(i,:);
    end


    AMplantCN(i2,:)=npp_plus_storage./Ndemand_plus_storage;
    AMplantCN(i2,isnan(AMplantCN(i2,:)))=0;

    AMCavailable(i2,:)=npp_plus_storage;
    AMNdeficit(i2,:)=Ndemand_plus_storage-AMfree(i2,:);
    xx=(AMNdeficit(i2,:)<=0);

    AMNacq(i2,xx)=0;
    AMCacq(i2,xx)=0;
    AMCgrowth(i2,xx)=AMCgrowth(i2,xx);

    xx=(AMNdeficit(i2,:)>0);
    AMCacq(i2,xx)=(AMCavailable(i2,xx)-AMfree(i2,xx).*AMplantCN(i,xx))./(1+AMplantCN(i,xx)./AMcost_acq(i2,xx));
    AMNacq(i2,xx)=min(AMNdeficit(i2,xx),AMCacq(i2,xx)./AMcost_acq(i2,xx));
    AMCgrowth(i2,xx)=AMCavailable(i2,xx)-(AMNacq(i2,xx).*AMcost_acq(i2,xx));

    AMtotal_Nacq_active(i2,xx)=AMCacq(i2,xx)./AMcost_active(i2,xx);
    AMtotal_Nacq_non_myc(i2,xx)=AMCacq(i2,xx)./AMcost_non_myc(i2,xx);
    AMtotal_Nacq_resorb(i2,xx)=AMCacq(i2,xx)./AMcost_resorb(i2,xx);
    AMtotal_Nacq_fix(i2,xx)=AMCacq(i2,xx)./AMcost_fixN(i2,xx);
    AMtotal_Nacq_storage(i2,xx)=AMCacq(i2,xx)./AMcost_storage(xx);


    if any(AMsoilN(i2,:)<0)
        error('soilN<0')
    end

    if exist('prescribed_AM_rhizoCflux','var')
        AMrhizoCflux(i2,:)=prescribed_AM_rhizoCflux(i2,:);
    else
        AMrhizoCflux(i2,:)=(AMtotal_Nacq_non_myc(i2,:)+AMtotal_Nacq_fix(i2,:)).*AMcost_acq(i2,:);
    end

    if exist('prescribed_AM_fungal_prod','var')
        AMfungalprod(i2,:)=prescribed_AM_fungal_prod(i2,:);
    else
        AMfungalprod(i2,:)=(AMtotal_Nacq_active(i2,:)).*AMcost_acq(i2,:);
    end

    % Add rhizo C flux to rhizosphere as fast C
    if do_exudation

        AM_rhizospheres.litterC(:,1)=AM_rhizospheres.litterC(:,1)+AMrhizoCflux(i2,:)';

             rhiz_fungi_frac=AM_rhizo_percent(i,:)'/100;
            AM_rhizospheres.litterC(:,1)=AM_rhizospheres.litterC(:,1)+AMfungalprod(i2,:)'.*rhiz_fungi_frac.*(1.0-litter_fungi_frac);
            AM_bulksoils.litterC(:,1)=AM_bulksoils.litterC(:,1)+AMfungalprod(i2,:)'.*(1.0-rhiz_fungi_frac).*(1.0-litter_fungi_frac);
            AM_litterlayers.litterC(:,1)=AM_litterlayers.litterC(:,1)+AMfungalprod(i2,:)'.*litter_fungi_frac;

            if exudate_CN>0
              AM_rhizospheres.litterN(:,1)=AM_rhizospheres.litterN(:,1)+AMfungalprod(i2,:)'.*rhiz_fungi_frac.*(1.0-litter_fungi_frac)/exudate_CN;
              AM_bulksoils.litterN(:,1)=AM_bulksoils.litterN(:,1)+AMfungalprod(i2,:)'.*(1.0-rhiz_fungi_frac).*(1.0-litter_fungi_frac)/exudate_CN;
              AM_litterlayers.litterN(:,1)=AM_litterlayers.litterN(:,1)+AMfungalprod(i2,:)'.*litter_fungi_frac/exudate_CN;
            end
    end


    AMtotal_Nacq_storage(i2,:)=min(AMtotal_Nacq_storage(i2,:),AMstorageN(i2,:));

    xx=AMtotal_Nacq_resorb(i2,:)>AMleafN2(1,:); %Truncate resorb if exceeds leaf N

    AMtotal_Nacq_resorb(i2,xx)=AMleafN2(1,xx);
    AMCgrowth(i2,xx)=AMCavailable(i2,xx)-(AMNacq(i2,xx).*AMcost_acq(i2,xx));


     %Update Storage N pool to include retranslocation
     if i<303
        AMstorageN(i2+1,:)=AMstorageN(i2,:)-AMtotal_Nacq_storage(i2,:);
     else
        AMstorageN(i2+1,:)=AMstorageN(i2,:)+AMtotal_Nacq_resorb(i2,:)+AMtotal_Nacq_active(i2,:)+AMtotal_Nacq_non_myc(i2,:)+AMtotal_Nacq_fix(i2,:);  %%%%Fix this for soil N uptake%%%%

        xx=AMleafN2(:)>0;

        AMresorbper(1,xx)=1-(AMtotal_Nacq_resorb(303,xx)./AMleafN2(1,xx));
        AMlitter_productionN(i2,xx)=(AMleafN(i,xx).*AMresorbper(1,xx));
        AMlitter_productionCN(i2,xx)=AMlitter_production(i,xx)./AMlitter_productionN(i2,xx);

        xx=AMleafN2(:)<=0;
        AMresorbper(1,xx)=0.0;
        AMlitter_productionN(i2,xx)=0.0;
        AMlitter_productionCN(i2,xx)=0.0;

    end

    %Add in ability to get N from storage in leaf_out period
    if i<303
      AMstorageNmob(i2,:)=min(AMstorageN(i2,:),AMtotal_Nacq_storage(i2,:));
      AMstorageN(i2+1,:)=AMstorageN(i2,:)-AMstorageNmob(i2,:);%Update storage N amount
    end

    if do_litter
      %Add litter to soil pools.
      if constant_litter_inputs
        litter_N_input=meanlitter_production(i)./repmat(AMlitter_productionCN(i2,:),[3,1]).*repmat(AM_litter_composition,[number,1])';
        litter_N_input(~isfinite(litter_N_input))=0.0;

      else
          litter_N_input=repmat(AMlitter_productionN(i2,:),[3,1]).*repmat(AM_litter_composition,[number,1])';
      end
      AM_litterlayers.litterN=AM_litterlayers.litterN+litter_N_input';

      if constant_litter_inputs
          litter_C_input=meanlitter_production(i).*repmat(AM_litter_composition,[number,1])';
      else
          litter_C_input=repmat(AMlitter_production(i,:),[3,1]).*repmat(AM_litter_composition,[number,1])';
      end

      AM_litterlayers.litterC=AM_litterlayers.litterC+litter_C_input';

      %Root litter
      if constant_litter_inputs
          rootlitter_N_input=repmat(mean(daily_froot_turnover(i,:))/1000.*(1-per_ECM(:,2)'/100)./AMrootCN,[3,1]).*repmat(AM_rootlitter_composition,[number,1])';
      else
          rootlitter_N_input=repmat(daily_froot_turnover(i,:)/1000.*(1-per_ECM(:,2)'/100)./AMrootCN,[3,1]).*repmat(AM_rootlitter_composition,[number,1])';
      end

      if constant_litter_inputs
          rootlitter_C_input=repmat(mean(daily_froot_turnover(i,:))/1000.*(1-per_ECM(:,2)'/100),[3,1]).*repmat(AM_rootlitter_composition,[number,1])';
      else
          rootlitter_C_input=repmat(daily_froot_turnover(i,:)/1000.*(1-per_ECM(:,2)'/100),[3,1]).*repmat(AM_rootlitter_composition,[number,1])';
      end


      if any(rootlitter_C_input>0)
%                 disp(i)
      end

      if root_inputs_only_rhiz
          AM_rhizospheres.litterC=AM_rhizospheres.litterC+rootlitter_C_input';
          AM_rhizospheres.litterN=AM_rhizospheres.litterN+rootlitter_N_input';
      else
          rf=repmat(AM_rhizo_percent(i,:)/100,[3,1]);
          AM_rhizospheres.litterC=AM_rhizospheres.litterC+rootlitter_C_input'.*rf';
          AM_bulksoils.litterC=AM_bulksoils.litterC+rootlitter_C_input'.*(1.0-rf)';
          AM_rhizospheres.litterN=AM_rhizospheres.litterN+rootlitter_N_input'.*rf';
          AM_bulksoils.litterN=AM_bulksoils.litterN+rootlitter_N_input'.*(1.0-rf)';
      end

    end

    check_cohorts(AM_rhizospheres);

    AMtotal_Nacq(i2,:)=AMNacq(i2,:)+AMfree(i2,:);
    AMtotal_Cgrowth(i2,:)=sum(AMCgrowth(i2,:));

    total_soil_N_uptake(i2,:)=ECMtotal_Nacq_active(i2,:)+ECMtotal_Nacq_non_myc(i2,:)+AMtotal_Nacq_active(i2,:)+AMtotal_Nacq_non_myc(i2,:);

    fac=min(mineralN'./total_soil_N_uptake(i2,:),1);



    total_soil_N_uptake(i2,fac<1)=mineralN(fac<1);

    %Truncate active uptake if exceeds soil N
    ECMtotal_Nacq_active(i2,:)= ECMtotal_Nacq_active(i2,:).*fac;
    ECMtotal_Nacq_non_myc(i2,:)= ECMtotal_Nacq_non_myc(i2,:).*fac;
    AMtotal_Nacq_active(i2,:)= AMtotal_Nacq_active(i2,:).*fac;
    AMtotal_Nacq_non_myc(i2,:)= AMtotal_Nacq_non_myc(i2,:).*fac;

    % N acquisition comes out of soil inorganic N pool
    AMsoilN(i2,:)=AMsoilN(i2,:) - AMtotal_Nacq_active(i2,:) - AMtotal_Nacq_non_myc(i2,:);
    if any(AMsoilN(i2,:)<0)
        error('soilN<0')
    end

    ECMsoilN(i2,:)=ECMsoilN(i2,:) - ECMtotal_Nacq_active(i2,:) - ECMtotal_Nacq_non_myc(i2,:);
    if any(ECMsoilN(i2,:)<0)
        error('soilN<0')
    end

    if(any(mineralN<0))
        error(mineralN<0)
    end

    total_mineralN(i2,:)=mineralN(:);

    mineralN(:)=(mineralN(:) - total_soil_N_uptake(i2,:)')*(1-mineral_N_lost_fraction);
    mineralN(:) = mineralN(:) + atmo_N_dep; % Add daily atmospheric N deposition


end % End of main loop


t_elapsed=toc;
fprintf('Time elapsed: %f seconds\n',t_elapsed)

%% Final stuff

% Save data for restarting simulation (e.g. spun up values)
if exist('write_restart_file','var')
    clear data_saved
    fprintf('Writing restart file to %s\n',write_restart_file);
    data_saved.AM_rhizospheres=AM_rhizospheres;
    data_saved.AM_mineral_N=AMsoilN(i2,:);
    data_saved.AM_bulksoil=AM_bulksoils;
    data_saved.AMstorageN=AMstorageN(i2,:);
    data_saved.AM_litterlayers=AM_litterlayers;
    data_saved.ECM_rhizospheres=ECM_rhizospheres;
    data_saved.ECM_mineral_N=ECMsoilN(i2,:);
    data_saved.ECM_bulksoil=ECM_bulksoils;
    data_saved.ECMstorageN=ECMstorageN(i2,:);
    data_saved.ECM_litterlayers=ECM_litterlayers;
    save(write_restart_file,'-struct','data_saved')
end

ECM_soil_outputs=add_outputs(ECM_rhiz_outputs,ECM_bulk_outputs);
ECM_total_outputs=add_outputs(ECM_soil_outputs,ECM_litter_outputs);

AM_soil_outputs=add_outputs(AM_rhiz_outputs,AM_bulk_outputs);
AM_total_outputs=add_outputs(AM_soil_outputs,AM_litter_outputs);

total_soil_outputs=add_outputs(AM_soil_outputs,ECM_soil_outputs);
total_outputs=add_outputs(AM_total_outputs,ECM_total_outputs);
total_litter_outputs=add_outputs(AM_litter_outputs,ECM_litter_outputs);

if do_exudation
    AM_rhiz_outputs_e=AM_rhiz_outputs;ECM_rhiz_outputs_e=ECM_rhiz_outputs;
    AM_bulk_outputs_e=AM_bulk_outputs;ECM_bulk_outputs_e=ECM_bulk_outputs;
    AM_total_outputs_e=AM_total_outputs;ECM_total_outputs_e=ECM_total_outputs;
    AM_litter_outputs_e=AM_litter_outputs;ECM_litter_outputs_e=ECM_litter_outputs;
    AM_soil_outputs_e=AM_soil_outputs;ECM_soil_outputs_e=ECM_soil_outputs;
    total_outputs_e=total_outputs; total_soil_outputs_e=total_soil_outputs;
else
    AM_rhiz_outputs_noe=AM_rhiz_outputs;ECM_rhiz_outputs_noe=ECM_rhiz_outputs;
    AM_bulk_outputs_noe=AM_bulk_outputs;ECM_bulk_outputs_noe=ECM_bulk_outputs;
    AM_total_outputs_noe=AM_total_outputs;ECM_total_outputs_noe=ECM_total_outputs;
    AM_litter_outputs_noe=AM_litter_outputs;ECM_litter_outputs_noe=ECM_litter_outputs;
    AM_soil_outputs_noe=AM_soil_outputs;ECM_soil_outputs_noe=ECM_soil_outputs;
    total_outputs_noe=total_outputs; total_soil_outputs_noe=total_soil_outputs;
end



%%%%%%%Variable Output Key%%%%%%%%%%%%
% ALL C or N units ARE IN KG/M2!!!!!
%AM or ECM_rhizopercent=percent of soil composed by AM or ECM roots
%AM or ECMCacq=total C spent on getting N from air, soil, leaves
%AM or ECMrhizoCflux=total C spent on getting N belowground
%AM or ECMCavailable=total C avaliable for growth or N acquisition
%AM or ECMfree= total N taken up from passive uptake or storage
%AM or ECMlitter_production= C in litter inputs
%AM or ECMlitter_productionCN=C to N ratio of litter
%AM or ECMlitter+productionN= N in litter inputs
%AM or ECMNacq= total N taken up by expending C
%AM or ECMNdeficit= N need but not taken up
%AM or ECMpassive= N taken up passively through transpiration stream
%AM or ECMrootCN= C to N of fine roots
%AM or ECMsoilN= soilN availble for uptake %%%%%This needs to be updated to
%reflect N taken out at the end of each loop
%AM or ECMstorageC= amount of C in storage
%AM or ECMstorageN= amount of N in storage
%AM or ECM_total_Nacq= total N taken up free and at cost
%AM or ECM_total_Nacq_fix= total N taken up from fixation
%AM or ECM_total_Nacq_active= total N taken up by mycorrhizal roots
%AM or ECM_total_Nacq_non_myc= total N taken up by non mycorrhizal roots
%AM or ECM_total_Nacq_non_myc= total N taken up from leaves upon senescence
%AM or ECMtotalNdemand= total amount of N to support plant nutrition
%daily_fine_root_turnover= amount of roots that die on a given timestep
%%%%%%This needs to be updated to be split into AM and ECM root litter
