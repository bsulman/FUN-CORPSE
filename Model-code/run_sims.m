clear
% rehash path

nyears=10;
do_litter=true;
litter_transfer_to_soil=0.1/365; %Transfer rate per day
initial_litter_fraction=0.16;
rhiz_fungi_frac=0.7;
litter_fungi_frac=0.2;
plots=1:45;plots_to_run=plots;
root_inputs_only_rhiz=false;
mineral_N_lost_fraction=0.1;
cost_param_mult = 1/2.5;

exudate_CN=30;

atmo_N_dep=6/100.^2/365; % kg/m2/day

constant_litter_inputs=false;
constant_NPP=true;

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

% Separate microbial C:N for AM and ECM based on Midgley et al (2016)
params_AM=params;
params_AM.CN_microb=6.1;

params_ECM=params;
params_ECM.CN_microb=8.5;


% Baseline simulation
NPP_mult = 1.0
rhiz_mult = 1.0
nyears=50;

data_file_out='../../Model-output/FUN-CORPSE-outputs-NPP-1.0.mat'

use_restart_file='../Restart-files/restart_exud_moreMB_3.mat'

do_exudation=false;
Daily_FUN_MODEX_Code_combined_N;

do_exudation=true;
Daily_FUN_MODEX_Code_combined_N;

totalrhizoCflux=ECMrhizoCflux+AMrhizoCflux;
totalfungalprod=ECMfungalprod+AMfungalprod;
AMNacq_noresorb=AMNacq-AMtotal_Nacq_resorb;
ECMNacq_noresorb=ECMNacq-ECMtotal_Nacq_resorb;
save(data_file_out,'total_outputs_e','total_outputs_noe','totalrhizoCflux','ECMrhizoCflux','AMrhizoCflux',...
    'AM_rhiz_outputs_e','AM_rhiz_outputs_noe','ECM_rhiz_outputs_e','ECM_rhiz_outputs_noe',...
    'AM_bulk_outputs_e','AM_bulk_outputs_noe','ECM_bulk_outputs_e','ECM_bulk_outputs_noe','per_ECM',...
    'ECMNacq_noresorb','AMNacq_noresorb','total_mineralN','totalfungalprod','ECMfungalprod','AMfungalprod','-v7')

nyears=10;

% Reversed ECM/AM plants simulation
NPP_mult = 1.0
rhiz_mult = 1.0
reverse_plant_myco = true
params_ECM.CN_microb=6.1;
params_AM.CN_microb=8.5;
nyears=50;

data_file_out='../../Model-output/FUN-CORPSE-outputs-NPP-1.0-myc-reversed.mat'


do_exudation=false;
Daily_FUN_MODEX_Code_combined_N;

do_exudation=true;
Daily_FUN_MODEX_Code_combined_N;

totalrhizoCflux=ECMrhizoCflux+AMrhizoCflux;
totalfungalprod=ECMfungalprod+AMfungalprod;
AMNacq_noresorb=AMNacq-AMtotal_Nacq_resorb;
ECMNacq_noresorb=ECMNacq-ECMtotal_Nacq_resorb;
save(data_file_out,'total_outputs_e','total_outputs_noe','totalrhizoCflux','ECMrhizoCflux','AMrhizoCflux',...
    'AM_rhiz_outputs_e','AM_rhiz_outputs_noe','ECM_rhiz_outputs_e','ECM_rhiz_outputs_noe',...
    'AM_bulk_outputs_e','AM_bulk_outputs_noe','ECM_bulk_outputs_e','ECM_bulk_outputs_noe','per_ECM',...
    'ECMNacq_noresorb','AMNacq_noresorb','total_mineralN','totalfungalprod','ECMfungalprod','AMfungalprod','-v7')

% Reset things for next simulations
reverse_myco_costs=false;
reverse_plant_myco = false;
params_AM.CN_microb=6.1;
params_ECM.CN_microb=8.5;
nyears=10;

%% NPP and litter elevated 20%
NPP_mult = 1.2;
rhiz_mult = 1.0;
litter_mult=1.2;

data_file_out='../../Model-output/FUN-CORPSE-outputs-NPP-1.2.mat'

do_exudation=false;
Daily_FUN_MODEX_Code_combined_N;
do_exudation=true;
Daily_FUN_MODEX_Code_combined_N;

totalrhizoCflux=ECMrhizoCflux+AMrhizoCflux;
totalfungalprod=ECMfungalprod+AMfungalprod;
AMNacq_noresorb=AMNacq-AMtotal_Nacq_resorb;
ECMNacq_noresorb=ECMNacq-ECMtotal_Nacq_resorb;
save(data_file_out,'total_outputs_e','total_outputs_noe','totalrhizoCflux','ECMrhizoCflux','AMrhizoCflux',...
    'AM_rhiz_outputs_e','AM_rhiz_outputs_noe','ECM_rhiz_outputs_e','ECM_rhiz_outputs_noe',...
    'AM_bulk_outputs_e','AM_bulk_outputs_noe','ECM_bulk_outputs_e','ECM_bulk_outputs_noe','per_ECM',...
    'ECMNacq_noresorb','AMNacq_noresorb','total_mineralN','totalfungalprod','ECMfungalprod','AMfungalprod','-v7')


%% NPP and litter elevated 20%, exud from baseline sim
NPP_mult = 1.2
rhiz_mult = 1.0
litter_mult=1.2;

x = load('../../Model-output/FUN-CORPSE-outputs-NPP-1.0.mat','AMrhizoCflux');
prescribed_AM_rhizoCflux=x.AMrhizoCflux;
x = load('../../Model-output/FUN-CORPSE-outputs-NPP-1.0.mat','ECMrhizoCflux');
prescribed_ECM_rhizoCflux=x.ECMrhizoCflux ;
x = load('../../Model-output/FUN-CORPSE-outputs-NPP-1.0.mat','AMfungalprod');
prescribed_AM_fungal_prod=x.AMfungalprod;
x = load('../../Model-output/FUN-CORPSE-outputs-NPP-1.0.mat','ECMfungalprod');
prescribed_ECM_fungal_prod=x.ECMfungalprod ;

data_file_out='../../Model-output/FUN-CORPSE-outputs-NPP-1.2-same-exud.mat'

do_exudation=false;
Daily_FUN_MODEX_Code_combined_N;
do_exudation=true;
Daily_FUN_MODEX_Code_combined_N;


totalrhizoCflux=ECMrhizoCflux+AMrhizoCflux;
totalfungalprod=ECMfungalprod+AMfungalprod;
AMNacq_noresorb=AMNacq-AMtotal_Nacq_resorb;
ECMNacq_noresorb=ECMNacq-ECMtotal_Nacq_resorb;
save(data_file_out,'total_outputs_e','total_outputs_noe','totalrhizoCflux','ECMrhizoCflux','AMrhizoCflux',...
    'AM_rhiz_outputs_e','AM_rhiz_outputs_noe','ECM_rhiz_outputs_e','ECM_rhiz_outputs_noe',...
    'AM_bulk_outputs_e','AM_bulk_outputs_noe','ECM_bulk_outputs_e','ECM_bulk_outputs_noe','per_ECM',...
    'ECMNacq_noresorb','AMNacq_noresorb','total_mineralN','totalfungalprod','ECMfungalprod','AMfungalprod','-v7')
