clear
% rehash path

nyears=50;
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

% Separate microbial C:N for AM and ECM based on Meghan's measurements
params_AM=params;
params_AM.CN_microb=6.1;

params_ECM=params;
params_ECM.CN_microb=8.5;


% Baseline simulation
NPP_mult = 1.0
rhiz_mult = 1.0

data_file_out='FUN-CORPSE-outputs-spinup.mat'

% This script was run repeatedly, saving new restart files, until simulations had reached steady state
use_restart_file='restart_noexud_moreMB_2.mat';
write_restart_file='restart_noexud_moreMB_3.mat'

do_exudation=false;
Daily_FUN_MODEX_Code_combined_N;

use_restart_file='restart_exud_moreMB_2.mat';
write_restart_file='restart_exud_moreMB_3.mat'
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
