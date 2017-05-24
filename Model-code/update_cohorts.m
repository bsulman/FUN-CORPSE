%Cohorts should be 1D array of structs with these fields:
    %Carbon pools
    %litterC  				%Carbon substrate. Relative amounts of C species are important for resp rate
    %litterN  				%Nitrogen substrate
    %protectedC 				%Aggregate and mineral complex C
    %protectedN 				%Aggregate and mineral complex C
    %livingMicrobeC
    %livingMicrobeN             					 %Current carbon mass of live microbes
    %CO2                       					 %Cumulative CO2 generated by decomposition
    %Rtot                        				 %Cumulative decomposition (includes double counting through microbial products)
    %Original C, for conservation checks
    %originalLitterC
    %originalLitterN         					  !Keep track for carbon balance check



function [cohorts,outputs] = update_cohorts(cohorts,T,theta,mineralN,dt,params)

% if length(cohorts) ~= length(mineralN)
%     error('Length of cohorts must equal length of mineralN')
% end

ncohorts=45;



mineralN=mineralN(:);

if T<0
    error('T<0')
end
if theta<0 || theta>1
    error('Theta must be between 0 and 1')
end
if any(mineralN<0)
    error('mineralN<0')
end

% global params;
totalResp=zeros(ncohorts,3);
totalN_decomposed=zeros(ncohorts,3);


et=params.et;
eup=repmat(params.eup,[ncohorts,1]);

cohorts_array=cohorts;
check_cohorts(cohorts_array);

totalN_start=sum(mineralN)+sum(sum(cohorts_array.litterN))+sum(sum(cohorts_array.protectedN))+sum(cohorts_array.livingMicrobeN);
totalC_start=sum(sum(cohorts_array.litterC))+sum(sum(cohorts_array.protectedC))+sum(cohorts_array.livingMicrobeC)+sum(cohorts_array.CO2);

%Calculate maximum potential C decomposition rate
potential_tempResp=Resp(cohorts_array.litterC,cohorts_array.livingMicrobeC,T,theta);

%Carbon loss cannot exceed size of pool
potential_tempResp(dt.*potential_tempResp > cohorts_array.litterC) = cohorts_array.litterC(dt.*potential_tempResp > cohorts_array.litterC)/dt;

%N associated with the C decomposition rate
pot_tempN_decomposed = potential_tempResp.*cohorts_array.litterN./cohorts_array.litterC; %kgN/m2/yr
pot_tempN_decomposed(cohorts_array.litterC<=0.0)=0.0;


%Maximum N immobilization rate
IMM_N_max=mineralN*params.max_immobilization_rate*365;   %kg/m2/yr
xx=IMM_N_max>mineralN/dt;
IMM_N_max(xx)=mineralN(xx)/dt;

microbeTurnover=max(0.0,(cohorts_array.livingMicrobeC-params.minMicrobeC*sum(cohorts_array.litterC,2))/params.Tmic);   %kg/m2/yr


% Carbon and nitrogen available for microbial growth
carbon_supply=sum(potential_tempResp.*repmat(params.eup,[ncohorts,1]),2);
nitrogen_supply=sum(pot_tempN_decomposed.*repmat(params.nup,[ncohorts,1]),2);

maintenance_resp=microbeTurnover*(1.0-params.et);
overflow_resp=zeros(ncohorts,1);

%Microbial growth
loc_Nlim=(carbon_supply - maintenance_resp > (nitrogen_supply+IMM_N_max)*params.CN_microb);
CN_imbalance_term=zeros(ncohorts,1);

% Growth is nitrogen limited, with not enough mineral N to support it with max immobilization
CN_imbalance_term(loc_Nlim) = -IMM_N_max(loc_Nlim);
cohorts_array.livingMicrobeC(loc_Nlim) = cohorts_array.livingMicrobeC(loc_Nlim) + dt.*((nitrogen_supply(loc_Nlim)+IMM_N_max(loc_Nlim))*params.CN_microb - microbeTurnover(loc_Nlim)*et);
cohorts_array.livingMicrobeN(loc_Nlim) = cohorts_array.livingMicrobeN(loc_Nlim) + dt.*(nitrogen_supply(loc_Nlim)+IMM_N_max(loc_Nlim) - microbeTurnover(loc_Nlim)*et/params.CN_microb);
overflow_resp(loc_Nlim)=carbon_supply(loc_Nlim)-maintenance_resp(loc_Nlim) - (nitrogen_supply(loc_Nlim)+IMM_N_max(loc_Nlim))*params.CN_microb;

loc_immob=(carbon_supply - maintenance_resp >= nitrogen_supply*params.CN_microb)&(carbon_supply - maintenance_resp < (nitrogen_supply+IMM_N_max)*params.CN_microb);
% Growth must be supported by immobilization of some mineral nitrogen, but is ultimately carbon limited
CN_imbalance_term(loc_immob) = -((carbon_supply(loc_immob)-maintenance_resp(loc_immob))/params.CN_microb - nitrogen_supply(loc_immob));
cohorts_array.livingMicrobeC(loc_immob) = cohorts_array.livingMicrobeC(loc_immob) + dt*(carbon_supply(loc_immob) - microbeTurnover(loc_immob));
cohorts_array.livingMicrobeN(loc_immob) = cohorts_array.livingMicrobeN(loc_immob) + dt*((carbon_supply(loc_immob)-maintenance_resp(loc_immob))/params.CN_microb - microbeTurnover(loc_immob)*et/params.CN_microb);

loc_Clim=~(loc_Nlim | loc_immob);
% Growth is carbon limited -- extra N is mineralized
cohorts_array.livingMicrobeC(loc_Clim) = cohorts_array.livingMicrobeC(loc_Clim) + dt*(carbon_supply(loc_Clim) - microbeTurnover(loc_Clim)) ;
cohorts_array.livingMicrobeN(loc_Clim) = cohorts_array.livingMicrobeN(loc_Clim) + dt*((carbon_supply(loc_Clim)-maintenance_resp(loc_Clim))/params.CN_microb - microbeTurnover(loc_Clim)*params.et/params.CN_microb);
CN_imbalance_term(loc_Clim) = nitrogen_supply(loc_Clim) - (carbon_supply(loc_Clim)-maintenance_resp(loc_Clim))/params.CN_microb;



tempResp=potential_tempResp;
tempN_decomposed=pot_tempN_decomposed;


%Good time to check if microbe C:N was preserved, and C and N are >= 0

%Microbial turnover
deadmic_C_produced=dt.*microbeTurnover*et;   %actual fraction of microbial turnover
cohorts_array.litterC(:,3)=cohorts_array.litterC(:,3)+deadmic_C_produced*(1-params.frac_turnover_slow);  %kg/m2
cohorts_array.litterC(:,2)=cohorts_array.litterC(:,2)+deadmic_C_produced*(params.frac_turnover_slow);
turnover_N_produced=dt.*microbeTurnover*et/params.CN_microb;
turnover_N_min=turnover_N_produced*params.frac_turnover_min;
turnover_N_slow=turnover_N_produced*params.frac_turnover_slow;
deadmic_N_produced=turnover_N_produced-turnover_N_min-turnover_N_slow;
cohorts_array.litterN(:,3)=cohorts_array.litterN(:,3)+deadmic_N_produced;  %kg/m2
cohorts_array.litterN(:,2)=cohorts_array.litterN(:,2)+turnover_N_slow;


%CO2 production and cumulative CO2 produced by cohort
CO2prod=dt.*(sum(tempResp.*(1.0-eup),2)+maintenance_resp+overflow_resp); %kg/m2
cohorts_array.CO2=cohorts_array.CO2+CO2prod;  %kg/m2

IMM_N_gross=zeros(ncohorts,1);
MINER_gross=zeros(ncohorts,1);

xx=find(CN_imbalance_term>0.0);
nup=repmat(params.nup,[ncohorts,1]);
IMM_N_gross(xx)=0.0;
MINER_gross(xx)=sum((1-nup(xx,:)).*tempN_decomposed(xx,:),2)+CN_imbalance_term(xx)+turnover_N_min(xx)/dt;

xx=find(CN_imbalance_term<=0.0);
IMM_N_gross(xx)=-CN_imbalance_term(xx);
MINER_gross(xx)=sum((1-nup(xx,:)).*tempN_decomposed(xx,:),2)+turnover_N_min(xx)/dt;



% Update the amount of organic C and N in the cohort after the decomposition process

cohorts_array.litterC=cohorts_array.litterC-dt.*tempResp;     %kg/m2
totalResp=totalResp+tempResp;             %kg/m2/yr

cohorts_array.litterN=cohorts_array.litterN-dt.*tempN_decomposed;     %kg/m2
totalN_decomposed=totalN_decomposed+tempN_decomposed;   %kg/m2/yr

IMM_Nprod=IMM_N_gross*dt;  		%kg/m2	!Gross immobilization of inorganic nitrogen for the cohort
MINERAL_prod=MINER_gross*dt;	%kg/m2	!Gross mineralization of organic nitrogen for the cohort

mineralN = mineralN + MINER_gross.*dt - IMM_N_gross.*dt;



%Update protected carbon
protectedCturnover = cohorts_array.protectedC./params.tProtected ;
protectedNturnover = cohorts_array.protectedN./params.tProtected  ;

newProtectedC=zeros(ncohorts,3);
prate=repmat(params.protection_rate,[ncohorts,1]);
newProtectedC(cohorts_array.litterC>0.0)=prate(cohorts_array.litterC>0.0).*cohorts_array.litterC(cohorts_array.litterC>0.0)*dt;

newProtectedC(newProtectedC>cohorts_array.litterC) =cohorts_array.litterC(newProtectedC>cohorts_array.litterC);
cohorts_array.protectedC = cohorts_array.protectedC + newProtectedC - dt.*protectedCturnover;

cohorts_array.litterC = cohorts_array.litterC - newProtectedC + dt.*protectedCturnover;

protected_produced=newProtectedC;  %kg/m2
protected_turnover_rate=protectedCturnover;  %kg/m2/dt


cohorts_array.Rtot=cohorts_array.Rtot+sum(totalResp,2);


newProtectedN=zeros(ncohorts,3);
newProtectedN(cohorts_array.litterN>0.0)=prate(cohorts_array.litterN>0.0).*cohorts_array.litterN(cohorts_array.litterN>0.0)*dt;


newProtectedN(newProtectedN>cohorts_array.litterN) =cohorts_array.litterN(newProtectedN>cohorts_array.litterN);
cohorts_array.protectedN = cohorts_array.protectedN + newProtectedN - dt.*protectedNturnover;
cohorts_array.litterN = cohorts_array.litterN - newProtectedN + dt.*protectedNturnover;

protected_N_produced=newProtectedN ; %!kg/m2
protected_N_turnover_rate=protectedNturnover ;% !kg/m2/dt

% cohorts=unzip_cohorts(cohorts_array);
cohorts=cohorts_array;

outputs.decomp=totalResp;
outputs.totalN_decomposed=totalN_decomposed;
outputs.protected_produced=protected_produced;
outputs.protected_N_produced=protected_N_produced;
outputs.N_immobilization=IMM_Nprod;
outputs.N_mineralization=MINERAL_prod;
outputs.mineralN=mineralN;
outputs.protected_turnover_rate=protected_turnover_rate;
outputs.protected_N_turnover_rate=protected_N_turnover_rate;
outputs.CO2prod=CO2prod;

% fprintf('Total C mismatch = %1.2e; Total N = %1.2e\n',totalC(cohort)-cohort.originalLitterC,...
%     totalN(cohort)+mineralN-cohort.originalLitterN)

totalN_end=sum(mineralN)+sum(sum(cohorts_array.litterN+cohorts_array.protectedN))+sum(cohorts_array.livingMicrobeN);
if abs(totalN_end-totalN_start)>1e-10
  fprintf('N diff = %g (%g%%)\n',totalN_end-totalN_start,(totalN_end-totalN_start)/totalN_start*100)
  error('N not conserved')
end

totalC_end=sum(sum(cohorts_array.litterC+cohorts_array.protectedC))+sum(cohorts_array.livingMicrobeC)+sum(cohorts_array.CO2);
if abs(totalC_end-totalC_start)>1e-10
  fprintf('C diff (end) = %g (%g%%)\n',totalC_end-totalC_start,(totalC_end-totalC_start)/totalC_start*100)
  error('C not conserved')
end

%Decomposition rate
function Resp = Resp(Ctotal,Chet,T,theta)
    %Chet        %heterotrophic (microbial) C, living microbes biomasa in the cohort
    %T,theta     %temperature (k), theta (fraction of 1.0)
    %Ctotal      %Substrate C (3-value vector)

%    global params




    if(theta==0.0)
        Resp=0.0;
        return
    end

    Resp=zeros(ncohorts,3);
    notzero=(sum(Ctotal,2)~=0 | Chet~=0);

    Chet2=repmat(Chet,[1,3]);
    kC=repmat(params.kC,[ncohorts,1]);
    Vmax=repmat(Vmax(T),[ncohorts,1]);

    Resp(notzero,:)=Vmax.*theta.^3.*(Ctotal(notzero,:)).*Chet2(notzero,:)./(Ctotal(notzero,:).*kC(notzero,:)+Chet2(notzero,:)).*(1-theta).^params.gas_diffusion_exp;

end


function Vmax=Vmax(T)
%     global params
    Tref=293.15;
    Rugas=8.314472;

     %Normalization value
    alpha=params.vmaxref./exp(-params.Ea/(Rugas.*Tref));
    Vmax=alpha.*exp(-params.Ea./(Rugas.*T));

end


function out=totalC(cohort)
    out=sum(cohort.litterC+cohort.protectedC)+cohort.livingMicrobeC+cohort.CO2;
end

function out=totalN(cohort)
    out=sum(cohort.litterN+cohort.protectedN)+cohort.livingMicrobeN;
end





function out=unzip_cohorts(cohorts)
    % put arrays back into array of cohorts


    for jj=1:ncohorts
        outstruct=struct;
        outstruct.litterC=cohorts.litterC(jj,:);
        outstruct.litterN=cohorts.litterN(jj,:);
        outstruct.protectedC=cohorts.protectedC(jj,:);
        outstruct.protectedN=cohorts.protectedN(jj,:);
        outstruct.livingMicrobeC=cohorts.livingMicrobeC(jj);
        outstruct.livingMicrobeN=cohorts.livingMicrobeN(jj);
        outstruct.CO2=cohorts.CO2(jj);
        outstruct.Rtot=cohorts.Rtot(jj);
        outstruct.originalLitterC=cohorts.originalLitterC(jj);
        outstruct.originalLitterN=cohorts.originalLitterN(jj);
        out(jj)=outstruct;
    end
end

end
