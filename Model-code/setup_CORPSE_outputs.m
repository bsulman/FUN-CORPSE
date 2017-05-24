function CORPSE_outputs = setup_CORPSE_outputs(ndays,nplots)

CORPSE_outputs.unprotectedC.fast=zeros(ndays,nplots);
CORPSE_outputs.unprotectedC.slow=zeros(ndays,nplots);
CORPSE_outputs.unprotectedC.deadmic=zeros(ndays,nplots);
CORPSE_outputs.protectedC.fast=zeros(ndays,nplots);
CORPSE_outputs.protectedC.slow=zeros(ndays,nplots);
CORPSE_outputs.protectedC.deadmic=zeros(ndays,nplots);
CORPSE_outputs.unprotectedN.fast=zeros(ndays,nplots);
CORPSE_outputs.unprotectedN.slow=zeros(ndays,nplots);
CORPSE_outputs.unprotectedN.deadmic=zeros(ndays,nplots);
CORPSE_outputs.protectedN.fast=zeros(ndays,nplots);
CORPSE_outputs.protectedN.slow=zeros(ndays,nplots);
CORPSE_outputs.protectedN.deadmic=zeros(ndays,nplots);
CORPSE_outputs.N_mineralization=zeros(ndays,nplots);
CORPSE_outputs.N_immobilization=zeros(ndays,nplots);
CORPSE_outputs.livingMicrobeC=zeros(ndays,nplots);
CORPSE_outputs.decomp.fast=zeros(ndays,nplots);
CORPSE_outputs.decomp.slow=zeros(ndays,nplots);
CORPSE_outputs.decomp.deadmic=zeros(ndays,nplots);
CORPSE_outputs.CO2prod=zeros(ndays,nplots);

end
