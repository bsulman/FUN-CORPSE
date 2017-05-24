%Multiply a cohort times a value
function out=multiply_cohort(cohort,value)
out.litterC=cohort.litterC*value; %gC of [fast,slow,dead microbe]
out.litterN=cohort.litterN*value; %gN of [fast,slow,dead microbe]
out.protectedC=cohort.protectedC*value;
out.protectedN=cohort.protectedN*value;
out.livingMicrobeC=cohort.livingMicrobeC*value;
out.livingMicrobeN=cohort.livingMicrobeN*value;
out.CO2=cohort.CO2*value;
out.Rtot=cohort.Rtot*value;
%These are for conservation checks
out.originalLitterC=cohort.originalLitterC*value;
out.originalLitterN=cohort.originalLitterN*value;

end

