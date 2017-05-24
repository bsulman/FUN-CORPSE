%Multiply a cohort times a value
function out=multiply_cohorts(cohort,value)
out.litterC=cohort.litterC.*repmat(value,[1,3]); %gC of [fast,slow,dead microbe]
out.litterN=cohort.litterN.*repmat(value,[1,3]); %gN of [fast,slow,dead microbe]
out.protectedC=cohort.protectedC.*repmat(value,[1,3]);
out.protectedN=cohort.protectedN.*repmat(value,[1,3]);
out.livingMicrobeC=cohort.livingMicrobeC.*value;
out.livingMicrobeN=cohort.livingMicrobeN.*value;
out.CO2=cohort.CO2.*value;
out.Rtot=cohort.Rtot.*value;
%These are for conservation checks
out.originalLitterC=cohort.originalLitterC.*value;
out.originalLitterN=cohort.originalLitterN.*value;

end

