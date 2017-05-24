%Add/subtract cohorts together
function out=add_cohorts(cohort1,cohort2)
out.litterC=cohort1.litterC + cohort2.litterC; %gC of [fast,slow,dead microbe]
out.litterN=cohort1.litterN + cohort2.litterN; %gN of [fast,slow,dead microbe]
out.protectedC=cohort1.protectedC + cohort2.protectedC;
out.protectedN=cohort1.protectedN + cohort2.protectedN;
out.livingMicrobeC=cohort1.livingMicrobeC + cohort2.livingMicrobeC;
out.livingMicrobeN=cohort1.livingMicrobeN + cohort2.livingMicrobeN;
out.CO2=cohort1.CO2 + cohort2.CO2;
out.Rtot=cohort1.Rtot + cohort2.Rtot;
%These are for conservation checks
out.originalLitterC=cohort1.originalLitterC + cohort2.originalLitterC;
out.originalLitterN=cohort1.originalLitterN + cohort2.originalLitterN;

end