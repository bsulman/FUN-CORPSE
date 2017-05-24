function outputs = record_CORPSE_outputs(outputs,cohort,otherdata,timepoint,plotnum)


outputs.unprotectedC.fast(timepoint,plotnum)=cohort.litterC(1);
outputs.unprotectedC.slow(timepoint,plotnum)=cohort.litterC(2);
outputs.unprotectedC.deadmic(timepoint,plotnum)=cohort.litterC(3);
outputs.protectedC.fast(timepoint,plotnum)=cohort.protectedC(1);
outputs.protectedC.slow(timepoint,plotnum)=cohort.protectedC(2);
outputs.protectedC.deadmic(timepoint,plotnum)=cohort.protectedC(3);
outputs.unprotectedN.fast(timepoint,plotnum)=cohort.litterN(1);
outputs.unprotectedN.slow(timepoint,plotnum)=cohort.litterN(2);
outputs.unprotectedN.deadmic(timepoint,plotnum)=cohort.litterN(3);
outputs.protectedN.fast(timepoint,plotnum)=cohort.protectedN(1);
outputs.protectedN.slow(timepoint,plotnum)=cohort.protectedN(2);
outputs.protectedN.deadmic(timepoint,plotnum)=cohort.protectedN(3);
outputs.N_mineralization(timepoint,plotnum)=otherdata.N_mineralization;
outputs.N_immobilization(timepoint,plotnum)=otherdata.N_immobilization;
outputs.livingMicrobeC(timepoint,plotnum)=cohort.livingMicrobeC;

outputs.decomp.fast(timepoint,plotnum)=otherdata.decomp(1);
outputs.decomp.slow(timepoint,plotnum)=otherdata.decomp(2);
outputs.decomp.deadmic(timepoint,plotnum)=otherdata.decomp(3);

outputs.N_decomposed.fast(timepoint,plotnum)=otherdata.totalN_decomposed(1);
outputs.N_decomposed.slow(timepoint,plotnum)=otherdata.totalN_decomposed(2);
outputs.N_decomposed.deadmic(timepoint,plotnum)=otherdata.totalN_decomposed(3);


end