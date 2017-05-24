function check_cohorts(cohorts)
    bad=false;
    if any(cohorts.litterC(:)<0) || any(cohorts.litterN(:)<0) || any(cohorts.protectedC(:)<0) || any(cohorts.protectedN(:)<0)
        bad=true
    end
    
    if any(cohorts.livingMicrobeC(:)<0) || any(cohorts.livingMicrobeN(:)<0) || any(cohorts.CO2(:)<0) || any(cohorts.Rtot(:)<0) || any(cohorts.originalLitterC(:)<0) || any(cohorts.originalLitterN(:)<0)
        bad=true
    end
    
    if any(~isfinite(cohorts.litterC(:))) || any(~isfinite(cohorts.litterN(:))) || any(~isfinite(cohorts.protectedC(:))) || any(~isfinite(cohorts.protectedN(:)))
        bad=true
    end
    
    if any(~isfinite(cohorts.livingMicrobeC(:))) || any(~isfinite(cohorts.livingMicrobeN(:))) || any(~isfinite(cohorts.CO2(:))) || any(~isfinite(cohorts.Rtot(:))) || any(~isfinite(cohorts.originalLitterC(:))) || any(~isfinite(cohorts.originalLitterN(:)))
        bad=true
    end
    
    if bad
        disp(cohorts)
        error('Cohort bad')
    end
        
end