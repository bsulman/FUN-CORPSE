function [totalsoilC,totalsoilN] = total_pool(outputs,include_protected)

if ~exist('include_protected','var')
    include_protected=true;
end

totalsoilC=outputs.unprotectedC.fast+outputs.unprotectedC.slow+outputs.unprotectedC.deadmic+outputs.livingMicrobeC;
    
totalsoilN=outputs.unprotectedN.fast+outputs.unprotectedN.slow+outputs.unprotectedN.deadmic+outputs.livingMicrobeC./8.0;
    
    
    if include_protected
        totalsoilC=totalsoilC+outputs.protectedC.fast+outputs.protectedC.slow+outputs.protectedC.deadmic;
        totalsoilN=totalsoilN+outputs.protectedN.fast+outputs.protectedN.slow+outputs.protectedN.deadmic;
    end
end