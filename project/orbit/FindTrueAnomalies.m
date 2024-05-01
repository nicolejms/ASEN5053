function [TA1f,TA2f] = FindTrueAnomalies(deltaTA, TA1, TA2)
    
    deltaTAall = [deltaTA; 360-deltaTA];
    TA1all = [TA1; -TA1; 360-TA1];
    TA2all = [TA2; -TA2; 360-TA2];
    tolerance = 1e-4;
    for i = 1:length(TA2all)
        for j = 1:length(TA1all)
            dta = TA2all(i) - TA1all(j);
            if abs(dta-deltaTA) < tolerance || abs(dta-(360-deltaTA)) < tolerance
                TA1f = TA1all(j);
                TA2f = TA2all(i);
                return;
            end
        end
    end
    
end
