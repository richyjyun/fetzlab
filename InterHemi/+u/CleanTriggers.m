function SL = CleanTriggers(SL)

% Cleans stimulation triggers by removing all that are not within an
% appropriate trial and by ensuring there is only one trigger per trial

for i = 1:length(SL)
    
    if(isempty(SL(i).trig1) || ~isempty(SL(i).Bad) || strcmp(SL(i).Condition,'nostim') || ...
            strcmp(SL(i).Condition,'tonic') || strcmp(SL(i).Condition,'Control') || strcmp(SL(i).Condition,'NaN'))
        continue;
    end
    
    if strcmp(SL(i).StimHemi,'R')
        if strcmp(SL(i).Condition(1:4),'Ipsi')
            Trials = SL(i).righttrials;
        else
            Trials = SL(i).lefttrials;
        end
    else
        if strcmp(SL(i).Condition(1:4),'Ipsi')
            Trials = SL(i).lefttrials;
        else
            Trials = SL(i).righttrials;
        end
    end

    good = ones(1,length(SL(i).trig1));
    stim = zeros(1,length(Trials));
    for j = 1:length(SL(i).trig1)
        ind = find(SL(i).trig1(j)>Trials(:,1)-50 & SL(i).trig1(j)<Trials(:,2)+50,1);%,'last');
        if(isempty(ind))
            good(j) = 0;
            continue;
        end
        if(stim(ind) == 0)
            stim(ind) = 1;
        else
            good(j) = 0;
        end
    end
    
    SL(i).trig1 = SL(i).trig1(logical(good));
    
end

end