function SL = AppendNormalizedDelayCursor(SL)
% For the single SL passed in, calculates the difference between the
% triggers and the peak of the movements (using accelerometer signal) and
% appends to the SL as a subfield. Hemi dictates which hemisphere is
% stimulated.
for i = 1:length(SL)
    if(isempty(SL(i).trig1)) % needs to have had a conditioning session
        continue;
    end
    
    Trials = [];
    if(any(char(SL(i).Condition)=='I'))
        Trials = SL(i).lefttrials;
        Success = SL(i).lefttrialsuccess;
    else
        Trials = SL(i).righttrials;
        Success = SL(i).righttrialsuccess;
    end
    
    delays = nan(1,length(SL(i).trig1));
    
    window = 1;
    window = floor(window*SL(i).fs);
    for j = 1:length(SL(i).trig1)
        temp = abs(Trials(:,1) - SL(i).trig1(j));
        trial_ind = find(temp == min(temp));
        if(Success(trial_ind))
            start = Trials(trial_ind,1);
            finish = min(start+window,length(SL(i).cursorintarget));
            ind = find(SL(i).cursorintarget(start:finish),1);
            ind = (SL(i).trig1(j) + str2num(SL(i).Stim_Delay)) - (ind+start);
            if(~isempty(ind))
                delays(j) = ind;
            end
        end
    end
    SL(i).NormDelay = nanmedian(delays);
end
end