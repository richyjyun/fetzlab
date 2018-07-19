function SL = AppendReactionTimesCursor(SL)
%
% For each SL index, adds two fields called rt_l and rt_r that are as long
% as left_trials and right_trials, respectively
%
% Finds reaction time as first instance cursor enters the target box for
% each trial

for i = 1:length(SL)
    
    if(isempty(SL(i).cursorintarget))
        continue;
    end
    
    window = 1000;
    rt_l = nan(1,length(SL(i).lefttrials));
    rt_r = nan(1,length(SL(i).righttrials));
    for j=1:max(length(SL(i).lefttrials),length(SL(i).righttrials))
        if(j <= length(SL(i).lefttrials) && SL(i).lefttrialsuccess(j))
            start  = SL(i).lefttrials(j,1);
            finish = min(SL(i).lefttrials(j,1)+window,length(SL(i).cursorintarget));
            time = find(SL(i).cursorintarget(start:finish),1);
            if(~isempty(time))
                rt_l(j) = time;
            end
        end
        if(j <= length(SL(i).righttrials) && SL(i).righttrialsuccess(j))
            start  = SL(i).righttrials(j,1);
            finish = min(SL(i).righttrials(j,1)+window,length(SL(i).cursorintarget));
            time = find(SL(i).cursorintarget(start:finish),1);
            if(~isempty(time))
                rt_r(j) = time;
            end
        end
    end
    
    SL(i).rts_l = rt_l;
    SL(i).rts_r = rt_r;
    
end
