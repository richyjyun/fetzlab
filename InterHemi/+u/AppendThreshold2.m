function SL = AppendThreshold2(SL)

for i = 1:length(SL)
    
    D = SL(i).Date;
    
    Session = char(D);
    disp(['Threshold Session ',Session])
    
    RFilter = u.FilterAccForMovementTimes(SL(i).accel_raw_r, SL(i).fs*SL.dwn, 'richardson');
    LFilter = u.FilterAccForMovementTimes(SL(i).accel_raw_l, SL(i).fs*SL.dwn, 'richardson');
    
    RSnips = [];
    LSnips = [];
    window = 0.6;
    window = floor(window.*SL(i).fs*SL.dwn);
    for j = 1:length(SL(i).righttrials)
        if(SL(i).righttrialsuccess(j))
            if((SL(i).righttrials(j,2)-window)>0 && SL(i).righttrials(j,2)< length(RFilter))
                ind = floor(SL(i).righttrials(j,2));%*SL(i).fs/1000);
                RSnips(end+1,:) = RFilter(ind-window:ind);%+window);
                if(max(RSnips(end,:))>1000)
                    RSnips(end,:) = [];
                end
            end
        end
    end
    for j = 1:length(SL(i).lefttrials)
        if(SL(i).lefttrialsuccess(j))
            if((SL(i).lefttrials(j,2)-window)>0 && SL(i).lefttrials(j,2)< length(LFilter))
                ind = floor(SL(i).lefttrials(j,2));%*SL(i).fs/1000);
                if(ind-window <= 0)
                    continue;
                end
                LSnips(end+1,:) = LFilter(ind-window:ind);
                if(max(LSnips(end,:))>1000)
                    LSnips(end,:) = [];
                end
            end
        end
    end
    %     threshold = 3/4*max((mean([RSnips;LSnips]))); % make threshold 2/3 of maximum value?
    %     SL(i).Threshold = threshold;
    m = [max(median(LSnips)),max(median(RSnips))];
    SL(i).Max = m;
end

end