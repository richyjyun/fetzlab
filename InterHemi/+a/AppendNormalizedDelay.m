function SL = AppendNormalizedDelay(SL,Save)
% For the single SL passed in, calculates the difference between the
% triggers and the calculated reaction time for those trials, and
% appends to the SL as a subfield. Hemi dictates which hemisphere is
% stimulated.
if(nargin < 2)
    Save = 0;
end

for i = 1:length(SL)
    if(isempty(SL(i).trig1)) % needs to have had a conditioning session
        continue;
    end
    
    D = SL(i).Date;
    Session = char(D);
    disp(['Normalized Delay Session ',Session])
    
    Filter = [];
    Trials = [];
    delays = [];
    
    hand = 0; % 0 for left hand, 1 for right hand
    if strcmp(SL(i).StimHemi,'L')
        if any(char(SL(i).Condition)=='C')
            hand = 1;
        elseif ~any(char(SL(i).Condition)=='I')
            SL(i).NormDelay = NaN;
            continue;
        end
    elseif strcmp(SL(i).StimHemi,'R')
        if any(char(SL(i).Condition)=='I')
            hand = 1;
        elseif ~any(char(SL(i).Condition)=='C')
            SL(i).NormDelay = NaN;
            continue;
        end
    else
        SL(i).NormDelay = NaN;
        continue;
    end
    
    if hand
        Filter = u.FilterAccForMovementTimes(SL(i).accel_raw_r, SL(i).fs, 'richardson');
        Trials = SL(i).righttrials;
        RT = SL(i).rts_r;
    else
        Filter = u.FilterAccForMovementTimes(SL(i).accel_raw_l, SL(i).fs, 'richardson');
        Trials = SL(i).lefttrials;
        RT = SL(i).rts_l;
    end
    
    %     if any(char(SL(i).Condition)=='C')
    %         Filter = u.FilterAccForMovementTimes(SL(i).accel_raw_r, SL(i).fs, 'richardson');
    %         Trials = SL(i).righttrials;
    %     elseif any(char(SL(i).Condition)=='I')
    %         Filter = u.FilterAccForMovementTimes(SL(i).accel_raw_l, SL(i).fs, 'richardson');
    %         Trials = SL(i).lefttrials;
    %     else
    %         continue;
    %     end
    
    delays = nan(1,length(SL(i).trig1));
    
    window = 0.6;
    window = floor(window*SL(i).fs);
    for j = 1:length(SL(i).trig1)
        
        norm = Trials(:,1)-((SL(i).trig1(j))+50);
        ind = find(norm>0,1)-1;
        
        if isempty(ind) || isnan(RT(ind))
            continue;
        end
        % Added in case trig1 happens to be slightly before the start of
        % trial
        if(abs(Trials(ind,2)-SL(i).trig1(j)) > abs(Trials(ind+1,1)-SL(i).trig1(j)))
            ind = ind+1;
        end

        if(strcmp(SL(i).Animal,'Igor'))
            delays(j) = SL(i).trig1(j) - (Trials(ind,1) + RT(ind));%find(snip == max(snip),1)*1000/SL(i).fs);
        else
            delays(j) = SL(i).trig1(j) + str2num(SL(i).Stim_Delay)/1000*SL(i).fs - (Trials(ind,1) + RT(ind));%find(snip == max(snip),1)*1000/SL(i).fs);
        end
        delays(j) = SL(i).trig1(j) - (Trials(ind,1));
    end
    SL(i).NormDelay = nanmedian(delays);
    SL(i).NormSE = nanstd(delays)/sqrt(sum(~isnan(delays)));
%     disp(num2str(sum(isnan(delays))/length(delays)));
    if(Save)
        SL.Delays = delays;
    end
end
end