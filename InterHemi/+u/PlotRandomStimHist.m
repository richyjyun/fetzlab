% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbiFinal.mat');

clearvars -except SL UbiSL IgorSL KatoSL

ContAll = {};
IpsiAll = {};
StimAllC = {};
CondAll = {};
StimAllI = {};
AllTrig = [];
IgorBadDates = [20120214;20120301;20120419;20120604];

for i = 1:length(SL)
    if (~isempty(SL(i).Bad) || strcmp(SL(i).Condition, 'Control')...
             || strcmp(SL(i).Condition, 'nostim')...
            || isempty(SL(i).trig1))
        continue;
    end
    
    if(~(strcmp(SL(i).Condition,'Contra_R') || strcmp(SL(i).Condition,'Ipsi_R')))
        continue;
    end


    if(strcmp(SL(i).StimHemi,'L'))
        ContRT = SL(i).rts_r;
        IpsiRT = SL(i).rts_l;
        ContT = SL(i).righttrials;
        IpsiT = SL(i).lefttrials;
    else
        ContRT = SL(i).rts_l;
        IpsiRT = SL(i).rts_r;
        ContT = SL(i).lefttrials;
        IpsiT = SL(i).righttrials;
    end
    
    if(strcmp(SL(i).Condition(1:4),'Ipsi'))
        Trial = IpsiT;
    else
        Trial = ContT;
    end
    
    % get all trials that had stimulation
    if(strcmp(SL(i).Animal,'Ubi')) %Kato trigger times are odd, so just assume all stim
        [StimT,idx,stimT] = u.getStimTrials(SL(i).trig1,Trial,50);
    else
%         CondStart = find(ContT(:,1)<SL(i).trig1(1),1,'last');
%         CondEnd = find(ContT(:,1)>SL(i).trig1(end),1);
%         StimT = CondStart:(CondEnd-1);
        [StimT,idx,stimT] = u.getStimTrials(SL(i).trig1,Trial,500);
    end
    trig = SL(i).trig1(idx);
    AllTrig(1,end+1) = length(SL(i).trig1);
    AllTrig(2,end) = length(StimT);
    
    
    % Set stim times
    if(strcmp(SL(i).Animal,'Ubi'))
        stimT = stimT + str2num(SL(i).Stim_Delay);
    end
    
    if(strcmp(SL(i).Condition(1:4),'Ipsi'))
        StimAllI{end+1} = stimT;
    else
        StimAllC{end+1} = stimT;
    end
    
end

%% Histogram
StimC = [StimAllC{:}];
StimI = [StimAllI{:}];

bins = 0:100:900;
figure; histogram(StimC,bins); hold on; histogram(StimI,bins);
title('Random Stim Delay Times')
ylabel('Count'); xlabel('Delay (ms)');
legend('Contra','Ipsi')



