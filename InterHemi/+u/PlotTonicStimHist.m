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
    
    d = str2num(SL(i).Date);
    
    if( ~((d>=20120423 && d<=20120507) || (d>=20111008 && d<=20111024)))
        continue;
    end
    
    if(~(strcmp(SL(i).Condition,'NaN') && isempty(SL(i).Notes) && length(SL(i).trig1)>1))
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
    
    [StimT,idx,stimT] = u.getStimTrials(SL(i).trig1,ContT,50);
    StimAllC{end+1} = stimT;
    
    
    [StimT,idx,stimT] = u.getStimTrials(SL(i).trig1,IpsiT,50);
    StimAllI{end+1} = stimT;
    
    
end

%% Histogram
StimC = [StimAllC{:}];
StimI = [StimAllI{:}];

bins = 0:100:900;
figure; histogram(StimC,bins); hold on; histogram(StimI,bins);
title('Tonic Stim Delay Times')
ylabel('Count'); xlabel('Delay (ms)');
legend('Contra','Ipsi')



