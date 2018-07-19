%% 3A. See u.PlotRTHistogram.m

%% 3B_1. Boxplots of RT, only considering trials w/ stim or trials immediately following stim

ContAll = [];
IpsiAll = [];

for i = 1:length(SL)
    if (~isempty(SL(i).Bad) || strcmp(SL(i).Condition, 'Control')...
            || strcmp(SL(i).Condition, 'NaN') || strcmp(SL(i).Condition, 'nostim')...
            || isempty(SL(i).trig1) || strcmp(SL(i).Condition(end),'R'))
        continue;
    end
    
    if(length(SL(i).Condition)<6 || ~strcmp(SL(i).Condition(1:6),'Contra'))
        continue;
    end
    
    if(str2num(SL(i).Condition(end))>0)
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
    
    % get all trials that had stimulation
    if(strcmp(SL(i).Animal,'Ubi'))
        StimT = u.getStimTrials(SL(i).trig1,ContT);
    else
        CondStart = find(ContT(:,2)<SL(i).trig1(1),1,'last');
        CondEnd = find(ContT(:,1)>SL(i).trig1(end),1);
        StimT = CondStart:(CondEnd-1);
    end
    
    % get pre and post stim RT 
    PreContRT = find(ContT(:,2)<SL(i).trig1(1),1,'last');
    PreContRT = ContRT(1:PreContRT);
    PreIpsiRT = find(IpsiT(:,2)<SL(i).trig1(1),1,'last');
    PreIpsiRT = ContRT(1:PreIpsiRT);
    PostContRT = find(ContT(:,1)>SL(i).trig1(end),1);
    PostContRT = ContRT(PostContRT:end);
    PostIpsiRT = find(IpsiT(:,1)>SL(i).trig1(end),1);
    PostIpsiRT = ContRT(PostIpsiRT:end);
    
    % put trials in order
    Trials = [ContT(:,1);IpsiT(:,1)];
    RT = [ContRT;IpsiRT];
    Stim = zeros(1,length(Trials));
    Stim(StimT) = 1; 
    Labels = [ones(1,length(ContT)),zeros(1,length(IpsiT))]; % 1 is contra, 0 is ipsi
    
    [Trials,order] = sort(Trials);
    RT = RT(order);
    Stim = Stim(order);
    Labels = Labels(order);
    
    % get all contra and ipsi within 1 trial of stim
    start = find(Stim,1); finish = find(Stim,1,'last')+1;
    CondContRT = []; CondIpsiRT = [];
    for t = start:finish
        if(~(Stim(t-1)))
            continue;
        end
        
        if(Labels(t))
            CondContRT(end+1) = RT(t);
        else
            CondIpsiRT(end+1) = RT(t);
        end
        
    end
    
    ContAll(end+1,1) = nanmedian(PreContRT);
    ContAll(end,2) = nanmedian(CondContRT);
    ContAll(end,3) = nanmedian(PostContRT);
    
    IpsiAll(end+1,1) = nanmedian(PreIpsiRT);
    IpsiAll(end,2) = nanmedian(CondIpsiRT);
    IpsiAll(end,3) = nanmedian(PostIpsiRT);
    
end

ContCond = ContAll(:,2)-ContAll(:,1);
IpsiCond = IpsiAll(:,2)-IpsiAll(:,1);
ContPost = ContAll(:,3)-ContAll(:,1);
IpsiPost = IpsiAll(:,3)-IpsiAll(:,1);

figure;%('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); 

subplot(2,2,1)
boxplot([ContAll(:,1),ContAll(:,2),ContAll(:,3)],'Labels',{'ContPre','ContCond','ContPost'})
title('Raw')
subplot(2,2,2)
boxplot([IpsiAll(:,1),IpsiAll(:,2),IpsiAll(:,3)],'Labels',{'IpsiPre','IpsiCond','IpsiPost'})
title('Raw')

subplot(2,1,2)
boxplot([ContCond,IpsiCond,ContPost,IpsiPost],'Labels',{'Contra Cond','Ipsi Cond','Contra Post','Ipsi Post'})
yl = ylim; title('Subtracting Pre')
hold on; plot([2.5,2.5],yl,'k--')



