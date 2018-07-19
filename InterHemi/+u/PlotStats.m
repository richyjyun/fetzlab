%% get control 

Cont = {};
Ipsi = {};
NTrials = [];

for i = 1:length(SL)
    cond = char(SL(i).Condition);
    
    if(~isempty(SL(i).Bad) || ~(strcmp(cond,'nostim') || strcmp(cond,'Control')...
            || (strcmp(cond,'NaN') && strcmp(SL(i).Animal,'Kato'))))
        continue;
    end
    
    if(str2num(SL(i).Date) == 20170405) %random delay stim
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

    % put trials in order
    Trials = [ContT(:,1);IpsiT(:,1)];
    RT = [ContRT;IpsiRT];
    Labels = [ones(1,length(ContT)),zeros(1,length(IpsiT))]; % 1 is contra, 0 is ipsi
    
    [Trials,order] = sort(Trials);
    RT = RT(order);
    Labels = Labels(order);
    
    %remove first 50 trials
    NTrials(end+1) = length(Trials);
    rm = 50;
    Trials = Trials(rm:end);
    RT = RT(rm:end);
    Labels = Labels(rm:end);
    
    bounds = [round(length(Trials)/3),round(2*length(Trials)/3)];
    C = find(Labels); I = find(Labels==0);
    
    Cont{1,end+1} = RT(C(C<bounds(1)));
    Cont{2,end} = RT(C(C>=bounds(1) & C< bounds(2)));
    Cont{3,end} = RT(C(C>=bounds(2)));
    
    Ipsi{1,end+1} = RT(I(I<bounds(1)));
    Ipsi{2,end} = RT(I(I>=bounds(1) & I< bounds(2)));
    Ipsi{3,end} =  RT(I(I>=bounds(2)));
end

% Consolidate the data
CondC = []; CondI = [];
for i = 1:size(Cont,2)
    CondC = [CondC;Cont{2,i} - nanmedian(Cont{1,i})];
    CondI = [CondI;Ipsi{2,i} - nanmedian(Ipsi{1,i})];
end

RT = [CondC;CondI]; 
idx = [ones(length(CondC),1);2*ones(length(CondI),1)];

% Plot
positions = sort([1:length(unique(idx))/2,(1:length(unique(idx))/2)+1/3]);
figure;
subplot(2,1,1);
boxplot(RT,idx, 'Notch','on', 'positions', positions,'Whisker',0,'Symbol','');
hold on; xl = xlim;
line(xl,[0,0],'linestyle','--','color',[0,0,0])

% Color boxes
color = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980]; % just using default matlab colors
color = repmat(color,length(unique(idx))/2,1); color = flipud(color); % need to flip because call is backwards
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.5); % set face colors
end

% legend for colors
c = get(gca, 'Children');
[~,leg] = legend(c(1:2), 'Contra', 'Ipsi' );
PatchInLegend = findobj(leg, 'type', 'patch');
set(PatchInLegend, 'facealpha', 0.5);

% Replotting to have box plot on top
boxplot(RT,idx, 'Notch','on', 'positions', positions,'Whisker',0,'Symbol','')

% labels
xticks([positions]);
xticklabels([sum(idx==1),sum(idx==2)]);
title([SL(i).Animal,'Control \DeltaRT']); ylabel('RT (ms)'); 

ylim([-150,150])

ControlRT = RT; Controlidx = idx; 

%% get data
ContAll = {};
IpsiAll = {};
StimAllC = {};
CondAll = {};
StimAllI = {};
AllTrig = [];

for i = 1:length(SL)
    if (~isempty(SL(i).Bad) || strcmp(SL(i).Condition, 'Control')...
            || strcmp(SL(i).Condition, 'NaN') || strcmp(SL(i).Condition, 'nostim')...
            || isempty(SL(i).trig1) || strcmp(SL(i).Condition(end),'R'))
        continue;
    end
    

    if(length(SL(i).Condition)<6 || ~strcmp(SL(i).Condition(1:6),'Contra'))
        continue;
    end
    
    if(strcmp(SL(i).Animal,'Ubi') && (str2num(SL(i).Date) < 20170226 || str2num(SL(i).Date) > 20170306))
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
    if(strcmp(SL(i).Animal,'Ubi')) %Kato trigger times are odd, so just assume all stim
        [StimT,idx,stimT] = u.getStimTrials(SL(i).trig1,ContT,50);
    else
%         CondStart = find(ContT(:,1)<SL(i).trig1(1),1,'last');
%         CondEnd = find(ContT(:,1)>SL(i).trig1(end),1);
%         StimT = CondStart:(CondEnd-1);
        [StimT,idx,stimT] = u.getStimTrials(SL(i).trig1,ContT,500);
    end
    trig = SL(i).trig1(idx);
    AllTrig(1,end+1) = length(SL(i).trig1);
    AllTrig(2,end) = length(StimT);
    
    % get pre and post stim RT
    PreContRT = find(ContT(:,2)<trig(1),1,'last');
    PreContRT = ContRT(1:PreContRT);
    PreIpsiRT = find(IpsiT(:,2)<trig(1),1,'last');
    PreIpsiRT = IpsiRT(1:PreIpsiRT);
    PostContRT = find(ContT(:,1)>trig(end),1);
    PostContRT = ContRT(PostContRT:end);
    PostIpsiRT = find(IpsiT(:,1)>trig(end),1);
    PostIpsiRT = IpsiRT(PostIpsiRT:end);
    
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
    
    
    % Set stim times
    if(strcmp(SL(i).Animal,'Ubi'))
        stimT = stimT + str2num(SL(i).Stim_Delay);
    end
    
    % get all contra and ipsi within 1 trial of stim
    % get ipsi trial stim time from GO
    start = find(Stim,1); finish = find(Stim,1,'last')+1;
    CondContRT = []; CondIpsiRT = []; StimI = []; StimC=[];
    since = 1;
    for t = start:finish
        if(~(Stim(t-since) && ~any(Stim(t-(since-1):t-1)))) %only looking at trials immediately after stim
            continue;
        end
        
%         if(~Stim(t))
%             continue;
%         end
        
        if(Labels(t))
            CondContRT(end+1) = RT(t);
            stimind = sum(Stim(1:t-since));
            StimC(end+1) = stimT(stimind);
        else
            CondIpsiRT(end+1) = RT(t);
            stimind = sum(Stim(1:t-since));
            StimI(end+1) = stimT(stimind);
        end
    end
    
    % Save RT
    ContAll{end+1,1} = PreContRT;
    ContAll{end,2} = CondContRT;
    ContAll{end,3} = PostContRT';
    
    IpsiAll{end+1,1} = PreIpsiRT;
    IpsiAll{end,2} = CondIpsiRT;
    IpsiAll{end,3} = PostIpsiRT';
    
    StimAllC{end+1} = StimC;
    StimAllI{end+1} = StimI;
    
    CondAll{end+1,1} = SL(i).Condition;
    CondAll{end,2} = SL(i).Stim_Delay;
    CondAll{end,3} = SL(i).StimHemi;
    
end

% consolidate data
Cstim = []; Istim = []; CRT = []; IRT = [];
CpreRT = cellfun(@nanmean,ContAll(:,1)); 
IpreRT = cellfun(@nanmean,IpsiAll(:,1)); 
bins = [-100,200; 200,500; 500,800];

nDays = zeros(1,size(bins,1)*2);
for i = 1:length(StimAllC)
    if(strcmp(SL(1).Animal,'Ubi') && i==14)
        continue;
    end
    
    stimT = StimAllC{i};
    if(strcmp(SL(1).Animal,'Ubi'))
        jitterT = 0;
    else
        jitterT = 0;
    end
    
    condRT = ContAll{i,2};
    Cstim = [Cstim,stimT];
    for b = 1:size(bins,1)
        if(any(stimT >= bins(b,1) & stimT < bins(b,2)))
            nDays(b*2-1) = nDays(b*2-1)+1;
        end
    end
    CRT = [CRT,condRT-CpreRT(i)];
    
    stimT = StimAllI{i};
    condRT = IpsiAll{i,2};
    Istim = [Istim,stimT];
    for b = 1:size(bins,1)
        if(any(stimT >= bins(b,1) & stimT < bins(b,2)))
            nDays(b*2) = nDays(b*2)+1;
        end
    end
    IRT = [IRT,condRT-IpreRT(i)];
end
% -100 to 100, 200-400, 400-600 and 700-900
% bins = -100:200:900;
labels = {};
for i = 1:size(bins,1)
    labels{i} = sprintf('%d-%d',bins(i,1),bins(i,2));
end

Cidx = zeros(1,length(Cstim));
Iidx = zeros(1,length(Istim));

for i = 1:size(bins,1)
    Cidx(Cstim >= bins(i,1) & Cstim < bins(i,2)) = i;
    Iidx(Istim >= bins(i,1) & Istim < bins(i,2)) = i;
end

% [~,~,Cidx] = histcounts(Cstim,bins);
% [~,~,Iidx] = histcounts(Istim,bins);

CRT(Cidx==0) = []; Cidx(Cidx==0) = []; Cidx = Cidx*2-1;
IRT(Iidx==0) = []; Iidx(Iidx==0) = []; Iidx = Iidx*2;

RT = [CRT,IRT]; idx = [Cidx,Iidx];


% Plot
positions = sort([1:length(unique(idx))/2,(1:length(unique(idx))/2)+1/3]);
subplot(2,1,2);
boxplot(RT,idx, 'Notch','on', 'positions', positions,'Whisker',0,'Symbol','');
hold on; xl = xlim;
line(xl,[0,0],'linestyle','--','color',[0,0,0])

% Color boxes
color = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980]; % just using default matlab colors
color = repmat(color,length(labels),1); color = flipud(color); % need to flip because call is backwards
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.5); % set face colors
end

% legend for colors
c = get(gca, 'Children');
[~,leg] = legend(c(1:2), 'Contra', 'Ipsi' );
PatchInLegend = findobj(leg, 'type', 'patch');
set(PatchInLegend, 'facealpha', 0.5);

% Replotting to have box plot on top
boxplot(RT,idx, 'Notch','on', 'positions', positions,'Whisker',0,'Symbol','')

% Set labels
labelpos = reshape(positions,2,length(positions)/2);
labelpos = labelpos(1,:)+diff(labelpos)/2;
set(gca,'xtick',labelpos)
set(gca,'xticklabel',labels)

% labels
title([SL(i).Animal,'Contra Exp \DeltaRT']); ylabel('RT (ms)'); xlabel('Stim Time (ms)')

ylim([-150,150])


%% Get stats
str = ''; alpha = 0.05;
for i = unique(idx)
    Cidx = 1; 
    if(mod(i,2) == 0)
        Cidx = 2;
    end
    [~,p] = ttest2(ControlRT(Controlidx==Cidx),RT(idx==i));
    sig = 'ns';
    if( p < alpha )
        sig = 's';
    end
    n = sum(idx==i);
    str = sprintf('%sBox %d: %0.3e     %s    %d\n',str,i,p,sig,n);
end
hold on; xl = xlim; yl = ylim;
text(xl(1),yl(1),str,'verticalalignment','bottom');







