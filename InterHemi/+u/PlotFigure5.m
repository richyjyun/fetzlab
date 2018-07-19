%% Figure 5. RT effects relative to time since last stim (relative to RP)
% A. CS_cont

ContAll = [];
IpsiAll = [];
ContTime = {}; ContTimeRT = {};
IpsiTime = {}; IpsiTimeRT = {};
CondAll = {};

for i = 1:length(SL)
    if (~isempty(SL(i).Bad) || strcmp(SL(i).Condition, 'Control')...
            || strcmp(SL(i).Condition, 'NaN') || strcmp(SL(i).Condition, 'nostim')...
            || isempty(SL(i).trig1) || strcmp(SL(i).Condition(end),'R'))
        continue;
    end
        
%     if(length(SL(i).Condition)<6 || ~strcmp(SL(i).Condition(1:6),'Contra'))
%         continue;
%     end
    
    if(length(SL(i).Condition)<4 || ~strcmp(SL(i).Condition(1:4),'Ipsi'))
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
    
    % get pre and post stim RT
    PreContRT = find(ContT(:,2)<SL(i).trig1(1),1,'last');
    PreContRT = ContRT(1:PreContRT);
    PreIpsiRT = find(IpsiT(:,2)<SL(i).trig1(1),1,'last');
    PreIpsiRT = IpsiRT(1:PreIpsiRT);
    PostContRT = find(ContT(:,1)>SL(i).trig1(end),1);
    PostContRT = ContRT(PostContRT:end);
    PostIpsiRT = find(IpsiT(:,1)>SL(i).trig1(end),1);
    PostIpsiRT = IpsiRT(PostIpsiRT:end);
    
    % get stim triggers
    if(strcmp(SL(i).Animal,'Ubi'))
        trig = SL(i).trig1+SL(i).Stim_Delay;
    else
        trig = SL(i).trig1;
    end
    
    % find time of last stim for contra trials
    CondContRT = []; CondContTime = [];
    for t = 1:length(ContT)
        RP = ContT(t,1)+ContRT(t); % can fix later to make it RP time rather than RT
        ind = find(trig < RP , 1,'last');
        if(isnan(RP) || isempty(ind))
            continue;
        end
        CondContTime(end+1) = trig(ind)-ContT(t,1);
        CondContRT(end+1) = ContRT(t);
    end
    
    % find time of last stim for ipsi trials
    CondIpsiRT = []; CondIpsiTime = [];
    for t = 1:length(IpsiT)
        RP = IpsiT(t,1)+IpsiRT(t); % can fix later to make it RP time rather than RT
        ind = find(trig < RP , 1,'last');
        if(isnan(RP) || isempty(ind))
            continue;
        end
        CondIpsiTime(end+1) = trig(ind)-IpsiT(t,1);
        CondIpsiRT(end+1) = IpsiRT(t);
    end
    
    ContTime{end+1} = CondContTime;
    ContTimeRT{end+1} = CondContRT;
    IpsiTime{end+1} = CondIpsiTime;
    IpsiTimeRT{end+1} = CondIpsiRT;
    
    ContAll(end+1,1) = nanmedian(PreContRT);
    ContAll(end,2) = nanmedian(CondContRT);
    ContAll(end,3) = nanmedian(PostContRT);
    
    IpsiAll(end+1,1) = nanmedian(PreIpsiRT);
    IpsiAll(end,2) = nanmedian(CondIpsiRT);
    IpsiAll(end,3) = nanmedian(PostIpsiRT);
    
    CondAll{end+1,1} = SL(i).Condition;
    CondAll{end,2} = SL(i).Stim_Delay;
    CondAll{end,3} = SL(i).StimHemi;
    CondAll{end,4} = SL(i).Date;
end

%% Scatter plot
figure; ax1 = subplot(2,1,1);
for i = 1:length(ContTime)
    if(strcmp(SL(1).Animal,'Ubi'))
        jitterT = 0;
    else
        jitterT = 0;
    end
    jitter = rand(1,length(ContTime{i}))*2*jitterT-jitterT; 
    hold on; scatter(ContTime{i}+jitter,ContTimeRT{i}-ContAll(i,1),'k.');
end
xlabel('Time of Last Stim (ms)'); ylabel('Reaction Time (ms)'); title('Contra Trials');
xlim([-2e4,0])

ax2 = subplot(2,1,2);
for i = 1:length(IpsiTime)
    if(strcmp(SL(1).Animal,'Ubi'))
        jitterT = 0;
    else
        jitterT = 0;
    end
    jitter = rand(1,length(ContTime{i}))*2*jitterT-jitterT; %100ms jitter
    hold on; scatter(IpsiTime{i},IpsiTimeRT{i}-IpsiAll(i,1),'k.');
end
xlabel('Time of Last Stim (ms)'); ylabel('Reaction Time (ms)'); title('Ipsi Trials');
xlim([-2e4,0])

linkaxes([ax1,ax2],'xy');

%% Binning 
CpreRT = nanmedian(ContAll(:,1));
IpreRT = nanmedian(IpsiAll(:,1));
Cstim = []; Istim = []; CRT = []; IRT = [];
for i = 1:length(ContTime)
    if(strcmp(SL(1).Animal,'Ubi') && i==14)
        continue;
    end
    
    if(strcmp(SL(1).Animal,'Ubi'))
        jitterT = 0;
    else
        jitterT = 0;
    end

    jitter = rand(1,length(ContTime{i}))*2*jitterT-jitterT; %100ms jitter
    Cstim = [Cstim,ContTime{i}+jitter];
    CRT = [CRT,ContTimeRT{i}];

    jitter = rand(1,length(IpsiTime{i}))*2*jitterT-jitterT; %100ms jitter
    Istim = [Istim,IpsiTime{i}+jitter];
    IRT = [IRT,IpsiTimeRT{i}];
end
% -100 to 100, 200-400, 400-600 and 700-900
figure;
% bins = -100:200:900;
bins = [-20000:2000:-2000;-18000:2000:0]';
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

RT = [CRT-CpreRT,IRT-IpreRT]; idx = [Cidx,Iidx];
bins = 1:((length(bins))*2);
nodata = setdiff(bins,idx);
if(~isempty(nodata))
    idx = [idx,nodata];
    RT = [RT,zeros(1,length(nodata))];
end

% Plot
positions = sort([1:length(unique(idx))/2,(1:length(unique(idx))/2)+1/3]);
boxplot(RT,idx, 'positions', positions);
hold on; xl = xlim;
line(xl,[0,0],'linestyle','--','color',[0,0,0])

% Set labels
labelpos = reshape(positions,2,length(positions)/2);
labelpos = labelpos(1,:)+diff(labelpos)/2;
set(gca,'xtick',labelpos)
set(gca,'xticklabel',labels)

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

% labels
title([SL(i).Animal,' \DeltaRT']); ylabel('RT (ms)'); xlabel('Stim Time (ms)')
%



