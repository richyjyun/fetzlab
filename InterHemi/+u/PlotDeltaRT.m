
clearvars -except SL UbiSL IgorSL KatoSL

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

%% Scatter Plot
CpreRT = cellfun(@(x) nanmean(x(100:end)),ContAll(:,1)); 
IpreRT = cellfun(@(x) nanmean(x(100:end)),IpsiAll(:,1));  
figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]);
for i = 1:length(StimAllC)
    if(strcmp(SL(1).Animal,'Ubi') && i==14)
        continue;
    end
    
    subplot(2,1,1);
    stimT = StimAllC{i};
    condRT = ContAll{i,2};
    
    jitter = 50;
    jitter = rand(1,length(stimT))*2*jitter-jitter; %100ms jitter
    
    scatter(stimT+jitter,condRT-CpreRT(i),'k.'); hold on;
    %         scatter(StimAll{i}+jitter,condRT,'k.'); hold on;
    xlabel('Time (ms)')
    ylabel('RT (ms)')
    title([SL(1).Animal,' Contra']);
    
    
    subplot(2,1,2);
    stimT = StimAllI{i};
    jitter = 50;
    jitter = rand(1,length(stimT))*2*jitter-jitter; %100ms jitter
    
    condRT = IpsiAll{i,2};
    %     scatter(stimT*ones(1,length(condRT))+jitter,condRT,'k.'); hold on;
    scatter(stimT+jitter,condRT-IpreRT(i),'k.'); hold on;
    xlabel('Time (ms)')
    ylabel('RT (ms)')
    title([SL(1).Animal,' Ipsi']);
    
end

CpreRT = cellfun(@nanmedian,ContAll(:,1)); CpreRT = mean(CpreRT);
IpreRT = cellfun(@nanmedian,IpsiAll(:,1)); IpreRT = mean(IpreRT);

subplot(2,1,1);% xlim([-100,800]);
xl = xlim; yl = ylim;
plot([CpreRT,CpreRT],yl,'k--');
plot(xl,[0,0],'k--');
subplot(2,1,2); % xlim([-100,800]);
xl = xlim; yl = ylim;
plot([CpreRT,CpreRT],yl,'k--');
plot(xl,[0,0],'k--');


%% Gut check to see % of stim times in each bin
BinCount = []; bins = [-100,200; 200,500; 500,800];

for i = 1:length(StimAll)
    stimT = StimAll{i};
    BinCount(1,i) = sum(stimT >= bins(1,1) & stimT < bins(1,2));
    BinCount(2,i) = sum(stimT >= bins(2,1) & stimT < bins(2,2));
    BinCount(3,i) = sum(stimT >= bins(3,1) & stimT < bins(3,2));
    
end

BinRatio = BinCount ./ sum(BinCount);
BinRatio = BinRatio';


%% Box Plot (cond post)
CRT = []; IRT = [];
CpreRT = cellfun(@nanmean,ContAll(:,1)); 
IpreRT = cellfun(@nanmean,IpsiAll(:,1)); 
for i = 1:length(StimAllC)
    if(strcmp(SL(1).Animal,'Ubi') && i==14)
        continue;
    end
    
    stimT = StimAllC{i};
    ind = stimT >= 500 & stimT < 800;
    condRT = ContAll{i,2}-CpreRT(i);
    CRT = [CRT,condRT(ind)];
    
    stimT = StimAllI{i};
    ind = stimT >= 500 & stimT < 800;
    condRT = IpsiAll{i,2}-IpreRT(i);
    IRT = [IRT,condRT(ind)];
end

CondRT = [CRT,IRT]; CondIdx = [ones(1,length(CRT)),2*ones(1,length(IRT))];

CRT = []; IRT = [];
for i = 1:length(StimAllC)
    if(strcmp(SL(1).Animal,'Ubi') && i==14)
        continue;
    end
    
    stimT = [StimAllC{i},StimAllI{i}];
    if(median(stimT) < 500 || median(stimT) > 800)
        continue;
    end
    
    condRT = ContAll{i,3}(50:end)-CpreRT(i);
    CRT = [CRT,condRT];
    
    condRT = IpsiAll{i,3}(50:end)-IpreRT(i);
    IRT = [IRT,condRT];
end

PostRT = [CRT,IRT]; PostIdx = [3*ones(1,length(CRT)),4*ones(1,length(IRT))];

RT = [CondRT,PostRT]; idx = [CondIdx,PostIdx];

% Plot
figure;
positions = sort([1:2,(1:2)+1/3]);
boxplot(RT,idx, 'Notch','on', 'positions', positions,'Whisker',0,'Symbol','');
hold on; xl = xlim;
line(xl,[0,0],'linestyle','--','color',[0,0,0])

% Color boxes
labels = {'Contra Cond','Ipsi Cond','Contra Post','Ipsi Post'};
color = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980]; % just using default matlab colors
color = repmat(color,2,1); color = flipud(color); % need to flip because call is backwards
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.5); % set face colors
end

% legend for colors
c = get(gca, 'Children');
[~,leg] = legend(c(1:2), 'Contra', 'Ipsi' );
PatchInLegend = findobj(leg, 'type', 'patch');
set(PatchInLegend, 'facealpha', 0.5);

ylim([-150,150])

% Plot vertical dashed line
x = positions(2) + (positions(3)-positions(2))/2;
yl = ylim;
plot([x,x],yl,'k--')

% Replotting to have box plot on top
boxplot(RT,idx, 'Notch','on', 'positions', positions,'Whisker',0,'Symbol','')

% Set labels
labelpos = positions;
set(gca,'xtick',labelpos)
set(gca,'xticklabel',labels)

% labels
title([SL(i).Animal,'\DeltaRT']); ylabel('RT (ms)'); xlabel('Stim Time (ms)')

ylim([-150,150])


%% Box Plot (post)
CRT = []; IRT = []; Cidx = []; Iidx = [];
CpreRT = cellfun(@nanmean,ContAll(:,1)); 
IpreRT = cellfun(@nanmean,IpsiAll(:,1)); 
bins = [-100,200; 200,500; 500,800];

for i = 1:length(StimAllC)
    if(strcmp(SL(1).Animal,'Ubi') && i==14)
        continue;
    end
    
    condRT = ContAll{i,3}(50:end)-CpreRT(i);
    idx = 0;
    stimT = mean(StimAllC{i});
    for b = 1:size(bins,1)
        if(stimT >= bins(b,1) && stimT < bins(b,2))
            idx = b*2-1;
            break;
        end
    end
    if(idx ~= 0)
        CRT = [CRT,condRT];
        Cidx = [Cidx,idx*ones(1,length(condRT))];
    end
    
    condRT = IpsiAll{i,3}(50:end)-IpreRT(i);
    idx = 0;
    stimT = mean(StimAllI{i});
    for b = 1:size(bins,1)
        if(stimT >= bins(b,1) && stimT < bins(b,2))
            idx = b*2;
            break;
        end
    end
    if(idx ~= 0)
        IRT = [IRT,condRT];
        Iidx = [Iidx,idx*ones(1,length(condRT))];
    end
end

RT = [CRT,IRT]; idx = [Cidx,Iidx]; %RT = [RT,0,0]; idx = [idx,5,6];

% Plot
figure;
positions = sort([1:size(bins,1),(1:size(bins,1))+1/3]);
boxplot(RT,idx, 'Notch','on', 'positions', positions,'Whisker',0,'Symbol','');
hold on; xl = xlim;
line(xl,[0,0],'linestyle','--','color',[0,0,0])

% Color boxes
labels = {'Prep','Move','Post'};
color = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980]; % just using default matlab colors
color = repmat(color,3,1); color = flipud(color); % need to flip because call is backwards
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.5); % set face colors
end

% legend for colors
c = get(gca, 'Children');
[~,leg] = legend(c(1:2), 'Contra', 'Ipsi' );
PatchInLegend = findobj(leg, 'type', 'patch');
set(PatchInLegend, 'facealpha', 0.5);

ylim([-150,150])

% Plot vertical dashed line
x = positions(2) + (positions(3)-positions(2))/2;
yl = ylim;

% Replotting to have box plot on top
boxplot(RT,idx, 'Notch','on', 'positions', positions,'Whisker',0,'Symbol','')

% Set labels
p1 = positions(1:2:end); p2 = positions(2:2:end);
labelpos = p1+(p2-p1)/2;
set(gca,'xtick',labelpos)
set(gca,'xticklabel',labels)

% labels
title([SL(i).Animal,' Post \DeltaRT']); ylabel('RT (ms)'); xlabel('Stim Time (ms)')

ylim([-150,150])


%% Box Plot
Cstim = []; Istim = []; CRT = []; IRT = [];
CpreRT = cellfun(@nanmean,ContAll(:,1)); 
IpreRT = cellfun(@nanmean,IpsiAll(:,1)); 
% bins = [-100,200; 200,500; 500,800];
bins = [0,250;250,500;500,750];

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
figure;
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
bins = 1:((length(bins))*2);
nodata = setdiff(bins,idx);
if(~isempty(nodata))
    idx = [idx,nodata];
    RT = [RT,zeros(1,length(nodata))];
end

% Plot
positions = sort([1:length(unique(idx))/2,(1:length(unique(idx))/2)+1/3]);
% subplot(2,1,1);
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
title([SL(i).Animal,'\DeltaRT']); ylabel('RT (ms)'); xlabel('Stim Time (ms)')

ylim([-150,150])


% %% Stats
% str = ''; sig = 0.05;
% 
% Cidx = [1,3,5];
% Iidx = [2,4,6];
% 
% Zstats = length(Cidx)+length(Iidx);
% for i = 1:length(Cidx)
%     Zstats(i) = signrank(RT(idx==Cidx(i)));
%     n = length(RT(idx==Cidx(i)));
%     stat = Zstats(i);
%     stat = stat*(nDays(Cidx(i)));
%     sstr = 'ns';
%     if(stat < sig)
%         sstr = 's';
%     end
%     str = sprintf('%sContra%d From Zero (%d) = %0.3e\t%s\n',str,i,n,stat,sstr);
% end
% 
% str = sprintf('%s\n',str);
% 
% for i = 1:length(Iidx)
%     Zstats(i) = signrank(RT(idx==Iidx(i)));
%     n = length(RT(idx==Iidx(i)));
%     stat = Zstats(i);
%     stat = stat*(nDays(Iidx(i)));
% 
%     sstr = 'ns';
%     if(stat < sig)
%         sstr = 's';
%     end
%     str = sprintf('%sIpsi%d From Zero (%d) = %0.3e \t%s\n',str,i,n,stat,sstr);
% end
% 
% str = sprintf('%s\n',str);
% 
% Compstats = length(Cidx);
% for i = 1:length(Cidx)
%     Compstats(i) = ranksum(RT(idx==Cidx(i)),RT(idx==Iidx(i)));
%     stat = Compstats(i);
%     stat = stat*(nDays(Cidx(i))+nDays(Iidx(i)));
% 
%     sstr = 'ns';
%     if(stat < sig)
%         sstr = 's';
%     end
%     str = sprintf('%sContra%d vs Ipsi%d = %0.3e\t%s\n',str,i,i,stat,sstr);
% end
% 
% str = sprintf('%s\n',str);
% 
% Cstats = nan(length(Cidx)-1);
% for i = 1:length(Cidx)
%     for j = i+1:length(Cidx)
%         Cstats(i,j) = ranksum(RT(idx==Cidx(i)),RT(idx==Cidx(j)));
%         stat = Cstats(i,j);
%         stat = stat*(nDays(Cidx(i))+nDays(Cidx(j)));
%         sstr = 'ns';
%         if(stat < sig)
%             sstr = 's';
%         end
%         str = sprintf('%sContra%d vs Contra%d = %0.3e\t%s\n',str,i,j,stat,sstr);
%     end
% end
% 
% str = sprintf('%s\n',str);
% 
% Istats = nan(length(Iidx)-1);
% for i = 1:length(Iidx)
%     for j = i+1:length(Iidx)
%         Istats(i,j) = ranksum(RT(idx==Iidx(i)),RT(idx==Iidx(j)));
%         stat = Istats(i,j);
%         stat = stat*(nDays(Iidx(i))+nDays(Iidx(j)));
% 
%         sstr = 'ns';
%         if(stat < sig)
%             sstr = 's';
%         end
%         str = sprintf('%sIpsi%d vs Ipsi%d = %0.3e   \t%s\n',str,i,j,stat,sstr);
%     end
% end
% str = sprintf('%s\n\n',str);
% 
% 
% fileID = fopen('C:\Users\richy.yun\Dropbox\repos\abogaard\efetz\RT manuscript\figures\Stats_Test.txt','a');
% fprintf(fileID,'%s Stats\n\n',SL(i).Animal);
% fprintf(fileID,str);
% fclose(fileID);


%% Stats 2. Run boxplot code first 
% % remove all that's outside the interquartile range
% RTs = []; idxs = [];
% for i = unique(idx)
%     temp = RT(idx==i);
%     q = quantile(temp,[0.25,0.75]);
%     temp(temp<q(1) | temp>q(2)) = [];
%     temp(isnan(temp)) = [];
%     RTs = [RTs,temp];
%     idxs = [idxs,i*ones(1,length(temp))];
% end


%% Stats
str = ''; sig = 0.05;

Cidx = [1,3,5];
Iidx = [2,4,6];

Zstats = length(Cidx)+length(Iidx);
for i = 1:length(Cidx)
    Zstats(i) = signrank(RT(idx==Cidx(i)));
    n = length(RT(idx==Cidx(i)));
    stat = Zstats(i);
    stat = stat*(nDays(Cidx(i)));
    sstr = 'ns';
    if(stat < sig)
        sstr = 's';
    end
    str = sprintf('%sContra%d From Zero (%d) = %0.3e   \t%s\n',str,i,n,stat,sstr);
end

str = sprintf('%s\n',str);

for i = 1:length(Iidx)
    Zstats(i) = signrank(RT(idx==Iidx(i)));
    n = length(RT(idx==Iidx(i)));
    stat = Zstats(i);
    stat = stat*(nDays(Iidx(i)));

    sstr = 'ns';
    if(stat < sig)
        sstr = 's';
    end
    str = sprintf('%sIpsi%d From Zero (%d) = %0.3e   \t%s\n',str,i,n,stat,sstr);
end

str = sprintf('%s\n',str);


naninds = isnan(RT);
RT(naninds) = []; idx(naninds) = [];


CvsI = [];
Cind = mod(idx,2)==1; Iind = mod(idx,2)==0;
CvsI(Cind) = 1; CvsI(Iind) = 2;

StimIdx = ceil(idx/2);

% [p,tbl,stats] = anovan(RT,{CvsI',StimIdx'},'varnames',{'CvsI','StimT'});%,'model','full');
% figure;
% results = multcompare(stats,'Dimension',[1,2],'CType','hsd');


% split up into contra vs ipsi? 
CRT = RT(CvsI == 1); CStim = StimIdx(CvsI == 1);
[p,tbl,stats] = kruskalwallis(CRT,CStim,'off');
results = multcompare(stats,'display','off');
counter = 1;
for i = 1:length(unique(CStim))
    for j = i+1:length(unique(CStim))
        p = results(counter,end);
        sstr = 'ns';
        if(p< sig)
            sstr = 's';
        end
        str = sprintf('%sContra%d vs Contra%d = %0.3e   \t%s\n',str,i,j,p,sstr);
        counter = counter+1;
    end
end


IRT = RT(CvsI == 2); IStim = StimIdx(CvsI == 2);
[p,tbl,stats] = kruskalwallis(IRT,IStim,'off');
results = multcompare(stats,'display','off');
str = sprintf('%s\n',str);

counter = 1;
for i = 1:length(unique(CStim))
    for j = i+1:length(unique(CStim))
        p = results(counter,end);
        sstr = 'ns';
        if(p< sig)
            sstr = 's';
        end
        str = sprintf('%sIpsi%d vs Ipsi%d = %0.3e   \t%s\n',str,i,j,p,sstr);
        counter = counter+1;
    end
end

subplot(2,1,2);
text(0,0,str,'Verticalalignment','bottom');
axis off;
