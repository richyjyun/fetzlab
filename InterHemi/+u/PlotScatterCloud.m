%% Scatter cloud of RT relative to stim times

% clear;
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbiFinal.mat');
% UbiSL = SL;
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgorFinal.mat');
% IgorSL = SL;
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKatoFinal.mat');
% KatoSL = SL;

clearvars -except SL UbiSL IgorSL KatoSL

ContAll = {};
IpsiAll = {};
StimAll = {};
CondAll = {};
StimAllI = {};

for i = 1:length(SL)
    if (~isempty(SL(i).Bad) || strcmp(SL(i).Condition, 'Control')...
            || strcmp(SL(i).Condition, 'NaN') || strcmp(SL(i).Condition, 'nostim')...
            || isempty(SL(i).trig1) || strcmp(SL(i).Condition(end),'R'))
        continue;
    end
    
    if(length(SL(i).Condition)<6 || ~strcmp(SL(i).Condition(1:6),'Contra'))
        continue;
    end
    
    % 100% stim for Ubi
    if(strcmp(SL(i).Animal,'Ubi') && (str2num(SL(i).Date) < 20170226 || str2num(SL(i).Date)> 20170306))
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
        [StimT,idx] = u.getStimTrials(SL(i).trig1,ContT);
    else
        CondStart = find(ContT(:,1)<SL(i).trig1(1),1,'last');
        CondEnd = find(ContT(:,1)>SL(i).trig1(end),1);
        StimT = CondStart:(CondEnd-1);
    end
    trig = SL(i).trig1;
    
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
    if(~strcmp(SL(i).Animal,'Ubi'))
        stimT = trig;
    else
        stimT = trig+str2num(SL(i).Stim_Delay); % have to add delay for Ubi
    end
    
    % get all contra and ipsi within 1 trial of stim
    % get ipsi trial stim time from GO
    start = find(Stim,1); finish = find(Stim,1,'last')+1;
    CondContRT = []; CondIpsiRT = []; StimI = []; StimC=[];
    for t = start:finish
        if(~(Stim(t-1))) %only looking at trials immediately after stim
            continue;
        end
        
        if(Labels(t))
            CondContRT(end+1) = RT(t);
        else
            CondIpsiRT(end+1) = RT(t);
            % Find most recent stim for ipsi trials
            ind = find(stimT < Trials(t,1),1,'last');
            if(isempty(ind))
                StimI(end+1) = nan;
            else
                StimI(end+1) = stimT(ind) - Trials(t,1);
            end
        end
    end
    
    % Save RT
    ContAll{end+1,1} = PreContRT;
    ContAll{end,2} = CondContRT;
    ContAll{end,3} = PostContRT;
    
    IpsiAll{end+1,1} = PreIpsiRT;
    IpsiAll{end,2} = CondIpsiRT;
    IpsiAll{end,3} = PostIpsiRT;
    
    % Find timing of stim relative to start of trials
    Stim = [];
    
    for s = 1:length(stimT)
        % if igor, exp that are not movement triggered have trigs at stim
        % time, do not need to add delay
        
        ind = find(ContT(:,1) < (stimT(s)+50),1,'last');
        if(isnan(ContRT(ind)))
            Stim(s) = nan;
        else
            Stim(s) = (stimT(s)-ContT(ind,1));
        end
        ind = find(IpsiT(:,2) > stimT(s) ,1);
    end
    
    StimAll{end+1} = Stim;
    StimAllI{end+1} = StimI;
    
    %     if(length(StimAll) == 6)
    %         keyboard
    %     end
    %
    CondAll{end+1,1} = SL(i).Condition;
    CondAll{end,2} = SL(i).Stim_Delay;
    CondAll{end,3} = SL(i).StimHemi;
    
end

%
%% Scatter Plot
figure;
for i = 1:length(StimAllC)
    if(strcmp(SL(1).Animal,'Ubi') && i==14)
        continue;
    end
    
    subplot(2,1,1);
    stimT = StimAllC{i};
    condRT = ContAll{i,2};
    
    jitter = 50;
    jitter = rand(1,length(stimT))*2*jitter-jitter; %100ms jitter
    
    scatter(stimT+jitter,condRT,'k.'); hold on;
    %         scatter(StimAll{i}+jitter,condRT,'k.'); hold on;
    
    subplot(2,1,2);
    stimT = StimAllI{i};
    jitter = 50;
    jitter = rand(1,length(stimT))*2*jitter-jitter; %100ms jitter
    
    condRT = IpsiAll{i,2};
    %     scatter(stimT*ones(1,length(condRT))+jitter,condRT,'k.'); hold on;
    scatter(stimT+jitter,condRT,'k.'); hold on;
end

CpreRT = cellfun(@nanmedian,ContAll(:,1)); CpreRT = mean(CpreRT);
IpreRT = cellfun(@nanmedian,IpsiAll(:,1)); IpreRT = mean(IpreRT);
% subplot(2,1,1); title('Contra RT'); %xlim([-100,800])
% xl = xlim; yl = ylim;
% plot([CpreRT,CpreRT],yl,'r','linewidth',2);
% plot(xl,[CpreRT,CpreRT],'r','linewidth',2);
% xlabel('Stim Time from GO (ms)'); ylabel('RT (ms)');
% subplot(2,1,2); title('Ipsi RT');  %xlim([-100,800])
% xl = xlim; yl = ylim;
% plot([CpreRT,CpreRT],yl,'r','linewidth',2);
% plot(xl,[IpreRT,IpreRT],'r','linewidth',2);
% xlabel('Stim Time from GO (ms)'); ylabel('RT (ms)');

% % See how many good RT there are per experiment
% good = cellfun(@(x) sum(~isnan(x)),StimAll);
% good = good./cellfun(@length,StimAll);
% good = good';



% figure;
% for i = 1:length(StimAllI)
%     jitterT = 10;
%     jitter = rand(1,length(StimAllI{i}))*2*jitterT-jitterT; %100ms jitter
%
%     scatter(StimAllI{i}+jitter, IpsiAll{i,2},'k.'); hold on;
% end
% xl = xlim;
% plot(xl,[IpreRT,IpreRT],'r','linewidth',2);

%% Bar Plot
Cstim = []; Istim = []; CRT = []; IRT = [];
CpreRT = cellfun(@nanmedian,ContAll(:,1)); CpreRT = mean(CpreRT);
IpreRT = cellfun(@nanmedian,IpsiAll(:,1)); IpreRT = mean(IpreRT);
for i = 1:length(StimAll)
    if(strcmp(SL(1).Animal,'Ubi') && i==14)
        continue;
    end
    
    stimT = nanmedian(StimAll{i});
    if(strcmp(SL(1).Animal,'Ubi'))
        jitterT = 0;
    else
        jitterT = 0;
    end
    
    condRT = ContAll{i,2};
    jitter = rand(1,length(condRT))*2*jitterT-jitterT; %100ms jitter
    Cstim = [Cstim,stimT+jitter];
    CRT = [CRT,condRT];
    
    stimT = nanmedian(StimAll{i});
    condRT = IpsiAll{i,2};
    jitter = rand(1,length(condRT))*2*jitterT-jitterT; %100ms jitter
    Istim = [Istim,stimT+jitter];
    IRT = [IRT,condRT];
end
% -100 to 100, 200-400, 400-600 and 700-900
% bins = -100:200:900;
bins = [-100,200; 200,500; 500,800];
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
figure; Bars = [];
nidx = length(unique(idx)); nidx = 1:(nidx/2);
positions = sort([((nidx)-1)*3+1,((nidx)-1)*3+2]);
color = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980]; % just using default matlab colors
for i = 1:length(positions)
    cind = 1;
    if mod(i,2) == 0
        cind = 2;
    end
    data = RT(idx==i);
    Bars(i) = bar(positions(i),nanmedian(data),'FaceColor',color(cind,:)); hold on;
    iqr = diff(prctile(data, [25 75]));
    err = 1.57 * iqr/sqrt(length(data));
    errorbar(positions(i),nanmedian(data),err,'k','Linewidth',1.5);
end

% legend for colors
legend( Bars(1:2),'Contra', 'Ipsi' );

% Set labels
labelpos = reshape(positions,2,length(positions)/2);
labelpos = labelpos(1,:)+diff(labelpos)/2;
set(gca,'xtick',labelpos)
set(gca,'xticklabel',labels)

% labels
title([SL(i).Animal,'\DeltaRT']); ylabel('RT (ms)'); xlabel('Stim Time (ms)')


%% Box Plot
Cstim = []; Istim = []; CRT = []; IRT = [];
CpreRT = cellfun(@nanmedian,ContAll(:,1)); 
IpreRT = cellfun(@nanmedian,IpsiAll(:,1)); 
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
    CRT = [CRT,condRT-CpreRT(i)];
    
    stimT = StimAllI{i};
    condRT = IpsiAll{i,2};
    Istim = [Istim,stimT];
    IRT = [IRT,condRT-IpreRT(i)];
end
% -100 to 100, 200-400, 400-600 and 700-900
figure;
% bins = -100:200:900;
bins = [-100,200; 200,500; 500,800];
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

%% Stats
str = ''; sig = 0.05;

Cidx = unique(Cidx);
Iidx = unique(Iidx);

Zstats = length(Cidx)+length(Iidx);
for i = 1:length(Cidx)
    [~,Zstats(i)] = ttest(RT(idx==Cidx(i)));
    n = length(RT(idx==Cidx(i)));
    stat = Zstats(i);
    sstr = 'ns';
    if(stat < sig)
        sstr = 's';
    end
    str = sprintf('%sContra%d From Zero (%d) = %0.3e\t%s\n',str,i,n,stat,sstr);
end

str = sprintf('%s\n',str);

for i = 1:length(Iidx)
    [~,Zstats(i)] = ttest(RT(idx==Iidx(i)));
    n = length(RT(idx==Iidx(i)));
    stat = Zstats(i);
    sstr = 'ns';
    if(stat < sig)
        sstr = 's';
    end
    str = sprintf('%sIpsi%d From Zero (%d) = %0.3e \t%s\n',str,i,n,stat,sstr);
end

str = sprintf('%s\n',str);

Compstats = length(Cidx);
for i = 1:length(Cidx)
    [~,Compstats(i)] = ttest2(RT(idx==Cidx(i)),RT(idx==Iidx(i)));
    stat = Compstats(i);
    sstr = 'ns';
    if(stat < sig)
        sstr = 's';
    end
    str = sprintf('%sContra%d vs Ipsi%d = %0.3e\t%s\n',str,i,i,stat,sstr);
end

str = sprintf('%s\n',str);

Cstats = nan(length(Cidx)-1);
for i = 1:length(Cidx)
    for j = i+1:length(Cidx)
        [~,Cstats(i,j)] = ttest2(RT(idx==Cidx(i)),RT(idx==Cidx(j)));
        stat = Cstats(i,j);
        sstr = 'ns';
        if(stat < sig)
            sstr = 's';
        end
        str = sprintf('%sContra%d vs Contra%d = %0.3e\t%s\n',str,i,j,stat,sstr);
    end
end

str = sprintf('%s\n',str);

Istats = nan(length(Iidx)-1);
for i = 1:length(Iidx)
    for j = i+1:length(Iidx)
        [~,Istats(i,j)] = ttest2(RT(idx==Iidx(i)),RT(idx==Iidx(j)));
        stat = Istats(i,j);
        sstr = 'ns';
        if(stat < sig)
            sstr = 's';
        end
        str = sprintf('%sIpsi%d vs Ipsi%d = %0.3e   \t%s\n',str,i,j,stat,sstr);
    end
end
str = sprintf('%s\n\n',str);


fileID = fopen('C:\Users\richy.yun\Dropbox\repos\abogaard\efetz\RT manuscript\figures\Stats_TTest.txt','a');
fprintf(fileID,'%s Stats\n\n',SL(i).Animal);
fprintf(fileID,str);
fclose(fileID);

% %% Mean +- Error
% figure; %subplot(3,1,3)
% Data = zeros(1,length(unique(idx)));
% positions = 1:(length(unique(idx))/2);
% positions = sort(repmat(positions,1,2));
% for i = 1:length(unique(idx))
%     if mod(i,2) == 0
%         c = color(1,:);
%     else
%         c = color(2,:);
%     end
%     data = RT(idx==i);
%     Data(i) = nanmean(data);
%     N = sum(~isnan(data));
%     error = nanstd(data)/sqrt(N);
%     errorbar(positions(i),nanmean(data),error,'o','color',c,'linewidth',1.5);
%     hold on;
%     text(positions(i),nanmean(data)+10,num2str(N));
% end
% plot(positions(1:2:length(Data)),Data(1:2:length(Data)),'color',color(2,:),'linewidth',1.5);
% plot(positions(2:2:length(Data)),Data(2:2:length(Data)),'color',color(1,:),'linewidth',1.5);
% set(gca,'xtick',1:(length(unique(idx))/2))
% set(gca,'xticklabel',labels);
% xlim([positions(1)-0.25,positions(end)+0.25])
% legend('Contra', 'Ipsi' );
% title([SL(i).Animal,'\DeltaRT']); ylabel('RT (ms)'); xlabel('Stim Time (ms)')

