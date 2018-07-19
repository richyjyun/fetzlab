% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgorFinal.mat');

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
    
    % get pre and post stim RT
    PreContInd = find(ContT(:,2)<trig(1),1,'last');
    PostContInd = find(ContT(:,1)>trig(end),1);
    PreIpsiInd = find(IpsiT(:,2)<trig(1),1,'last');
    PostIpsiInd = find(IpsiT(:,1)>trig(end),1);
    
    PreContRT = ContRT(1:PreContInd);
    CondContRT = ContRT(PreContInd+1:PostContInd-1);
    PostContRT = ContRT(PostContInd:end);
    PreIpsiRT = ContRT(1:PreIpsiInd);
    CondIpsiRT = ContRT(PreIpsiInd+1:PostIpsiInd-1);
    PostIpsiRT = ContRT(PostIpsiInd:end);
    
    
    % Save RT
    ContAll{end+1,1} = PreContRT;
    ContAll{end,2} = CondContRT;
    ContAll{end,3} = PostContRT;
    
    IpsiAll{end+1,1} = PreIpsiRT;
    IpsiAll{end,2} = CondIpsiRT;
    IpsiAll{end,3} = PostIpsiRT;
    
    
    CondAll{end+1,1} = SL(i).Condition;
    CondAll{end,2} = SL(i).Stim_Delay;
    CondAll{end,3} = SL(i).StimHemi;
    CondAll{end,4} = SL(i).Date;
end


%% Boxplots (Cond)
RT = []; Groups = []; counter = 1;
for i = 1:size(ContAll,1)
    RT = [RT;ContAll{i,2} - nanmean(ContAll{i,1})];
    Groups = [Groups;counter*ones(length(ContAll{i,2}),1)];
    counter = counter+1;
    
    RT = [RT;IpsiAll{i,2} - nanmean(IpsiAll{i,1})];
    Groups = [Groups;counter*ones(length(IpsiAll{i,2}),1)];
    counter = counter+1;
end

positions = unique(Groups);
positions = 1:((length(positions)/2)*3); 
positions(3:3:end) = [];
figure;
boxplot(RT,Groups, 'Notch','on', 'positions', positions,'Whisker',0,'Symbol','');
hold on; xl = xlim;
line(xl,[0,0],'linestyle','--','color',[0,0,0])

% Color boxes
color = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980]; % just using default matlab colors
color = repmat(color,length(positions)/2,1); color = flipud(color); % need to flip because call is backwards
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.5); % set face colors
end

% plot again 
boxplot(RT,Groups, 'Notch','on', 'positions', positions,'Whisker',0,'Symbol','');

% set labels
labels = 1:(length(positions)/2);
xt = positions(1:2:end) + 0.5;
xticks(xt); xticklabels(labels);

% labels
title([SL(i).Animal,' Tonic Cond \DeltaRT']); ylabel('RT (ms)'); xlabel('Experiment Days')

ylim([-150,150])

%% Boxplots (Post)
RT = []; Groups = []; counter = 1;
for i = 1:size(ContAll,1)
    RT = [RT;ContAll{i,3} - nanmean(ContAll{i,1})];
    Groups = [Groups;counter*ones(length(ContAll{i,3}),1)];
    counter = counter+1;
    
    RT = [RT;IpsiAll{i,3} - nanmean(IpsiAll{i,1})];
    Groups = [Groups;counter*ones(length(IpsiAll{i,3}),1)];
    counter = counter+1;
end

positions = unique(Groups);
positions = 1:((length(positions)/2)*3); 
positions(3:3:end) = [];
figure;
boxplot(RT,Groups, 'Notch','on', 'positions', positions,'Whisker',0,'Symbol','');
hold on; xl = xlim;
line(xl,[0,0],'linestyle','--','color',[0,0,0])

% Color boxes
color = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980]; % just using default matlab colors
color = repmat(color,length(positions)/2,1); color = flipud(color); % need to flip because call is backwards
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.5); % set face colors
end

% plot again 
boxplot(RT,Groups, 'Notch','on', 'positions', positions,'Whisker',0,'Symbol','');

% set labels
labels = 1:(length(positions)/2);
xt = positions(1:2:end) + 0.5;
xticks(xt); xticklabels(labels);

% labels
title([SL(i).Animal,' Tonic Post \DeltaRT']); ylabel('RT (ms)'); xlabel('Experiment Days')

ylim([-150,150])

