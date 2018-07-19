% 
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbiFinal.mat');
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbiTime.mat');

clearvars -except SL SLNeuroTime

dates = extractfield(SL,'Date');
Cont = struct('C',[],'I',[]); Ipsi = struct('C',[],'I',[]);
StimAll = struct('C',[],'I',[]);

CondAll = {};

% CHANGE PER ANIMAL
% Ubi
Lchns = [31:32]; Rchns = [6,7];
% Igor
% Lchns = [2:3]; Rchns = [10,11];
%Kato 
% Lchns = [9,10]; Rchns = [13,14];

for n = 1:length(SLNeuroTime)
    i = find(strcmp(dates,SLNeuroTime(n).Date)); trig = SL(i).trig1;
    D = str2num(SL(i).Date);
    
    if (~isempty(SL(i).Bad) || strcmp(SL(i).Condition, 'Control')...
            || strcmp(SL(i).Condition, 'NaN') || strcmp(SL(i).Condition, 'nostim')...
            || isempty(SL(i).trig1) || strcmp(SL(i).Condition(end),'R'))
        continue;
    end
    
    
    if(length(SL(i).Condition)<6 || ~strcmp(SL(i).Condition(1:6),'Contra'))
        continue;
    end
    
    %     % 100% stim for Ubi
    %     if(strcmp(SL(i).Animal,'Ubi') && (str2num(SL(i).Date) >= 20170226 && str2num(SL(i).Date) <= 20170306))
    %         continue;
    %     end
    
    
    % Remove stim channels
    Lchn = Lchns; Rchn = Rchns;
    if(strcmp(SL(i).Animal,'Ubi'))
        StimSites = strsplit(SL(i).Stim_Loc,'/');
        StimSites{2} = [StimSites{1}(1:4),StimSites{2}];
    else
        StimSites = SL(i).Stim_Loc;
    end
    for s = 1:length(StimSites)
        site = find(cell2mat(cellfun(@(x) strcmp(x,StimSites{s}), SLNeuroTime(n).chnm,'UniformOutput',0)));
        if(ismember(site,Lchn))
            Lchn(find(Lchn==site)) = [];
        elseif(ismember(site,Rchn))
            Rchn(find(Rchn==site)) = [];
        end
    end
    
    
    if(strcmp(SL(i).StimHemi,'L'))
        temp = SLNeuroTime(n).Rbeta(:,Rchn); temp = cellfun(@transpose,temp,'Uniformoutput',false);
        ContRT.C = cellfun(@(col) vertcat(col{:}), num2cell(temp, 1), 'UniformOutput', false);
        temp = SLNeuroTime(n).Rbeta(:,Lchn); temp = cellfun(@transpose,temp,'Uniformoutput',false);
        ContRT.I = cellfun(@(col) vertcat(col{:}), num2cell(temp, 1), 'UniformOutput', false);
        temp = SLNeuroTime(n).Lbeta(:,Rchn); temp = cellfun(@transpose,temp,'Uniformoutput',false);
        IpsiRT.C = cellfun(@(col) vertcat(col{:}), num2cell(temp, 1), 'UniformOutput', false);
        temp = SLNeuroTime(n).Lbeta(:,Lchn); temp = cellfun(@transpose,temp,'Uniformoutput',false);
        IpsiRT.I= cellfun(@(col) vertcat(col{:}), num2cell(temp, 1), 'UniformOutput', false);
        ContT = SL(i).righttrials;
        IpsiT = SL(i).lefttrials;
        CRT = SL(i).rts_r;
        IRT = SL(i).rts_l;
    else
        temp = SLNeuroTime(n).Lbeta(:,Rchn); temp = cellfun(@transpose,temp,'Uniformoutput',false);
        ContRT.C = cellfun(@(col) vertcat(col{:}), num2cell(temp, 1), 'UniformOutput', false);
        temp = SLNeuroTime(n).Lbeta(:,Lchn); temp = cellfun(@transpose,temp,'Uniformoutput',false);
        ContRT.I = cellfun(@(col) vertcat(col{:}), num2cell(temp, 1), 'UniformOutput', false);
        temp = SLNeuroTime(n).Rbeta(:,Rchn); temp = cellfun(@transpose,temp,'Uniformoutput',false);
        IpsiRT.C = cellfun(@(col) vertcat(col{:}), num2cell(temp, 1), 'UniformOutput', false);
        temp = SLNeuroTime(n).Rbeta(:,Lchn); temp = cellfun(@transpose,temp,'Uniformoutput',false);
        IpsiRT.I= cellfun(@(col) vertcat(col{:}), num2cell(temp, 1), 'UniformOutput', false);
        ContT = SL(i).lefttrials;
        IpsiT = SL(i).righttrials;
        CRT = SL(i).rts_l;
        IRT = SL(i).rts_r;
    end
    
    % get all trials that had stimulation
    if(strcmp(SL(i).Animal,'Ubi')) %Kato trigger times are odd, so just assume all stim
        [StimT,idx,stimT] = u.getStimTrials(SL(i).trig1,ContT,50);
    else
        [StimT,idx,stimT] = u.getStimTrials(SL(i).trig1,ContT,500);
    end
    trig = SL(i).trig1(idx);
    
    % put trials in order
    Trials = [ContT(:,1);IpsiT(:,1)];
    Stim = zeros(1,length(Trials));
    Stim(StimT) = 1;
    Labels = [ones(1,length(ContT)),zeros(1,length(IpsiT))]; % 1 is contra, 0 is ipsi
    RT = [CRT;IRT];
    
    [Trials,order] = sort(Trials);
    Stim = Stim(order);
    Labels = Labels(order);
    RT = RT(order);

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
        
        if(isnan(RT(t)) || Stim(t))
            continue;
        end
        
        if(Labels(t))
            stimind = sum(Stim(1:t-since));
            StimC(end+1) = stimT(stimind);
        else
            stimind = sum(Stim(1:t-since));
            StimI(end+1) = stimT(stimind);
        end
    end
   
    % save stim timings
    StimAll.C{end+1} = StimC;
    StimAll.I{end+1} = StimI;
    
    % find median pre
    CPreInd = find(Labels(1:start-1)); IPreInd = find(~Labels(1:start-1));
    CPre.C = cellfun(@(x) nanmedian(x(CPreInd,:)),ContRT.C,'UniformOutput',false);
    CPre.I = cellfun(@(x) nanmedian(x(CPreInd,:)),ContRT.I,'UniformOutput',false);
    IPre.C = cellfun(@(x) nanmedian(x(IPreInd,:)),IpsiRT.C,'UniformOutput',false);
    IPre.I = cellfun(@(x) nanmedian(x(IPreInd,:)),IpsiRT.I,'UniformOutput',false);

    % save appropriate values
    Cind = find(ContT(:,1)>trig(end),1); Iind = find(IpsiT(:,1)>trig(end),1);
    Cont.C{end+1} = cellfun(@(x,y) x(Cind:end,:)-y,ContRT.C,CPre.C,'UniformOutput',false);
    Cont.I{end+1} = cellfun(@(x,y) x(Cind:end,:)-y,ContRT.I,CPre.I,'UniformOutput',false);
    Ipsi.C{end+1} = cellfun(@(x,y) x(Iind:end,:)-y,IpsiRT.C,IPre.C,'UniformOutput',false);
    Ipsi.I{end+1} = cellfun(@(x,y) x(Iind:end,:)-y,IpsiRT.I,IPre.I,'UniformOutput',false);
    
%     Cont.C{end+1} = cellfun(@(x) x(Cind,:),ContRT.C,'UniformOutput',false);
%     Cont.I{end+1} = cellfun(@(x) x(Cind,:),ContRT.I,'UniformOutput',false);
%     Ipsi.C{end+1} = cellfun(@(x) x(Iind,:),IpsiRT.C,'UniformOutput',false);
%     Ipsi.I{end+1} = cellfun(@(x) x(Iind,:),IpsiRT.I,'UniformOutput',false);
    
    
    % save conditions
    CondAll{end+1,1} = SL(i).Condition;
    CondAll{end,2} = SL(i).Stim_Delay;
    CondAll{end,3} = SL(i).StimHemi;
    
end


%% Box plots
bins = [-100,200; 200,500; 500,800];
Beta = struct('C',[],'I',[]);  %for each hemisphere
Group = struct('C',[],'I',[]);

for d = 1:length(Cont.C)
    bcounter = 1;
    
    for t = 1:2 % type of trial
        if(t==1)
            Trial = Cont;
            side = 'C';
        else
            Trial = Ipsi;
            side = 'I';
        end
        for i = 1:3 %epoch (pre/move/post)
            
            stimT = nanmean(StimAll.(side){d});
            bshift = 0;
            for b = 1:size(bins,1)
                if(stimT >= bins(b,1) && stimT < bins(b,2))
                    bshift = b;
                end
            end
            if(bshift == 0)
                continue;
            end
%             if(bshift ==1)
%                 keyboard;
%             end
            temp = cellfun(@(x) x(:,i),Trial.C{d},'UniformOutput',false);
            temp = [temp{:}]; temp = temp(:);
            Beta.C = [Beta.C;temp];
            Group.C = [Group.C;((bshift-1)*6+bcounter)*ones(length(temp),1)];
            
            temp = cellfun(@(x) x(:,i),Trial.I{d},'UniformOutput',false);
            temp = [temp{:}]; temp = temp(:);
            Beta.I = [Beta.I;temp];
            Group.I = [Group.I;((bshift-1)*6+bcounter)*ones(length(temp),1)];
            bcounter = bcounter+1;
        end
    end
end

% set variables
positions = 1:20; positions([7,14]) = [];
color1 = [0, 0.4470, 0.7410]; color2 = [0.8500, 0.3250, 0.0980]; % just using default matlab colors
color = [repmat(color1,3,1); repmat(color2,3,1)];
color = repmat(color,3,1); color = flipud(color); 
xt = [3.5,10.5,17.5];
xtl = {'Pre','Move','Post'};

% Plot
figure; subplot(2,1,1);
boxplot(Beta.C,Group.C,'Notch','on', 'positions', positions,'Whisker',0,'Symbol','');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.5); % set face colors
end
title('Contra Exp, Beta, Contra Hemi'); 
hold on; xl = xlim;
plot(xl,[0,0],'k--')
% xticks(xt); xticklabels(xtl);

c = get(gca, 'Children');
[~,leg] = legend(c([2,5]), 'Contra', 'Ipsi' );
PatchInLegend = findobj(leg, 'type', 'patch');
set(PatchInLegend, 'facealpha', 0.5);
boxplot(Beta.C,Group.C,'Notch','on', 'positions', positions,'Whisker',0,'Symbol','');
ylim([-3,3]);

subplot(2,1,2);
boxplot(Beta.I,Group.I,'Notch','on', 'positions', positions,'Whisker',0,'Symbol','');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.5); % set face colors
end
title('Contra Exp, Beta, Ipsi Hemi'); 
hold on; xl = xlim;
plot(xl,[0,0],'k--')
boxplot(Beta.I,Group.I,'Notch','on', 'positions', positions,'Whisker',0,'Symbol','');
ylim([-3,3])

% xticks(xt); xticklabels(xtl);
