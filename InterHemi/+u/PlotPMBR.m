% SLNeuroTime = u.getLFPTime(SL);
% save('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbiNeuroTime','SLNeuroTime','-v7.3');
% 
% SL = IgorSL;
% SLNeuroTime = u.getLFPTime(SL);
% save('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgorNeuroTime','SLNeuroTime','-v7.3');
% 
% SL = KatoSL;
% SLNeuroTime = u.getLFPTime(SL);
% save('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKatoNeuroTime','SLNeuroTime','-v7.3');

%% 'RT' IS PMBR MAGNITUDE. JUST EASIER TO KEEP FROM BEFORE

SLNeuroTime = Original;

dates = extractfield(SL,'Date');
Cont = struct('C',[],'I',[]); Ipsi = struct('C',[],'I',[]); 
PreContRT = []; PreIpsiRT = []; PostContRT = []; PostIpsiRT = [];

StimAll = {};
CondAll = {};
StimAllI = {};

% CHANGE PER ANIMAL
% Ubi
Lchn = [17:32]; Rchn = [1:16];
% Igor
% Lchn = [1:9]; Rchn = [10,11];

ContRT = []; IpsiRT = [];

for n = 1:length(SLNeuroTime)
    
    i = find(strcmp(dates,SLNeuroTime(n).Date)); trig = SL(i).trig1;
    D = str2num(SL(i).Date);
    if(~strcmp(SL(i).Condition(1:6),'Contra'))
        continue;
    end
    
    if (~isempty(SL(i).Bad) || strcmp(SL(i).Condition, 'Control')...
            || strcmp(SL(i).Condition, 'NaN') || strcmp(SL(i).Condition, 'nostim')...
            || isempty(SL(i).trig1) || strcmp(SL(i).Condition(end),'R'))
        continue;
    end
    
    if(length(SL(i).Condition)<6 || ~strcmp(SL(i).Condition(1:6),'Contra'))
        continue;
    end
    
    %         % 100% stim for Ubi
    %         if(str2num(SL(i).Date) < 20170226 || str2num(SL(i).Date)> 20170306)
    %             continue;
    %         end
    
    % Remove stim channels
    if(strcmp(SL(i).Animal,'Ubi'))
        StimSites = strsplit(SL(i).Stim_Loc,'/');
        StimSites{2} = [StimSites{1}(1:4),StimSites{2}];
    else
        StimSites = SL(i).Stim_Loc;
    end
    for s = 1:length(StimSites)
        site = find(cell2mat(cellfun(@(x) strcmp(x,StimSites{s}), SLNeuroTime(n).chnm,'UniformOutput',0)));
        if(~isempty(site))
            sz = cellfun(@length,SLNeuroTime(n).Rbeta(:,site));
            SLNeuroTime(n).Rbeta(:,site) = {nan(1,sz(1));nan(1,sz(2));nan(1,sz(3))};
            sz = cellfun(@length,SLNeuroTime(n).Lbeta(:,site));
            SLNeuroTime(n).Lbeta(:,site) = {nan(1,sz(1));nan(1,sz(2));nan(1,sz(3))};
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
        [StimT,idx] = u.getStimTrials(SL(i).trig1,ContT);
    else
        CondStart = find(ContT(:,1)<SL(i).trig1(1),1,'last');
        CondEnd = find(ContT(:,1)>SL(i).trig1(end),1);
        StimT = CondStart:(CondEnd-1);
    end
    trig = SL(i).trig1;
    
    % get pre and post stim 
    PreContInd = find(ContT(:,2)<trig(1),1,'last');
    PreContRT.C = cellfun(@(x) x(1:PreContInd,:),ContRT.C,'Uniformoutput',false);
    PreContRT.I = cellfun(@(x) x(1:PreContInd,:),ContRT.I,'Uniformoutput',false);
    PreIpsiInd = find(IpsiT(:,2)<trig(1),1,'last');
    PreIpsiRT.C = cellfun(@(x) x(1:PreIpsiInd,:),IpsiRT.C,'Uniformoutput',false);
    PreIpsiRT.I = cellfun(@(x) x(1:PreIpsiInd,:),IpsiRT.I,'Uniformoutput',false);
    
    PostContInd = find(ContT(:,1)>trig(end),1);
    PostContRT.C = cellfun(@(x) x(PostContInd:end,:,:),ContRT.C,'Uniformoutput',false);
    PostContRT.I = cellfun(@(x) x(PostContInd:end,:,:),ContRT.I,'Uniformoutput',false);
    PostIpsiInd = find(IpsiT(:,1)>trig(end),1);
    PostIpsiRT.C = cellfun(@(x) x(PostIpsiInd:end,:,:),IpsiRT.C,'Uniformoutput',false);
    PostIpsiRT.I = cellfun(@(x) x(PostIpsiInd:end,:,:),IpsiRT.C,'Uniformoutput',false);
    
    % put trials in order
    Trials = [ContT(:,1);IpsiT(:,1)];
    
    
    % Put PMBR in order
    CHemi = [ContRT.C;IpsiRT.C];
    CHemi = cellfun(@(col) vertcat(col{:}), num2cell(CHemi, 1), 'UniformOutput', false);
    CHemi = cellfun(@(col) horzcat(col{:}), num2cell(CHemi, 2), 'UniformOutput', false);
    IHemi = [ContRT.I;IpsiRT.I];
    IHemi = cellfun(@(col) vertcat(col{:}), num2cell(IHemi, 1), 'UniformOutput', false);
    IHemi = cellfun(@(col) horzcat(col{:}), num2cell(IHemi, 2), 'UniformOutput', false);
    
    Stim = zeros(1,length(Trials));
    Stim(StimT) = 1;
    Labels = [ones(1,length(ContT)),zeros(1,length(IpsiT))]; % 1 is contra, 0 is ipsi
    
    [Trials,order] = sort(Trials);
    CHemi = CHemi{1}(order,:);
    IHemi = IHemi{1}(order,:);

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
    CTrial = struct('C',[],'I',[]); ITrial = struct('C',[],'I',[]);
    for t = start:finish
        if(~Stim(t-1))%(~(Stim(t) || Stim(t-1)))
            continue;
        end
        
        if(Labels(t)) % if contra trial
            CTrial.C(end+1,:) = CHemi(t,:);
            CTrial.I(end+1,:) = IHemi(t,:);
        else
            ITrial.C(end+1,:) = CHemi(t,:);
            ITrial.I(end+1,:) = IHemi(t,:);
        end
    end
    
    % Save All
    Cont.C{end+1} = CTrial.C - cellfun(@nanmedian,PreContRT.C);
    Cont.I{end+1} = CTrial.I - cellfun(@nanmedian,PreContRT.I);
    Ipsi.C{end+1} = ITrial.C - cellfun(@nanmedian,PreIpsiRT.C);
    Ipsi.I{end+1} = ITrial.I - cellfun(@nanmedian,PreIpsiRT.I);
    
    % Find timing of stim relative to start of trials
    Stim = [];
    
    for s = 1:length(stimT)
        % if igor, exp that are not movement triggered have trigs at stim
        % time, do not need to add delay
        
        ind = find(ContT(:,1) < (stimT(s)+50),1,'last');
        if(isnan(CRT))
            Stim(s) = nan;
        else
            Stim(s) = (stimT(s)-ContT(ind,1));
        end
        ind = find(IpsiT(:,2) > stimT(s) ,1);
    end
    
    StimAll{end+1} = Stim;
    
    %     if(length(StimAll) == 6)
    %         keyboard
    %     end
    %
    CondAll{end+1,1} = SL(i).Condition;
    CondAll{end,2} = SL(i).Stim_Delay;
    CondAll{end,3} = SL(i).StimHemi;
    
end

% Cont.C = cellfun(@(x) x(:),Cont.C,'Uniformoutput',false);
% Ipsi.C = cellfun(@(x) x(:),Ipsi.C,'Uniformoutput',false);
% Cont.I = cellfun(@(x) x(:),Cont.I,'Uniformoutput',false);
% Ipsi.I = cellfun(@(x) x(:),Ipsi.I,'Uniformoutput',false);

Cont.C = cellfun(@(x) nanmedian(x,2),Cont.C,'Uniformoutput',false);
Ipsi.C = cellfun(@(x) nanmedian(x,2),Ipsi.C,'Uniformoutput',false);
Cont.I = cellfun(@(x) nanmedian(x,2),Cont.I,'Uniformoutput',false);
Ipsi.I = cellfun(@(x) nanmedian(x,2),Ipsi.I,'Uniformoutput',false);


%% Scatter Plots
for f = 1:2
    if(f==1)
        data = Cont;
        str = 'Contra';
    else
        data = Ipsi;
        str = 'Ipsi';
    end
    figure;
    for i = 1:length(StimAll)
        if(strcmp(SL(1).Animal,'Ubi') && i==14)
            continue;
        end
        
        subplot(2,1,1);
        stimT = nanmedian(StimAll{i});
        if(strcmp(SL(1).Animal,'Ubi'))
            jitterT = 100;
        else
            jitterT = 50;
        end
        d = data.C{i};
        jitter = rand(1,length(d))*2*jitterT-jitterT; %100ms jitter
        scatter(stimT*ones(1,length(d))+jitter,d,'k.'); hold on;
        
        subplot(2,1,2);
        d = data.I{i};
        jitter = rand(1,length(d))*2*jitterT-jitterT; %100ms jitter
        scatter(stimT*ones(1,length(d))+jitter,d,'k.'); hold on;
    end
    subplot(2,1,1);
    title([str,' Trial, Contra Hemi']);
    xlabel('Stim Time (ms)')
    ylabel('\DeltaPMBR');
    
    subplot(2,1,2);
    title([str,' Trial, Ipsi Hemi']);
    xlabel('Stim Time (ms)')
    ylabel('\DeltaPMBR');
end


%% Box Plots
for f = 1:2
    if(f==1)
        data = Cont;
        str = 'Contra';
    else
        data = Ipsi;
        str = 'Ipsi';
    end
    figure;
    Cstim = []; Istim = []; CRT = []; IRT = [];
    for i = 1:length(StimAll)
        if(strcmp(SL(1).Animal,'Ubi') && i==14)
            continue;
        end
        
        stimT = nanmedian(StimAll{i});
        if(strcmp(SL(1).Animal,'Ubi'))
            jitterT = 100;
        else
            jitterT = 50;
        end
        
        condRT = data.C{i};
        jitter = rand(1,length(condRT))*2*jitterT-jitterT; %100ms jitter
        Cstim = [Cstim,stimT+jitter];
        CRT = [CRT;condRT];
        
        stimT = nanmedian(StimAll{i});
        condRT = data.I{i};
        jitter = rand(1,length(condRT))*2*jitterT-jitterT; %100ms jitter
        Istim = [Istim,stimT+jitter];
        IRT = [IRT;condRT];
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
    
    RT = [CRT;IRT]; idx = [Cidx,Iidx];
    bins = 1:((size(bins,1))*2);
    nodata = setdiff(bins,idx);
    if(~isempty(nodata))
        idx = [idx,nodata];
        RT = [RT;zeros(length(nodata),1)];
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
    title([SL(i).Animal,', ',str,' Trials, \DeltaPMBR']); ylabel('PMBR (a.u.)'); xlabel('Stim Time (ms)')
    
    % Plot again
    positions = sort([1:length(unique(idx))/2,(1:length(unique(idx))/2)+1/3]);
    boxplot(RT,idx, 'positions', positions);
    
end


%% Stats
for f = 1:2
    if(f==1)
        data = Cont;
        trial = 'Contra';
    else
        data = Ipsi;
        trial = 'Ipsi';
    end
    Cstim = []; Istim = []; CRT = []; IRT = [];
    for i = 1:length(StimAll)
        if(strcmp(SL(1).Animal,'Ubi') && i==14)
            continue;
        end
        
        stimT = nanmedian(StimAll{i});
        if(strcmp(SL(1).Animal,'Ubi'))
            jitterT = 100;
        else
            jitterT = 50;
        end
        
        condRT = data.C{i};
        jitter = rand(1,length(condRT))*2*jitterT-jitterT; %100ms jitter
        Cstim = [Cstim,stimT+jitter];
        CRT = [CRT;condRT];
        
        stimT = nanmedian(StimAll{i});
        condRT = data.I{i};
        jitter = rand(1,length(condRT))*2*jitterT-jitterT; %100ms jitter
        Istim = [Istim,stimT+jitter];
        IRT = [IRT;condRT];
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
    
    RT = [CRT;IRT]; idx = [Cidx,Iidx];
    bins = 1:((length(bins))*2);
    nodata = setdiff(bins,idx);
    if(~isempty(nodata))
        idx = [idx,nodata];
        RT = [RT,zeros(1,length(nodata))];
    end
    
%     % Remain in "whisker range"
%     w = 1.5;
%     q1 = quantile(RT,0.25);
%     q3 = quantile(RT,0.75);
%     lower = q1-w*(q3-q1);
%     upper = q3+w*(q3-q1);
%     inbound = RT > lower & RT < upper;
%     RT = RT(inbound);
%     idx = idx(inbound);
    
    str = ''; sig = 0.05;
    
    Cidx = unique(Cidx);
    Iidx = unique(Iidx);   
    
    Zstats = length(Cidx)+length(Iidx);
    for i = 1:length(Cidx)
        Zstats(i) = signrank(RT(idx==Cidx(i)));
        stat = Zstats(i);
        sstr = 'ns';
        if(stat < sig)
            sstr = 's';
        end
        str = sprintf('%sContra%d From Zero = %0.3e\t%s\n',str,i,stat,sstr);
    end
    
    str = sprintf('%s\n',str);
    
    for i = 1:length(Iidx)
        Zstats(i) = signrank(RT(idx==Iidx(i)));
        stat = Zstats(i);
        sstr = 'ns';
        if(stat < sig)
            sstr = 's';
        end
        str = sprintf('%sIpsi%d From Zero = %0.3e \t%s\n',str,i,stat,sstr);
    end
    
    str = sprintf('%s\n',str);
    
    Compstats = length(Cidx);
    for i = 1:length(Cidx)
        Compstats(i) = ranksum(RT(idx==Cidx(i)),RT(idx==Iidx(i)));
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
            Cstats(i,j) = ranksum(RT(idx==Cidx(i)),RT(idx==Cidx(j)));
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
            Istats(i,j) = ranksum(RT(idx==Iidx(i)),RT(idx==Iidx(j)));
            stat = Istats(i,j);
            sstr = 'ns';
            if(stat < sig)
                sstr = 's';
            end
            str = sprintf('%sIpsi%d vs Ipsi%d = %0.3e   \t%s\n',str,i,j,stat,sstr);
        end
    end
    str = sprintf('%s\n\n',str);
    
    
    fileID = fopen(['C:\Users\richy.yun\Dropbox\repos\abogaard\efetz\RT manuscript\figures\StatsPMBR',trial,'.txt'],'a');
    fprintf(fileID,'%s Stats\n\n',SL(i).Animal);
    fprintf(fileID,str);
    fclose(fileID);
end