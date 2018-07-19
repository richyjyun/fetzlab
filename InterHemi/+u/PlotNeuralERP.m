function [left,right,trig,ccep] = PlotNeuralERP(SL)
%
% Plots neural data ERPs
% First figure - window of data synced to beginning of left trials
% Second figure - window of data synced to beginning of right trials
% Third figure - window of data in trials with stimulation compared to
% those during the conditioing period without stimulation
%
% RJY March 27 2017

left = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]);
right = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]);
trig = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]);
ccep = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]);


% Loading in the data. They should have the same fs.
if(strcmp(SL.Animal,'Ubi'))
    D = SL.Date;
    S = SL.Session_Guger_Train;
    Session = [char(D),'_',char(S(2))];
    trainfile = [Session,'.f32'];
    gugfile = Session;
    if(~exist(trainfile) || ~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
        return;
    end
    
    StimChns = strsplit(SL.Stim_Loc,'/');
    n = char(StimChns(1));
    StimChns(2) = cellstr([n(1:4),char(StimChns(2))]);
    
    % down sampling rate
    dwn = 10;
    [~, trig1, ~, lefttrials, righttrials, ~, lefttrialsuccess, righttrialsuccess] = u.LoadTrainUbi(trainfile,dwn);
    [data, fs, chnm, ~] = u.LoadGug(gugfile, dwn);
    
elseif(strcmp(SL.Animal,'Kato'))
    D = SL.Date;
    Sessions = strsplit(char(SL.Sessions),'_');
    dwn = 5;
    trig1 = []; lefttrials = []; righttrials = []; lefttrialsuccess = []; righttrialsuccess = []; data = []; last = 0;
    for i = 1:3
        trainfile = [char(D),'_',char(Sessions(i)),'.f32'];
        gugfile = [char(D),'_',char(Sessions(i))];
        if(~exist([gugfile,'.i16']) || ~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
            return;
        end
        if(~exist(trainfile))
            disp(['making f32 file for ' gugfile]);
            cd('F:\Dropbox\repos\abogaard\efetz\DualAccelRt\code');
            err = utils.trainalign2(gugfile, 1);
        end
        
        [acc, t1, lt, rt, ~, lts, rts] = u.LoadTrainKato(trainfile,dwn);
        [d, fs, chnm, ~] = u.LoadGug(gugfile, dwn);     
        
        lefttrials = [lefttrials;last+lt];
        righttrials = [righttrials;last+rt];
        lefttrialsuccess = [lefttrialsuccess;lts];
        righttrialsuccess = [righttrialsuccess;rts];
        data = [data;d];
        if(i == 2)
%             if(strcmp(SL.Condition(end),'M'))
%                 trig1 = u.getKatoTrig(gugfile,fs,1);
%                 trig2 = u.getKatoTrig(gugfile,fs,2);
%             end
            if(any(SL.Condition=='I'))
                trig1 = last+rt(:,1);
            elseif(any(SL.Condition=='C'))
                trig1 = last+lt(:,1);
            end
        end
        last = last+length(acc);
    end
    StimChns = SL.Stim_Loc;
elseif(strcmp(SL.Animal,'Igor'))
    D = SL.Date;
    Sessions = strsplit(char(SL.Sessions),'_');
    data = [];
    for i = 1:3
        gugfile = [char(D),'_',char(Sessions(i))];
        if( ~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
            return;
        end
        [d, fs, chnm] = u.ImportIgorNeuralData(gugfile);
        
        data = [data;d];
        
    end
    StimChns = SL.Stim_Loc; trig1 = SL.trig1; 
    lefttrials = SL.lefttrials; righttrials = SL.righttrials;
    lefttrialsuccess = SL.lefttrialsuccess; righttrialsuccess = SL.righttrialsuccess;
    dwn = 1;
end

% Sometimes there are different number of trials than before due to
% downsampling. Just accounting for that.
if(~strcmp(SL.Animal,'Igor'))
    temp = lefttrials*1000/(fs*dwn);
    check = zeros(length(lefttrials),1);
    for i = 1:length(SL.lefttrials)
        norm = abs(temp(:,1)-SL.lefttrials(i,1));
        ind = find(norm == min(norm));
        check(ind) = 1;
    end
    lefttrials = lefttrials(find(check),:);
    lefttrialsuccess = lefttrialsuccess(find(check));
    
    temp = righttrials*1000/(fs*dwn);
    check = zeros(1,length(righttrials));
    for i = 1:length(SL.righttrials)
        norm = abs(temp(:,1)-SL.righttrials(i,1));
        ind = find(norm == min(norm));
        check(ind) = 1;
    end
    righttrials = righttrials(find(check),:);
    righttrialsuccess = righttrialsuccess(find(check));
end

% Fixing the triggers for this day (Ubi) in case it wasn't fixed earlier
if(strcmp(SL.Date,'20170127'))
    trig1(1:11) = [];
end

% % For plotting around triggers
% trials = [];
% if strcmp(SL.Condition(1:4),'Ipsi')
%     trials = lefttrials;
% elseif strcmp(SL.Condition(1:6),'Contra')
%     trials = righttrials;
% end


% Band pass filter to extract TRPs
% bpf = [1 10]; % chL/R transform: 1. bandpass filter (Hz) (around beta wave range)
% [bbpf,abpf] = butter(1,bpf/(fs/2)); % 1st order bandpass (see daq_sapi_*)
% Filter = filtfilt(bbpf,abpf,double(data));

Filter = data;
clear data;

IpsiNeural = []; ContNeural = [];
bounds = [trig1(1)*SL.fs/1000/dwn,trig1(end)*SL.fs/1000/dwn];
ctr = 1;
Channels = [];
for i = 1:size(Filter,2)
    % if doesn't have L or R, or matches a stimulation channel, skip
    if(~(any(char(chnm(i)) == 'L') || any(char(chnm(i)) == 'R') || any(char(chnm(i)) == 'E') || any(char(chnm(i)) == 'I')) ||...
            strcmp(chnm(i),char(StimChns(1))) || strcmp(chnm(i),char(StimChns(2))))
        continue;
    end
    
    Channels(end+1) = i;
    
    window = 1.2;
    window = floor(window*fs);
    back = 1/3;
    
    % plots based on Ipsi trials
    figure(left); subplot(6,5,ctr);
    set(gca, 'fontsize', 7)
    
    if(strcmp(SL.StimHemi,'R'))
        ipsitrials = righttrials;
        ipsisuccess = righttrialsuccess;
        ipsirt = SL.rts_r;
        conttrials = lefttrials;
        contsuccess = lefttrialsuccess;
        contrt = SL.rts_l;
    elseif(strcmp(SL.StimHemi,'L'))
        ipsitrials = lefttrials;
        ipsisuccess = lefttrialsuccess;
        ipsirt = SL.rts_l;
        conttrials = righttrials;
        contsuccess = righttrialsuccess;
        contrt = SL.rts_r;
    else
        continue;
    end
    
    % save left or right cond to plot against trigger
    Nostim = [];
    
    pre = [];
    cond = [];
    post = [];
    for j = 1:length(ipsitrials)
        if(isnan(ipsirt(j)))
            continue;
        end
        if(ipsisuccess(j))
            start = floor(ipsitrials(j,1)*SL.fs/1000/dwn); fin = floor(ipsitrials(j,2)*SL.fs/1000/dwn);
            if((start+floor(ipsirt(j))-window)<0 || (start+floor(ipsirt(j))+window)>length(Filter))
                continue;
            end
            if(fin < bounds(1))
                pre(end+1,:) = Filter((start+floor(ipsirt(j))-window):(start+floor(ipsirt(j))+window),i);
            elseif(start>(bounds(2)+50))
                post(end+1,:) = Filter((start+floor(ipsirt(j))-window):(start+floor(ipsirt(j))+window),i);
            elseif(~any((trig1./dwn)>start-50 & (trig1./dwn)<fin))
                cond(end+1,:) = Filter((start+floor(ipsirt(j))-window):(start+floor(ipsirt(j))+window),i);
            end
        end
    end
    
    t = (-window:window).*1000./fs;
    plot(t,median((pre-median(pre,2))), 'color', 'k', 'linewidth', 1.5), hold on;
    if ~isempty(cond)
        plot(t,median((cond-median(cond,2))), 'color', 'b', 'linewidth', 1.5, 'linestyle', ':'), hold on;
    end
    plot(t,median((post-median(post,2))), 'color', 'r', 'linewidth', 1.5, 'linestyle', '--'), hold off;
    xlim([t(1),t(end)]); title(char(chnm(i)),'fontsize', 9);
    yl = ylim;
    if(strcmp(SL.Animal,'Ubi') && ~isempty(cond))
        IpsiNeural(end+1,:) = [sum(abs(median((pre-median(pre,2))))),sum(abs(median((cond-median(cond,2))))),sum(abs(median((post-median(post,2)))))];
    end
%     if(length(median(pre))==length(median(cond)) && length(median(cond)) == length(median(post)))
%         cor(1) = corr(median(pre)',median(cond)');
%         cor(2) = corr(median(pre)',median(post)');
%         cor(3) = corr(median(cond)',median(post)');
%         text(t(1)+20,yl(2)/2,{num2str(cor(1)),num2str(cor(2)),num2str(cor(3))},'fontsize',7);
%     end
    
    if strcmp(SL.Condition(1:4),'Ipsi')
        Nostim = cond;
    end
    
    % plots based on right trials
    figure(right); subplot(6,5,ctr);
    set(gca, 'fontsize', 7)
    
    pre = [];
    cond = [];
    post = [];
    for j = 1:length(conttrials)
        if(isnan(contrt(j)))
            continue;
        end
        if(contsuccess(j))
            start = floor(conttrials(j,1)*SL.fs/1000/dwn); fin = floor(conttrials(j,2)*SL.fs/1000/dwn);
            if((start+floor(contrt(j))-window)<0 || (start+floor(contrt(j))+window)>length(Filter))
                continue;
            end
            if(fin < bounds(1))
                pre(end+1,:) = Filter((start+floor(contrt(j))-window):(start+floor(contrt(j))+window),i);
            elseif(start>(bounds(2)+50))
                post(end+1,:) = Filter((start+floor(contrt(j))-window):(start+floor(contrt(j))+window),i);
            elseif(~any((trig1./dwn)>start-50 & (trig1./dwn)<fin))
                cond(end+1,:) = Filter((start+floor(contrt(j))-window):(start+floor(contrt(j))+window),i);
            end
            %             if((start-window*back)<0 || (start+window)>length(Filter))
            %                 continue;
            %             end
            %             if(fin < bounds(1))
            %                 pre(end+1,:) = Filter((start-floor(window*back)):(start+window),i);
            %             elseif(start>(bounds(2)+50))
            %                 post(end+1,:) = Filter((start-floor(window*back)):(start+window),i);
            %             elseif(~any((trig1./dwn)>start-50 & (trig1./dwn)<fin))
            %                 cond(end+1,:) = Filter((start-floor(window*back)):(start+window),i);
            %             end
        end
    end
    
    plot(t,median((pre-median(pre,2))), 'color', 'k', 'linewidth', 1.5), hold on;
    if ~isempty(cond)
        plot(t,median((cond-median(cond,2))), 'color', 'b', 'linewidth', 1.5, 'linestyle', ':'), hold on;
    end
    plot(t,median((post-median(post,2))), 'color', 'r', 'linewidth', 1.5, 'linestyle', '--'), hold off;
    xlim([t(1),t(end)]); title(char(chnm(i)),'fontsize', 9);
    yl = ylim;
    if(strcmp(SL.Animal,'Ubi') && ~isempty(cond))
        ContNeural(end+1,:) = [sum(abs(median((pre-median(pre,2))))),sum(abs(median((cond-median(cond,2))))),sum(abs(median((post-median(post,2)))))];
    end
%     if(length(mean(pre))==length(mean(cond)) && length(mean(cond)) == length(mean(post)))
%         cor(1) = corr(mean(pre)',mean(cond)');
%         cor(2) = corr(mean(pre)',mean(post)');
%         cor(3) = corr(mean(cond)',mean(post)');
%         text(t(1)+20,yl(2)/2,{num2str(cor(1)),num2str(cor(2)),num2str(cor(3))},'fontsize',7);
%     end
    
    %     window = 1.2;
    %     window = floor(window*fs);
    %     back = 1/3;
    
    if strcmp(SL.Condition(1:6),'Contra')
        Nostim = cond;
    end
    
    % plot around triggers (if contra, plot against right trials without
    % stim during conditioning period, and same with ipsi)
    figure(trig); subplot(6,5,ctr);
    set(gca, 'fontsize', 7)
    
%     if strcmp(SL.Condition(1:4),'Ipsi')
%         trials = ipsitrials;
%     elseif strcmp(SL.Condition(1:6),'Contra')
%         trials = conttrials;
%     end
    
    snips = zeros(length(trig1),window*2+1);
    for j = 1:length(trig1)
        %         norm = abs(trials(:,1)-trig1(j));
        %         ind = find(norm == min(norm));
        %         start = floor(trials(ind,1)/dwn); fin = floor(trials(ind,2)/dwn);
        %         if((fin-window)<0 || (fin+window*back)>length(Filter))
        %             continue;
        %         end
        %         snips(end+1,:) = Filter(fin-window:fin+floor(window*back),i);
        if(floor(trig1(j)*SL.fs/1000/dwn)-window<0 || floor(trig1(j)*SL.fs/1000/dwn)+window > length(Filter))
            continue;
        end
        snips(j,:) = Filter(floor(trig1(j)*SL.fs/1000/dwn)-window:floor(trig1(j)*SL.fs/1000/dwn)+window,i)';
    end
    
    t2 = (-window:window).*1000./fs;
    plot(t2,median(snips-median(snips,2)), 'color', 'k', 'linewidth', 1.5)
    %     if(~isempty(Nostim) && strcmp(SL.Animal,'Ubi'))
    %         hold on; plot(t,mean(Nostim)-mean(mean(Nostim)),'color','r', 'linewidth', 1.5, 'linestyle', ':'); hold off;
    %     end
    xlim([t2(1),t2(end)]); title(char(chnm(i)),'fontsize', 9);
    stimHemi = SL.StimHemi;
    yl = ylim;
    if(strcmp(SL.Animal,'Ubi'))
        if(any(char(chnm(i)) == stimHemi))
            if ((yl(2)-yl(1))>200)
                ylim([-100,100]);
            end
        end
    end
    if(strcmp(SL.Animal,'Kato') || strcmp(SL.Animal,'Igor'))
        if ((yl(2)-yl(1))>200)
            ylim([-100,100]);
        end
    end
    
    % plot CCEPs
    figure(ccep); subplot(6,5,ctr);
    set(gca, 'fontsize', 7)
    
    snips = [];
    snips1 = [];
    for j = 1:floor(length(trig1)/2)
        %         norm = abs(trials(:,1)-trig1(j));
        %         ind = find(norm == min(norm));
        %         start = floor(trials(ind,1)/dwn); fin = floor(trials(ind,2)/dwn);
        if(floor(trig1(j)*SL.fs/1000/dwn)-window<0 || floor(trig1(j)*SL.fs/1000/dwn)+window>length(Filter))
            continue;
        end
        snips(end+1,:) = Filter(floor(trig1(j)*SL.fs/1000/dwn)-window:floor(trig1(j)*SL.fs/1000/dwn)+window,i);
    end
    for j = ceil(length(trig1)/2):length(trig1)
        %         norm = abs(trials(:,1)-trig1(j));
        %         ind = find(norm == min(norm));
        %         start = floor(trials(ind,1)/dwn); fin = floor(trials(ind,2)/dwn);
        if(floor(trig1(j)*SL.fs/1000/dwn)-window<0 || floor(trig1(j)*SL.fs/1000/dwn)+window>length(Filter))
            continue;
        end
        snips1(end+1,:) = Filter(floor(trig1(j)*SL.fs/1000/dwn)-window:floor(trig1(j)*SL.fs/1000/dwn)+window,i);
    end
    
    %     if(strcmp(SL.Animal,'Ubi'))
    plot(t2,median(snips-median(snips,2)), 'color', 'k', 'linewidth', 1.5); hold on;
    plot(t2,median(snips1-median(snips1,2)), 'color', 'r', 'linewidth', 1.5, 'linestyle', '--'); hold off;
    %     end
    xlim([t2(1),t2(end)]); title(char(chnm(i)),'fontsize', 9);
    yl = ylim;
    stimHemi = SL.StimHemi;
    if(strcmp(SL.Animal,'Ubi'))
        if(any(char(chnm(i)) == stimHemi))
            if ((yl(2)-yl(1))>200)
                ylim([-100,100]);
            end
        end
    end
    if(strcmp(SL.Animal,'Kato') || strcmp(SL.Animal,'Igor'))
        if ((yl(2)-yl(1))>200)
            ylim([-100,100]);
        end
    end
    

    %     if(strcmp(SL.Animal,'Ubi'))
    %         % subtract trials without stimulation to see clear evoked response
    %         % without TRPs involved
    %         half = floor(sum((trials(:,1)>trig1(1) & trials(:,1)<trig1(floor(length(trig1)/2))+50))+1-floor(length(trig1)/2));
    %         if(~exist('half') || isnan(half) || half == 0 || half == length(Nostim))
    %             continue;
    %         end
    %         if (~isempty(Nostim) && size(Nostim,1)>half)
    %             first = (mean(snips)-mean(mean(snips)))-(mean(Nostim(1:half,:))-mean(mean(Nostim(1:half,:))));
    %             second = (mean(snips1)-mean(mean(snips1)))-(mean(Nostim(half+1:end,:))-mean(mean(Nostim(half+1:end,:))));
    %             plot(t2,first, 'color', 'k', 'linewidth', 1.5); hold on;
    %             plot(t2,second, 'color', 'r', 'linewidth', 1.5, 'linestyle', '--'); hold off;
    %         else
    %             plot(t2,mean(snips)-mean(mean(snips)), 'color', 'k', 'linewidth', 1.5); hold on;
    %             plot(t2,mean(snips1)-mean(mean(snips1)), 'color', 'r', 'linewidth', 1.5, 'linestyle', '--'); hold off;
    %         end
    %         xlim([t2(1),t2(end)]); title(char(chnm(i)),'fontsize', 9);
    %         if(any(char(chnm(i)) == stimHemi))
    %             yl = ylim;
    %             if ((yl(2)-yl(1))>100)
    %                 ylim([-50,50]);
    %             end
    %         end
    %     end
    
    
    ctr = ctr+1;
end

figure(left); a = axes; t1 = title(['Ipsi Trial ERPs; ', char(SL.Condition),'; StimHemi ', SL.StimHemi]);
a.Visible = 'off'; t1.Visible = 'on';
figure(right); a = axes; t1 = title(['Contra Trial ERPs; ', char(SL.Condition),'; StimHemi ', SL.StimHemi]);
a.Visible = 'off'; t1.Visible = 'on';
figure(trig); a = axes; t1 = title(['Stim Cond Trials; ', char(SL.Condition),'; StimHemi ', SL.StimHemi]);
a.Visible = 'off'; t1.Visible = 'on';
figure(ccep); a = axes; t1 = title(['Early vs Late Stim Cond Trials; ', char(SL.Condition),'; StimHemi ', SL.StimHemi]);
a.Visible = 'off'; t1.Visible = 'on';

if(strcmp(SL.Animal,'Ubi'))
    SL.IpsiNeural = IpsiNeural;
    SL.ContNeural = ContNeural;
    SL.Channels = chnm(Channels);
end

end
