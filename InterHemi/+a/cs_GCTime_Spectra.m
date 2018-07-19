close all; clear; pack

%% Ubi
fname = fullfile('F:\Dropbox\repos\abogaard\efetz\U\packets', 'ubi-Coh_GCTime_Windowed.ps');
if(exist(fname))
    delete(fname);
end
load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat')

%% Igor
% fname = fullfile('F:\Dropbox\repos\abogaard\efetz\U\packets', 'igor-Coh_GCTime_Windowed.ps');
% if(exist(fname))
%     delete(fname);
% end
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgor.mat')

SessionList = SL;

GC_all = struct([]); % Contains all GC data

for day = 1:length(SessionList)
    SL = SessionList(day);
    disp(['Session ', SL.Date,': ',num2str(round(day/length(SessionList)*100,1)),'%']);
    if(strcmp(SL.Bad,'1') || strcmp(SL.Condition,'nostim') || strcmp(SL.Condition,'tonic') || strcmp(SL.Condition(end),'R'))
        continue;
    end
    
    if(str2num(SL.Date)>=20170226 && str2num(SL.Date)<=20170306) % days with 100% stim, ignore for the time being
        continue;
    end
    
    %% Load Data
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
        [accel, trig1, ~, lefttrials, righttrials, ~, lefttrialsuccess, righttrialsuccess] = u.LoadTrainUbi(trainfile,dwn);
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
        SL.lefttrials = SL.lefttrials*fs/1000; SL.righttrials = SL.righttrials*fs/1000;
        SL.trig1 = SL.trig1*fs/1000; trig1 = SL.trig1;
        lefttrialsuccess = SL.lefttrialsuccess; righttrialsuccess = SL.righttrialsuccess;
        dwn = 1;
    end
    
    if(str2num(SL.Date) == 20170127)
        trig1(1:11) = [];
    end
    
    GC_all(end+1).Session = SL.Date;
    
    % fo = 60;  q = 35; bw = (fo/(fs/2))/q; % comb notch filter for 60Hz noise
    % [BB,A] = iircomb(fs/fo,bw,'notch'); % Note type flag 'notch'
    Filter = data; clear data;
    
    % Re-calculate reaction time using upsampled, aligned data. Everything
    % is already aligned with Igor
    if(~strcmp(SL.Animal,'Igor'))
        SL.accel_raw_r = accel(:,2);
        SL.accel_raw_l = accel(:,1);
        SL.lefttrials = lefttrials;
        SL.righttrials = righttrials;
        SL.lefttrialsuccess = lefttrialsuccess;
        SL.righttrialsuccess = righttrialsuccess;
        SL.trig1 = trig1;
        SL.fs = fs;
          
        SL.accelfs = fs; %needed for calculating reaction time correctly. Also for plotexperiment later.
        
        SL = a.AppendReactionTimes(SL);
        SL = a.AppendNormalizedDelay(SL);
    end
    GC_all(end).Delay = SL.NormDelay;
    GC_all(end).Condition = SL.Condition;
    
    %% Define triggers and window
    % lefttrialRT = lefttrials(:,1) + SL.rts_l./1000.*(fs*dwn);
    % lefttrialsuccess(isnan(SL.rts_l)) = 0;
    ofs = 250; %RP starts ~250ms before reaction time calculation (at least for Ubi)
    
    if(strcmp(SL.Condition,'Control')) % if it's a control session, put fake triggers in the middle 1/2 of experiment
        trials = [lefttrials(:,1);righttrials(:,1)];
        [~,ind] = sort(trials);
        trig1(1) = trials(ind(floor(length(ind)/4)));
        trig1(2) = trials(ind(floor(3*length(ind)/4)));
        SL.trig1 = trig1;
    end
    
    lefttrials = SL.lefttrials+(SL.rts_l-ofs)/1000*SL.fs; lefttrials = lefttrials(~isnan(SL.rts_l)& lefttrialsuccess,:);
    trigPreL = lefttrials((lefttrials(:,1)<trig1(1)));
    trigPostL = lefttrials((lefttrials(:,1)>trig1(end)));
    trigCondL = lefttrials(lefttrials(:,1)>=trig1(1) & lefttrials(:,1)<=trig1(end),:);
    GC_all(end).RTL(1) = nanmedian(SL.rts_l(SL.lefttrials(:,1)<trig1(1)));
    GC_all(end).RTL(2) = nanmedian(SL.rts_l(SL.lefttrials(:,1)>=trig1(1) & SL.lefttrials(:,1)<=trig1(end)));
    GC_all(end).RTL(3) = nanmedian(SL.rts_l(SL.lefttrials(:,1)>trig1(end)));
    remove = [];
    if(strcmp(SL.Condition(1),'I') && strcmp(SL.Animal,'Ubi'))
        for i = 1:length(trig1)
            norm = abs(trigCondL(:,1)-trig1(i));
            ind = find(norm == min(norm));
            remove(end+1) = ind;
        end
    end
    trigCondL(remove,:) = [];
    trigCondL = trigCondL(:,1);%+SL.rts_l(length(trigPreL)+1:length(trigPreL)+length(trigCondL));
    shift =0;% nanmedian(SL.rts_l);
    
    righttrials = SL.righttrials+(SL.rts_r-ofs)/1000*SL.fs; righttrials = righttrials(~isnan(SL.rts_r) & righttrialsuccess,:);
    trigPreR = righttrials(righttrials(:,1)<trig1(1),1);
    trigPostR = righttrials(righttrials(:,1)>trig1(end),1);
    trigCondR = righttrials(righttrials(:,1)>=trig1(1) & righttrials(:,1)<=trig1(end),:);
    GC_all(end).RTR(1) = nanmedian(SL.rts_r(SL.righttrials(:,1)<trig1(1)));
    GC_all(end).RTR(2) = nanmedian(SL.rts_r(SL.righttrials(:,1)>=trig1(1) & SL.righttrials(:,1)<=trig1(end)));
    GC_all(end).RTR(3) = nanmedian(SL.rts_r(SL.righttrials(:,1)>trig1(end)));
    remove = [];
    if(strcmp(SL.Condition(1),'C')&& strcmp(SL.Animal,'Ubi'))
        for i = 1:length(trig1)
            norm = abs(trigCondR(:,1)-trig1(i));
            ind = find(norm == min(norm));
            remove(end+1) = ind;
        end
    end
    trigCondR(remove,:) = [];
    trigCondR = trigCondR(:,1);
    
    % window - time in s
    window = 0.5;
    
    %% Plot Experiment
    h = u.PlotExperiment(SL);
    print(h, '-dpsc2', fname, '-append'); close(h);
    
    %% Do GC
    tic
    % Define parameters
    movingwin = [0.25,0.01];
    win = round(movingwin(1)*fs); wstep = round(movingwin(2)*fs);
    
    params.tapers = [3,5]; params.Fs = fs; params.fpass = [9,100]; params.trialave = 1;
    
    % specify channels to look at, for control as well
    if(strcmp(SL.Animal,'Ubi'))
        if(length(SL.Stim_Loc) < 5)
            lchn = 24;
        else
            stim = str2num(SL.Stim_Loc(5));
            if(stim == 3), lchn = 24; else, lchn = 22; end;
        end
        
        chnsR = [5,6,11,12,15,16];
        chnsL = [lchn,lchn+1,[9,10,13,14,15,16]+17];
    elseif(strcmp(SL.Animal,'Igor'))
        chnsR = [17,18];
        chnsL = 3:8;
    end
    
    % parameters for GC
    params.regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default) tried LWR, takes much longer but gives same results
    params.icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
    
    params.morder    = round(fs*0.05); % 50ms, from PNAS paper
    params.momax     = params.morder*3;     % maximum model order for model order estimation. Default takes way too long, gives same results. Apparently often the same as model order, but 3* to give it more
    
    params.acmaxlags = [];   % maximum autocovariance lags (empty for automatic calculation)
    
    params.statt     = 'F';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
    params.alpha_sig     = 0.005;   % significance level for significance test
    params.mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
    
    params.fs        = fs;    % sample rate (Hz)
    params.fres      = [];     % frequency resolution (empty for automatic calculation)
    
    params.nvars     = 2;     % number of variables. only doing pairwise
    
    params.frange    = [9,100]; %range of frequencies to look at
        
    % do GC for all trial types
    for trial = 1:2
        Trial_all = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(Trial_all,'visible','off');
        Trial_avgl2r = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(Trial_avgl2r,'visible','off');
        Trial_avgr2l = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(Trial_avgr2l,'visible','off');
        yl1 = [0,0]; yl2 = [0,0]; yl3 = [0,0];
        
        for epoch = 1:3
            % Pick trials and window, define indices
            if trial == 1 % left trials
                if epoch == 1 % pre
                    trig = trigPreL;
                    epochtype = 'Pre';
                elseif epoch == 2
                    trig = trigCondL;
                    epochtype= 'Cond';
                else
                    trig = trigPostL;
                    epochtype = 'Post';
                end
                trialtype = 'Left';
            else
                if epoch == 1 % pre
                    trig = trigPreR;
                    epochtype = 'Pre';
                elseif epoch == 2
                    trig = trigCondR;
                    epochtype= 'Cond';
                else
                    trig = trigPostR;
                    epochtype = 'Post';
                end
                trialtype = 'Right';
            end   
            
            if(~strcmp(SL.Animal,'Ubi') && epoch == 2)
                if((strcmp(SL.Condition(1),'I') && strcmp(SL.StimHemi,'L')) ||...
                        (strcmp(SL.Condition(1),'C') && strcmp(SL.StimHemi,'R')))
                    hand = 1; % left hand trials have stim
                else
                    hand = 2; % right hand trials have stim
                end
                if(trial == hand)
                    continue;
                end
            end
            
            inds = floor(-0.5*fs:1:1.5*fs); %inds = floor(-window*fs:1:2*window*fs);
            
            trialinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
            trialinds(:,floor(trialinds(1,:))<=0) = []; trialinds(:,floor(trialinds(end,:))>length(Filter)) = [];
            
            params.ntrials   = length(trig);     % number of trials
            params.nobs      = length(inds);   % number of observations per trial
            
            % Set variable in correct format for GC
            X = [];
            gc_time = [];
            
            for i = 1:length(chnsR)
                d = Filter(:,chnsR(i));
                X(end+1,:,:) = u.meanSubtract(d(floor(trialinds)),params);
            end
            for i = 1:length(chnsL)
                d = Filter(:,chnsL(i));
                X(end+1,:,:) = u.meanSubtract(d(floor(trialinds)),params);
            end
            
            % Do GC in windows
            winstart = 1:wstep:size(X,2)-win+1; winmid=winstart+round(win/2);
            
            startadd = nan;
            for k = 1:length(winstart)
                disp(['Session ', SL.Date,': ',num2str(round(day/length(SessionList)*100,1)),'%']);
                disp([trialtype,', ',epochtype])
                disp([num2str(k),'/',num2str(length(winstart))]);
                if(isnan(startadd))
                    Lbound = winstart(k); Rbound = Lbound+win-1;
                else
                    Lbound = winstart(startadd); Rbound = winstart(k)+win-1;
                end
                x = X(:,Lbound:Rbound,:);
                [gc,p,s] = a.GC_time(x,params);
                gc_time(k,:,:) = gc;
                % if causality diverges, add more to the time span
                if(isnan(gc))
                    startadd = k;
                else
                    startadd = nan;
                end
            end
            
            GC_all(end).([trialtype,epochtype]) = gc_time;
            
            % Readiness potential to plot with GC
            RP1 = squeeze(X(length(chnsR),:,:)); RP1 = mean(RP1,2);
            RP2 = squeeze(X(length(chnsR)+1,:,:)); RP2 = mean(RP2,2);
            
            % Plot the changes in causality over time for all channels
            set(0, 'CurrentFigure', Trial_all); 
            for i = 1:length(chnsR)
                figure;
                for j = 1:length(chnsL)
%                     subplot(3,2,(epoch-1)*2+1); hold on;
                    subplot(2,1,1); hold on;
                    x = winmid(1:size(gc_time,1))/fs; y = gc_time(:,i,length(chnsR)+j);
                    plot(x(~isnan(y)),y(~isnan(y))); title([epochtype,' L2R']); xlim([x(1),x(end)]);
                    yl = ylim;
                    if(yl(1)<yl1(1)|| yl1(1)==0)
                        yl1(1) = yl(1);
                    end
                    if(yl(2)>yl1(2))
                        yl1(2) = yl(2);
                    end
%                     subplot(3,2,epoch*2); hold on;
                    subplot(2,1,2); hold on;
                    y = gc_time(:,length(chnsR)+j,i);
                    plot(x(~isnan(y)),y(~isnan(y))); title([epochtype,' R2L']); xlim([x(1),x(end)]);
                    yl = ylim;
                    if(yl(1)<yl1(1))
                        yl1(1) = yl(1);
                    end
                    if(yl(2)>yl1(2))
                        yl1(2) = yl(2);
                    end
                end
            end

            
            % Plot the averages of above
            set(0, 'CurrentFigure', Trial_avgl2r); subplot(3,1,epoch);
            x = winmid(1:size(gc_time,1))/fs; y = squeeze(nanmean(nanmean(gc_time(:,1:length(chnsR),length(chnsR)+1:end),2),3));
            yyaxis left; plot(x(~isnan(y)),y(~isnan(y)));
            yl = ylim;
            if(yl(1)<yl2(1)|| yl2(1)==0)
                yl2(1) = yl(1);
            end
            if(yl(2)>yl2(2))
                yl2(2) = yl(2);
            end
            yyaxis right; plot((1:length(RP1))/fs,RP1); hold on; plot((1:length(RP2))/fs,RP2);  xlim([x(1),x(end)]);
            title([trialtype,' ',epochtype,' L2R']);
            
            set(0, 'CurrentFigure', Trial_avgr2l); subplot(3,1,epoch);
            y = squeeze(nanmean(nanmean(gc_time(:,length(chnsR)+1:end,1:length(chnsR)),2),3));
            yyaxis left; plot(x(~isnan(y)),y(~isnan(y)));
            yl = ylim;
            if(yl(1)<yl3(1) || yl3(1)==0)
                yl3(1) = yl(1);
            end
            if(yl(2)>yl3(2))
                yl3(2) = yl(2);
            end
            yyaxis right; plot((1:length(RP1))/fs,RP1); hold on; plot((1:length(RP2))/fs,RP2); xlim([x(1),x(end)]);
            title([trialtype,' ',epochtype,' R2L']);
            
        end
        
        % set ylim to be equal within a page, and print figures out to file
        set(0, 'CurrentFigure', Trial_all); 
        for i = 1:6
            subplot(3,2,i); ylim(yl1);
        end
        ax = axes; t1 = title([trialtype,' Trials_']); ax.Visible = 'off'; t1.Visible = 'on';
        print(Trial_all, '-dpsc2', fname, '-append'); close(Trial_all);
        
        for i = 1:3
            set(0, 'CurrentFigure', Trial_avgl2r); 
            subplot(3,1,i); yyaxis left; ylim(yl2);
            set(0, 'CurrentFigure', Trial_avgr2l); 
            subplot(3,1,i); yyaxis left; ylim(yl3);
        end
        
        print(Trial_avgl2r, '-dpsc2', fname, '-append'); close(Trial_avgl2r);
        print(Trial_avgr2l, '-dpsc2', fname, '-append'); close(Trial_avgr2l);
        
    end
    
    toc
    
%     save('F:\Dropbox\repos\abogaard\efetz\U\Mat_Data\GC_timeRPost.mat','gc_time','-v7.3');

    
%     %% Plot gc and p values (GC_time plot not working for some reason)
%     Trial = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(Trial,'visible','off');
%     clim = [0,0];
%     for i = 1:3
%         for j = 1:3
%             epoch = '';
%             if j ==1
%                 epoch = 'Pre';
%             elseif j == 2
%                 epoch = 'Cond';
%             else
%                 epoch = 'Post';
%             end
%             
%             subplot(6,3,(i-1)*3+j);
%             d = LGC((i-1)*3+j).gc;
%             d(1:length(chnsR),1:length(chnsR)) = 0;    d(length(chnsR)+1:end,length(chnsR)+1:end) = 0;
%             plot_pw(d); title(['LTrial, Window ',num2str(i),', ',epoch],'fontsize',7); colorbar;
%             set(gca,'XTickLabel',chnm([chnsR,chnsL])); xtickangle(45);
%             set(gca,'YTickLabel',chnm([chnsR,chnsL]));
%             if(diff(caxis)>diff(clim))
%                 clim = caxis;
%             end
%             
%             subplot(6,3,9+(i-1)*3+j);
%             d = RGC((i-1)*3+j).gc;
%             d(1:length(chnsR),1:length(chnsR)) = 0;    d(length(chnsR)+1:end,length(chnsR)+1:end) = 0;
%             plot_pw(d); title(['RTrial, Window ',num2str(i),', ',epoch],'fontsize',7); colorbar;
%             set(gca,'XTickLabel',chnm([chnsR,chnsL])); xtickangle(45);
%             set(gca,'YTickLabel',chnm([chnsR,chnsL]));
%             if(diff(caxis)>diff(clim))
%                 clim = caxis;
%             end
%         end
%     end
%     for i = 1:18
%         figure(Trial); subplot(6,3,i); caxis(clim);
%     end
%     print(Trial, '-dpsc2', fname, '-append'); close(Trial);
%     
%     %% Plot differences between cond/post and pre
%     LTrial = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(LTrial,'visible','off');
%     clim = [0,0];
%     % Create colormap that is red for negative, green for positive,
%     % and a chunk inthe middle that is black.
%     greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
%     redColorMap = [linspace(1, 0, 124), zeros(1, 132)];
%     colorMap = [redColorMap; greenColorMap; zeros(1, 256)]';
%     % Apply the colormap.
%     colormap(colorMap);
%     for i = 1:3
%         subplot(3,2,(i-1)*2+1);
%         d = LGC((i-1)*3+2).gc-LGC((i-1)*3+1).gc;    % Difference between Cond and Pre
%         d(1:length(chnsR),1:length(chnsR)) = 0;     % Removing R2R GC
%         d(length(chnsR)+1:end,length(chnsR)+1:end) = 0;   % Removing L2L GC
%         d(LGC((i-1)*3+2).s == 0) = 0; d(LGC((i-1)*3+1).s == 0) = 0;    % Removing statistically insignificant GC
%         imagesc(d); n = size(d,1);      % Plot in grid
%         axis('square'); xlabel('from'); ylabel('to');
%         set(gca,'XTick',1:n); set(gca,'XTickLabel',chnm([chnsR,chnsL])); xtickangle(45);
%         set(gca,'YTick',1:n); set(gca,'YTickLabel',chnm([chnsR,chnsL])); title(['Window ',num2str(i),' Cond-Pre'])
%         colorbar;
%         if(diff(caxis)>diff(clim))  % Normalizing colorbar across all
%             clim = caxis;
%         end
%         l2r = d(1:length(chnsR),length(chnsR)+1:end); r2l = d(length(chnsR)+1:end,1:length(chnsR));
%         GC_all(end).LDiffCond(i,:) = [mean2(l2r),mean2(r2l)];   % Store the mean values
%         text(length(chnsR)+ceil(length(chnsL)/2),length(chnsR)+1,num2str(mean2(l2r),'%10.2e'),'color','white');
%         text(length(chnsR)+1,length(chnsR)+ceil(length(chnsL)/2),num2str(mean2(r2l),'%10.2e'),'color','white');
%         
%         subplot(3,2,(i-1)*2+2);
%         d = LGC((i-1)*3+3).gc-LGC((i-1)*3+1).gc;    % Same as above but Post - Pre
%         d(LGC((i-1)*3+3).s == 0) = 0; d(LGC((i-1)*3+1).s == 0) = 0;
%         d(1:length(chnsR),1:length(chnsR)) = 0;
%         d(length(chnsR)+1:end,length(chnsR)+1:end) = 0;
%         imagesc(d); n = size(d,1);
%         axis('square'); xlabel('from'); ylabel('to');
%         set(gca,'XTick',1:n); set(gca,'XTickLabel',chnm([chnsR,chnsL])); xtickangle(45);
%         set(gca,'YTick',1:n); set(gca,'YTickLabel',chnm([chnsR,chnsL])); title(['Window ',num2str(i),' Post-Pre'])
%         colorbar;
%         if(diff(caxis)>diff(clim))
%             clim = caxis;
%         end
%         l2r = d(1:length(chnsR),length(chnsR)+1:end); r2l = d(length(chnsR)+1:end,1:length(chnsR));
%         GC_all(end).LDiffPost(i,:) = [mean2(l2r),mean2(r2l)];
%         text(length(chnsR)+ceil(length(chnsL)/2),length(chnsR)+1,num2str(mean2(l2r),'%10.2e'),'color','white');
%         text(length(chnsR)+1,length(chnsR)+ceil(length(chnsL)/2),num2str(mean2(r2l),'%10.2e'),'color','white');
%         
%     end
%     clim = [-max(abs(clim)),max(abs(clim))];
%     for i = 1:6
%         subplot(3,2,i)
%         caxis(clim);
%     end
%     ax = axes; t1 = title('Left Trials'); ax.Visible = 'off'; t1.Visible = 'on'; 
%     print(LTrial, '-dpsc2', fname, '-append'); close(LTrial);
%     
%     % Same as above but for right trials
%     % Could make into function, but eh
%     RTrial = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(RTrial,'visible','off');
%     clim = [0,0]; colormap(colorMap);
%     for i = 1:3
%         subplot(3,2,(i-1)*2+1);
%         d = RGC((i-1)*3+2).gc-RGC((i-1)*3+1).gc;
%         d(1:length(chnsR),1:length(chnsR)) = 0;
%         d(length(chnsR)+1:end,length(chnsR)+1:end) = 0;
%         d(LGC((i-1)*3+2).s == 0) = 0; d(LGC((i-1)*3+1).s == 0) = 0;
%         imagesc(d); n = size(d,1);
%         axis('square'); xlabel('from'); ylabel('to');
%         set(gca,'XTick',1:n); set(gca,'XTickLabel',chnm([chnsR,chnsL])); xtickangle(45);
%         set(gca,'YTick',1:n); set(gca,'YTickLabel',chnm([chnsR,chnsL])); title(['Window ',num2str(i),' Cond-Pre'])
%         colorbar;
%         if(diff(caxis)>diff(clim))
%             clim = caxis;
%         end
%         l2r = d(1:length(chnsR),length(chnsR)+1:end); r2l = d(length(chnsR)+1:end,1:length(chnsR));
%         GC_all(end).RDiffCond(i,:) = [mean2(l2r),mean2(r2l)];
%         text(length(chnsR)+ceil(length(chnsL)/2),length(chnsR)+1,num2str(mean2(l2r),'%10.2e'),'color','white');
%         text(length(chnsR)+1,length(chnsR)+ceil(length(chnsL)/2),num2str(mean2(r2l),'%10.2e'),'color','white');
%         
%         subplot(3,2,(i-1)*2+2);
%         d = RGC((i-1)*3+3).gc-RGC((i-1)*3+1).gc;
%         d(1:length(chnsR),1:length(chnsR)) = 0;
%         d(length(chnsR)+1:end,length(chnsR)+1:end) = 0;
%         d(LGC((i-1)*3+3).s == 0) = 0; d(LGC((i-1)*3+1).s == 0) = 0;
%         imagesc(d); n = size(d,1);
%         axis('square'); xlabel('from'); ylabel('to');
%         set(gca,'XTick',1:n); set(gca,'XTickLabel',chnm([chnsR,chnsL])); xtickangle(45);
%         set(gca,'YTick',1:n); set(gca,'YTickLabel',chnm([chnsR,chnsL])); title(['Window ',num2str(i),' Post-Pre'])
%         colorbar;
%         if(diff(caxis)>diff(clim))
%             clim = caxis;
%         end
%         l2r = d(1:length(chnsR),length(chnsR)+1:end); r2l = d(length(chnsR)+1:end,1:length(chnsR));
%         GC_all(end).RDiffPost(i,:) = [mean2(l2r),mean2(r2l)];
%         text(length(chnsR)+ceil(length(chnsL)/2),length(chnsR)+1,num2str(mean2(l2r),'%10.2e'),'color','white');
%         text(length(chnsR)+1,length(chnsR)+ceil(length(chnsL)/2),num2str(mean2(r2l),'%10.2e'),'color','white');
%     end
%     clim = [-max(abs(clim)),max(abs(clim))];
%     for i = 1:6
%         subplot(3,2,i)
%         caxis(clim);
%     end
%     ax = axes; t1 = title('Right Trials'); ax.Visible = 'off'; t1.Visible = 'on';
%     print(RTrial, '-dpsc2', fname, '-append'); close(RTrial);

end

save('F:\Dropbox\repos\abogaard\efetz\U\Mat_Data\GC_time.mat','GC_all','-v7.3');

% %% Overall Plotting (Differences vs RT or Delay)
% LTrial = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); %set(RTrial,'visible','off');
% for i = 1:length(GC_all)
%     if(~strcmp(GC_all(i).Condition(1:6),'Contra'))
%         continue;
%     end
%     disp(GC_all(i).Session);
%     for win = 1:3
%         pre = GC_all(i).LTrial((win-1)*3+1);
%         for epoch = 2:3
%             data = GC_all(i).LTrial((win-1)*3+epoch);
%             gc = data.gc - pre.gc;
%             gc(data.s == 0) = NaN; gc(pre.s == 0) = NaN;
%             l2r = nanmean(nanmean(gc(1:length(chnsR),length(chnsR)+1:end)));%:length(chnsR)
%             r2l = nanmean(nanmean(gc(length(chnsR)+1:end,1:length(chnsR))));
%             disp([num2str(l2r),', ',num2str(r2l)])
%             subplot(3,2,(win-1)*2+epoch-1)
%             hold on; 
% %             scatter(GC_all(i).Delay,l2r,'k') 
%             scatter(GC_all(i).RTL(epoch)-GC_all(i).RTL(1),l2r,'k')
%             hold on; 
% %             scatter(GC_all(i).Delay,r2l,'r') 
%             scatter(GC_all(i).RTL(epoch)-GC_all(i).RTL(1),r2l,'r')
%         end
%     end
% end

