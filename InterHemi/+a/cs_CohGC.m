%% Print out relevant coherence and causality data to see if it's affecting anything.
% For looking into specificities of each day
% Prints a packet with the experiment, neural snips and spectrograms over
% time, coherence and spectral power, and granger causality between
% hemispheres during either trial. Currently uses a full window of 
% -500-1500ms from reaction time, 500ms before reaction time, 500 ms after, 
% and 500-1000ms from reaction time. 

close all; clear; pack

fname = fullfile('F:\Dropbox\repos\abogaard\efetz\U\packets', 'ubi-Coh_GC.ps');
% Delete previous ps file so it doesn't keep appending
if(exist(fname))
    delete(fname);
end

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat')
SessionList = SL;

for day = 1:length(SessionList)
    
    SL = SessionList(day);
    disp(['Session ', SL.Date,': ',num2str(round(day/length(SessionList)*100,1)),'%']);
    if(strcmp(SL.Bad,'1') || isempty(SL.trig1) || strcmp(SL.Condition,'nostim') || strcmp(SL.Condition,'tonic') || strcmp(SL.Condition(end),'R'))
        continue;
    end
    
    if(str2num(SL.Date)>=20170226 && str2num(SL.Date)<=20170306) % days with 100% stim, ignore for the time being
        continue;
    end
    
    %% Load Data
    D = SL.Date;
    S = SL.Session_Guger_Train;
    Session = [char(D),'_',char(S(2))];
    trainfile = [Session,'.f32'];
    gugfile = Session;
    if(~exist(trainfile) || ~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
        continue;
    end
    
%     StimChns = strsplit(SL.Stim_Loc,'/');
%     n = char(StimChns(1));
%     StimChns(2) = cellstr([n(1:4),char(StimChns(2))]);
    
    % down sampling rate
    dwn = 10;
    [accel, trig1, ~, lefttrials, righttrials, ~, lefttrialsuccess, righttrialsuccess] = u.LoadTrainUbi(trainfile,dwn);
    [data, fs, chnm, ~] = u.LoadGug(gugfile, dwn);
    
    if(str2num(SL.Date) == 20170127)
        trig1(1:11) = [];
    end
    
    fo = 60;  q = 35; bw = (fo/(fs/2))/q; % comb notch filter for 60Hz noise
    [B,A] = iircomb(fs/fo,bw,'notch'); % Note type flag 'notch'
    Filter = data;%filtfilt(B,A,double(data));
    
    % Re-calculate reaction time using upsampled, aligned data
    SL.accel_raw_r = accel(:,2);
    SL.accel_raw_l = accel(:,1);
    SL.lefttrials = lefttrials;
    SL.righttrials = righttrials;
    SL.lefttrialsuccess = lefttrialsuccess;
    SL.righttrialsuccess = righttrialsuccess;
    SL.trig1 = trig1;
    SL.fs = fs;
    
    SL.accelfs = fs;
    
    SL = a.AppendReactionTimes(SL);
    SL = a.AppendNormalizedDelay(SL);
    
    %% Define triggers and window
    % lefttrialRT = lefttrials(:,1) + SL.rts_l./1000.*(fs*dwn);
    % lefttrialsuccess(isnan(SL.rts_l)) = 0;
    ofs = 250; %RP starts ~200ms before reaction time calculation
    
    if(isempty(trig1))
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
    remove = [];
    if(strcmp(SL.Condition(1),'I'))
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
    remove = [];
    if(strcmp(SL.Condition(1),'C') && ~strcmp(SL.Condition,'Control'))
        for i = 1:length(trig1)
            norm = abs(trigCondR(:,1)-trig1(i));
            ind = find(norm == min(norm));
            remove(end+1) = ind;
        end
    end
    trigCondR(remove,:) = [];
    trigCondR = trigCondR(:,1);
    
    % window - time in ms after trigger
    % back - time before trigger as a fraction of window
    window = .5; 

    % Define initial channels
    rchn = 6;
    stim = str2num(SL.Stim_Loc(5));
    if(stim == 3), lchn = 24; else, lchn = 22; end;
    
    %% Plot Experiment
    h = u.PlotExperiment(SL);
    print(h, '-dpsc2', fname, '-append'); close(h);
    
    for win = 1:4
        if win == 1
            inds = floor(-window*fs:1:2*window*fs);
        elseif win == 2
            inds = floor(-window*fs:1:0);
        elseif win == 3
            inds = floor(0:1:window*fs);
        else
            inds = floor(window*fs:1:2*window*fs);
        end
        
        %% Plot spectra for the two hemispheres for sanity
        params.tapers = [3,5]; params.Fs = fs; params.fpass = [9,100]; params.trialave = 1;
        movingwin = [0.25,0.01];
        
        clim = [0,0];
        Spec = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(Spec,'visible','off');
        for i = 1:3 %pre/cond/post
            for j = 1:2 %Ltrial/RTrial
                if j == 1
                    if i ==1, trig = trigPreL; elseif i == 2, trig = trigCondL; else, trig = trigPostL; end
                else
                    if i ==1, trig = trigPreR; elseif i == 2, trig = trigCondR; else, trig = trigPostR; end
                end
                trialinds = repmat((trig'), length(inds), 1) + repmat(inds(:), 1, size(trig,1));
                trialinds(:,floor(trialinds(1,:))<=0) = []; trialinds(:,floor(trialinds(end,:))>length(data)) = [];
                for k = 1:2 %Lhemi/Rhemi
                    if k == 1, chn = lchn; else, chn = rchn; end
                    
                    label = [];
                    if(j == 1 && k == 1)
                        if(i == 1), label = 'Pre'; elseif(i==2), label = 'Cond'; else, label = 'Post'; end
                    end
                    
                    D = Filter(:,chn);
                    Snips = u.meanSubtract(D(floor(trialinds)),params);
                    
                    [S,t,f] = mtspecgramc(Snips,movingwin,params);
                    subplot(6,4,(i-1)*8+(j-1)*2+k+4);
                    imagesc(t,f,log(S'))
                    set(gca,'YDir','normal')
                    if(diff(caxis)>diff(clim))
                        clim = caxis;
                    end
                    xl = xlim; ylabel(label);
                    
                    subplot(6,4,(i-1)*8+(j-1)*2+k);
                    hold on; plot((1:size(Snips,1))./fs,Snips); plot((1:size(Snips,1))/(fs),mean(Snips'),'r','LineWidth',2);
                    plot((1:size(Snips,1))/(fs),median(Snips'),'k','LineWidth',2); hold off;
                    ylabel(label); xlim([0,length(Snips)/fs]); 
                    if(strcmp(SL.Animal,'Ubi'))
                        ylim([-5,5])
                    end
                    tr = []; hem = [];
                    if(i==1)
                        if(j==1), tr = 'LTrial '; else, tr = 'RTrial '; end
                        if(k==1), hem = 'LHemi'; else, hem = 'RHemi'; end
                    end
                    title([tr,hem]);
                end
            end
        end
        for i = 1:4
            subplot(6,4,4+i); caxis(clim); subplot(6,4,8+i); caxis(clim); subplot(6,4,12+i); caxis(clim);
        end
        set(0, 'CurrentFigure', Spec); ax = axes; t1 = title([char(chnm(lchn)),', ',char(chnm(rchn)),'_ ']); ax.Visible = 'off'; t1.Visible = 'on';
        print(Spec, '-dpsc2', fname, '-append'); close(Spec);
%         set(Spec,'visible','on');

        %% Get snips from both hemispheres for both trial types across the three parts of the session
        params.tapers = [3,5]; params.Fs = fs; params.fpass = [9,100]; params.trialave = 1;
        movingwin = [0.25,0.01];
        
%         rchn = 5;
%         stim = str2num(SL.Stim_Loc(5));
%         if(stim == 3), lchn = 24; else, lchn = 22; end;
        
        FiltRHemi = Filter(:,rchn);;%filtfilt(bbpf,abpf,double(d));
        FiltLHemi = Filter(:,lchn);;%filtfilt(bbpf,abpf,double(d2));
        
        Coh = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(Coh,'visible','off');
        for i = 1:3
            switch i
                case 1
                    trigL = trigPreL;
                    trigR = trigPreR;
                    label = 'Pre';
                    col = 'k';
                case 2
                    trigL = trigCondL;
                    trigR = trigCondR;
                    label = 'Cond';
                    col = 'b:';
                case 3
                    trigL = trigPostL;
                    trigR = trigPostR;
                    label = 'Post';
                    col = 'r--';
            end
            
            trialindsL = repmat((trigL'), length(inds), 1) + repmat(inds(:), 1, size(trigL,1));
            trialindsR = repmat((trigR'), length(inds), 1) + repmat(inds(:), 1, size(trigR,1));
            trialindsL(:,floor(trialindsL(1,:))<=0) = []; trialindsL(:,floor(trialindsL(end,:))>length(FiltLHemi)) = [];
            trialindsR(:,floor(trialindsR(1,:))<=0) = []; trialindsR(:,floor(trialindsR(end,:))>length(FiltLHemi)) = [];
            
            RSnips_LTrial = u.meanSubtract(FiltRHemi(floor(trialindsL)),params);
            LSnips_LTrial = u.meanSubtract(FiltLHemi(floor(trialindsL)),params);
            
            [C,phi,S12,S1,S2,f,]=coherencyc(RSnips_LTrial,LSnips_LTrial,params);
            
            subplot(4,2,1); hold on;
            plot(f,C,col,'LineWidth',1.5); xlim(params.fpass); hold off;
            title('Left Trials Coherence and Spectra')
            
            subplot(4,2,i*2+1); hold on;
            plot(f,S1,'r','LineWidth',1.5);  plot(f,S2,'k','LineWidth',1.5); hold off;
            xlim(params.fpass); if(i==1), legend(char(chnm(rchn)),char(chnm(lchn))); end
            ylabel(label);
            set(gca, 'YScale', 'log');
            
            RSnips_RTrial = u.meanSubtract(FiltRHemi(floor(trialindsR)),params);
            LSnips_RTrial = u.meanSubtract(FiltLHemi(floor(trialindsR)),params);
            
            [C,phi,S12,S1,S2,f,]=coherencyc(RSnips_RTrial,LSnips_RTrial,params);
            
            subplot(4,2,2); hold on;
            plot(f,C,col,'LineWidth',1.5); xlim(params.fpass); hold off;
            title('Right Trials Coherence and Spectra')
            subplot(4,2,(i+1)*2); hold on;
            plot(f,S1,'r','LineWidth',1.5); plot(f,S2,'k','LineWidth',1.5); xlim(params.fpass); hold off;
            set(gca, 'YScale', 'log'); ylabel(label);
            
        end
        set(0, 'CurrentFigure', Coh); ax = axes; t1 = title('Epochs Coherence_ '); ax.Visible = 'off'; t1.Visible = 'on';
        print(Coh, '-dpsc2', fname, '-append'); close(Coh);
        
        
        %% Granger Causality - pairwise for all left to right channels and vice versa
        % define params
        chnsR = [5,6];%,11,12];
        chnsL = [lchn,lchn+1,[13,14,15,16]+17];
        
        params.regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default) tried LWR, takes much longer but gives same results
        params.icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
        
        params.morder    = fs*0.05; % 50ms, from PNAS paper
        params.momax     = 50;     % maximum model order for model order estimation
        
        params.acmaxlags = [];   % maximum autocovariance lags (empty for automatic calculation)
        
        params.tstat     = 'F';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
        params.alpha     = 0.005;   % significance level for significance test
        params.mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
        
        params.fs        = fs;    % sample rate (Hz)
        params.fres      = [];     % frequency resolution (empty for automatic calculation)
        
        params.nvars     = 2;     % number of variables. only doing pairwise
        
        params.frange    = [9,100]; %range of frequencies to look at
        
        % Left trials
        R2L = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(R2L,'visible','off');
        L2R = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(L2R,'visible','off');
        
        for t = 1:3
            X = [];
            
            switch t
                case 1
                    trig = trigPreL;
                    line = 'k';
                case 2
                    trig = trigCondL;
                    line = 'b:';
                case 3
                    trig = trigPostL;
                    line = 'r--';
            end
            
            trialinds = repmat((trig'), length(inds), 1) + repmat(inds(:), 1, size(trig,1));
            trialinds(:,floor(trialinds(1,:))<=0) = []; trialinds(:,floor(trialinds(end,:))>length(data)) = [];
            
            params.ntrials   = length(trig);     % number of trials
            params.nobs      = length(inds);   % number of observations per trial
            
            for i = 1:length(chnsR)
                d = Filter(:,chnsR(i));%filtfilt(b,a,double(data(:,chnsR(i)))); %comb notch
                Snips = u.meanSubtract(d(floor(trialinds)),params);
                X(1,:,:) = Snips;
                for j = 1:length(chnsL)
                    d2 = Filter(:,chnsL(j));%filtfilt(b,a,double(data(:,chnsL(j)))); %comb notch
                    Snips2 = u.meanSubtract(d2(floor(trialinds)),params);
                    X(2,:,:) = Snips2;
                    
                    disp([num2str(i),num2str(j)]);
                    
                    [gc,f] = a.GC(X,params); % custom helper function for using MVGC toolbox
%                     disp([num2str(round( (((j/J+(i-1))/I+(t-1))/(2*T) + (win-1))/3 * 100,1)),'%']);
                    
                    % plot
                    set(0, 'CurrentFigure', R2L)
                    subplot(length(chnsL),length(chnsR),(j-1)*(length(chnsR))+i)
                    hold on; plot(f,gc(1,:),line,'Linewidth',1.5); xlim(params.frange); hold off;
                    if(j == 1), title(char(chnm(chnsR(i))),'Fontsize',8); end
                    if(i == 1), ylabel(char(chnm(chnsL(j))),'Fontsize',8); end
                    
                    set(0, 'CurrentFigure', L2R)
                    subplot(length(chnsL),length(chnsR),(j-1)*(length(chnsR))+i)
                    hold on; plot(f,gc(2,:),line,'Linewidth',1.5); xlim(params.frange); hold off;
                    if(j == 1), title(char(chnm(chnsR(i))),'Fontsize',8); end
                    if(i == 1), ylabel(char(chnm(chnsL(j))),'Fontsize',8); end
                    
                end
            end
        end
        set(0, 'CurrentFigure', R2L); ax = axes; t1 = title('Left Trial R2L_ '); ax.Visible = 'off'; t1.Visible = 'on';
        set(0, 'CurrentFigure', L2R); ax = axes; t1 = title('Left Trial L2R_ '); ax.Visible = 'off'; t1.Visible = 'on';
        
        print(L2R, '-dpsc2', fname, '-append'); print(R2L, '-dpsc2', fname, '-append'); close(L2R); close(R2L);
        
        % Right trials
        R2L = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(R2L,'visible','off');
        L2R = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(L2R,'visible','off');
        
        for t = 1:3
            X = [];
            
            switch t
                case 1
                    trig = trigPreR;
                    line = 'k';
                case 2
                    trig = trigCondR;
                    line = 'b:';
                case 3
                    trig = trigPostR;
                    line = 'r--';
            end
            
            trialinds = repmat((trig'), length(inds), 1) + repmat(inds(:), 1, size(trig,1));
            trialinds(:,floor(trialinds(1,:))<=0) = []; trialinds(:,floor(trialinds(end,:))>length(data)) = [];
            
            params.ntrials   = length(trig);     % number of trials
            params.nobs      = length(inds);   % number of observations per trial
            
            for i = 1:length(chnsR)
                d = Filter(:,chnsR(i));%filtfilt(b,a,double(data(:,chnsR(i)))); %comb notch
                Snips = u.meanSubtract(d(floor(trialinds)),params);
                X(1,:,:) = Snips;
                for j = 1:length(chnsL)
                    d2 = Filter(:,chnsL(j));%filtfilt(b,a,double(data(:,chnsL(j)))); %comb notch
                    Snips2 = u.meanSubtract(d2(floor(trialinds)),params);
                    X(2,:,:) = Snips2;
                    
                    disp([num2str(i),num2str(j)]);
                    
                    I = length(chnsR); J = length(chnsL); T = 3;
                    [gc,f] = a.GC(X,params); % custom helper function for using MVGC toolbox
%                     disp([num2str(round( (((j/J+(i-1))/I+(t-1))/(2*T) + 0.5 + (win-1))/3 * 100,1)),'%']);
                    
                    % plot
                    set(0, 'CurrentFigure', R2L)
                    subplot(length(chnsL),length(chnsR),(j-1)*(length(chnsR))+i)
                    hold on; plot(f,gc(1,:),line,'Linewidth',1.5); xlim(params.frange); hold off;
                    if(j == 1), title(char(chnm(chnsR(i))),'Fontsize',8); end
                    if(i == 1), ylabel(char(chnm(chnsL(j))),'Fontsize',8); end
                    
                    set(0, 'CurrentFigure', L2R)
                    subplot(length(chnsL),length(chnsR),(j-1)*(length(chnsR))+i)
                    hold on; plot(f,gc(2,:),line,'Linewidth',1.5); xlim(params.frange); hold off;
                    if(j == 1), title(char(chnm(chnsR(i))),'Fontsize',8); end
                    if(i == 1), ylabel(char(chnm(chnsL(j))),'Fontsize',8); end
                    
                end
            end
        end
        set(0, 'CurrentFigure', R2L); ax = axes; t1 = title('Right Trial R2L_ '); ax.Visible = 'off'; t1.Visible = 'on';
        set(0, 'CurrentFigure', L2R); ax = axes; t1 = title('Right Trial L2R_ '); ax.Visible = 'off'; t1.Visible = 'on';
        
        print(L2R, '-dpsc2', fname, '-append'); print(R2L, '-dpsc2', fname, '-append'); close(L2R); close(R2L);
    end
end
