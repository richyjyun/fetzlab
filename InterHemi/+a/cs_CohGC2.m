%% Print out relevant coherence and causality data to see if it's affecting anything.
% For looking into comparisons between days. 
% Trying to distill more of the data into simple points.
% Creates a packet of the experiment, snips and spectrograms over time,
% coherence and spectral power, and comparisons of integral in the
% beta/gamma range across the epochs. 
close all; clear; pack

tic
fname = fullfile('F:\Dropbox\repos\abogaard\efetz\U\packets', 'ubi-Coh_GC2.ps');
% Delete previous ps file so it doesn't keep appending
if(exist(fname))
    delete(fname);
end

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat')
SessionList = SL; 
B = []; G = []; % Beta causality (day,window,epoch,trial type,l2r vs r2l), Gamma causality
Sessions = cell(0); % Sessions
RT = []; % reaction time (day,epoch,trial type)
delays = []; % normalized delay, not recalculated, just taken from Session list. 
Contra = [];
power_all = [];

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
    
    Sessions(end+1,1) = cellstr(Session);
    delays(end+1) = SL.NormDelay;
    if(strcmp(SL.Condition(1),'C'))
        Contra(end+1) = 1;
    else
        Contra(end+1) = 0;
    end
    
    StimChns = strsplit(SL.Stim_Loc,'/');
    n = char(StimChns(1));
    StimChns(2) = cellstr([n(1:4),char(StimChns(2))]);
    
    % down sampling rate
    dwn = 10;
    [accel, trig1, ~, lefttrials, righttrials, ~, lefttrialsuccess, righttrialsuccess] = u.LoadTrainUbi(trainfile,dwn);
    [data, fs, chnm, ~] = u.LoadGug(gugfile, dwn);
        
    if(str2num(SL.Date) == 20170127)
        trig1(1:11) = [];
    end
    
    fo = 60;  q = 35; bw = (fo/(fs/2))/q; % comb notch filter for 60Hz noise
    [BB,A] = iircomb(fs/fo,bw,'notch'); % Note type flag 'notch'
    Filter = data;%filtfilt(BB,A,double(data));
    
    % Re-calculate reaction time using upsampled, aligned data
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
    
    %% Define triggers and window
    % lefttrialRT = lefttrials(:,1) + SL.rts_l./1000.*(fs*dwn);
    % lefttrialsuccess(isnan(SL.rts_l)) = 0;
    ofs = 250; %RP starts ~250ms before reaction time calculation (at least for Ubi)
    
    lefttrials = SL.lefttrials+(SL.rts_l-ofs)/1000*SL.fs; lefttrials = lefttrials(~isnan(SL.rts_l)& lefttrialsuccess,:);
    trigPreL = lefttrials((lefttrials(:,1)<trig1(1)));
    trigPostL = lefttrials((lefttrials(:,1)>trig1(end)));
    trigCondL = lefttrials(lefttrials(:,1)>=trig1(1) & lefttrials(:,1)<=trig1(end),:);
    RT(end+1,1,1) = nanmedian(SL.rts_l(SL.lefttrials(:,1)>=trig1(1) & SL.lefttrials(:,1)<=trig1(end)))...
        -nanmedian(SL.rts_l(SL.lefttrials(:,1)<trig1(1)));
    RT(end,2,1) = nanmedian(SL.rts_l(SL.lefttrials(:,1)>trig1(end)))...
        -nanmedian(SL.rts_l(SL.lefttrials(:,1)<trig1(1)));
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
    RT(end,1,2) = nanmedian(SL.rts_r(SL.righttrials(:,1)>=trig1(1) & SL.righttrials(:,1)<=trig1(end)))...
        -nanmedian(SL.rts_r(SL.righttrials(:,1)<trig1(1)));
    RT(end,2,2) = nanmedian(SL.rts_r(SL.righttrials(:,1)>trig1(end)))...
        -nanmedian(SL.rts_r(SL.righttrials(:,1)<trig1(1)));
    remove = [];
    if(strcmp(SL.Condition(1),'C'))
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
    
    beta = [14,30]; 
    gamma = [85,100];
    
    %% Plot Experiment
    h = u.PlotExperiment(SL);
    print(h, '-dpsc2', fname, '-append'); close(h);
    
    LTrial = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(LTrial,'visible','off');
    RTrial = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(RTrial,'visible','off');

    for win = 1:4
        if win == 1
            inds = floor(-window*fs:1:2*window*fs);
            line = 'k';
        elseif win == 2
            inds = floor(-window*fs:1:0);
            line = 'b:';
        elseif win == 3
            inds = floor(0:1:window*fs);
            line = 'r--';
        else
            inds = floor(window*fs:1:2*window*fs);
            line = 'g-.';
        end
        
        %% Plot spectra for the two hemispheres for sanity
        params.tapers = [3,5]; params.Fs = fs; params.fpass = [9,100]; params.trialave = 1;
        movingwin = [0.25,0.01];
        rchn = 6;
        stim = str2num(SL.Stim_Loc(5));
        if(stim == 3), lchn = 24; else, lchn = 22; end;
        
        clim = [0,0];
        Spec = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(Spec,'visible','off');
        for i = 1:3 %pre/cond/post
            for j = 1:2 %Ltrial/RTrial
                if j == 1
                    if i ==1, trig = trigPreL; elseif i == 2, trig = trigCondL; else, trig = trigPostL; end
                else
                    if i ==1, trig = trigPreR; elseif i == 2, trig = trigCondR; else, trig = trigPostR; end
                end
                trialinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
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

        %% Get spectral power and coherence for both trial types across the three epochs of the session
        params.tapers = [3,5]; params.Fs = fs; params.fpass = [9,100]; params.trialave = 1;
        movingwin = [0.25,0.01];
        
        rchn = 5;
        stim = str2num(SL.Stim_Loc(5));
        if(stim == 3), lchn = 24; else, lchn = 22; end;
        
        FiltRHemi = Filter(:,rchn);;%filtfilt(bbpf,abpf,double(d));
        FiltLHemi = Filter(:,lchn);;%filtfilt(bbpf,abpf,double(d2));
        
        betaP = []; gammaP = [];
        
        Coh = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(Coh,'visible','off');
        Phase = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(Phase,'visible','off');
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
            
            trialindsL = repmat(trigL', length(inds), 1) + repmat(inds(:), 1, size(trigL,1));
            trialindsR = repmat(trigR', length(inds), 1) + repmat(inds(:), 1, size(trigR,1));
            trialindsL(:,floor(trialindsL(1,:))<=0) = []; trialindsL(:,floor(trialindsL(end,:))>length(FiltLHemi)) = [];
            trialindsR(:,floor(trialindsR(1,:))<=0) = []; trialindsR(:,floor(trialindsR(end,:))>length(FiltLHemi)) = [];
            
            RSnips_LTrial = u.meanSubtract(FiltRHemi(floor(trialindsL)),params);
            LSnips_LTrial = u.meanSubtract(FiltLHemi(floor(trialindsL)),params);
            
            [C,phi,S12,S1,S2,f,]=coherencyc(RSnips_LTrial,LSnips_LTrial,params);
            
            [~,ind] = max(C(f>=beta(1) & f<=beta(2)));
            lbound = find(f>=beta(1),1)-1;
            betaP(1,i,1) = f(lbound+ind); %frequency of max coherence
            betaP(1,i,2) = phi(lbound+ind)/(2*pi*betaP(1,i,1))*1000; %phase shift in ms
            
            [~,ind] = max(C(f>=gamma(1) & f<=gamma(2)));
            lbound = find(f>=gamma(1),1)-1;
            gammaP(1,i,1) = f(lbound+ind); %frequency of max coherence
            gammaP(1,i,2) = phi(lbound+ind)/(2*pi*gammaP(1,i,1))*1000; %phase shift in ms
            
            set(0,'CurrentFigure',Coh);
            subplot(4,2,1); hold on;
            plot(f,C,col,'LineWidth',1.5); xlim(params.fpass); hold off;
            title('Left Trials Coherence and Spectra')            
            subplot(4,2,i*2+1); hold on;
            plot(f,S1,'r','LineWidth',1.5);  plot(f,S2,'k','LineWidth',1.5); hold off;
            xlim(params.fpass); if(i==1), legend(char(chnm(rchn)),char(chnm(lchn))); end
            ylabel(label);
            set(gca, 'YScale', 'log');
            
            set(0,'CurrentFigure',Phase);
            subplot(3,2,i*2-1); hold on; plot(f,phi./(2*pi*f')); plot(f,(phi+2*pi)./(2*pi*f'));
            
            
            RSnips_RTrial = u.meanSubtract(FiltRHemi(floor(trialindsR)),params);
            LSnips_RTrial = u.meanSubtract(FiltLHemi(floor(trialindsR)),params);
            
            [C,phi,S12,S1,S2,f,]=coherencyc(RSnips_RTrial,LSnips_RTrial,params);
            
            [~,ind] = max(C(f>=beta(1) & f<=beta(2)));
            lbound = find(f>=beta(1),1)-1;
            betaP(2,i,1) = f(lbound+ind); %frequency of max coherence
            betaP(2,i,2) = phi(lbound+ind)/(2*pi*betaP(2,i,1))*1000; %phase shift in ms
            [~,ind] = max(C(f>=gamma(1) & f<=gamma(2)));
            lbound = find(f>=gamma(1),1)-1;
            gammaP(2,i,1) = f(lbound+ind); %frequency of max coherence
            gammaP(2,i,2) = phi(lbound+ind)/(2*pi*gammaP(2,i,1))*1000; %phase shift in ms
            
            set(0,'CurrentFigure',Coh);
            subplot(4,2,2); hold on;
            plot(f,C,col,'LineWidth',1.5); xlim(params.fpass); hold off;
            title('Right Trials Coherence and Spectra')
            subplot(4,2,(i+1)*2); hold on;
            plot(f,S1,'r','LineWidth',1.5); plot(f,S2,'k','LineWidth',1.5); xlim(params.fpass); hold off;
            ylabel(label);
            set(gca, 'YScale', 'log');
                       
            set(0,'CurrentFigure',Phase);
            subplot(3,2,i*2); plot(f,phi./(2*pi*f'));
            
        end
        set(0, 'CurrentFigure', Coh); subplot(4,2,1); yl = ylim; hold on;
        text(30,yl(2)*2/3,{[num2str(round(betaP(1,1,1),1)),'Hz; ',num2str(round(betaP(1,1,2),1)),'ms   ',num2str(round(gammaP(1,1,1),1)),'Hz; ',num2str(round(gammaP(1,1,2),1)),'ms'],...
            [num2str(round(betaP(1,2,1),1)),'Hz; ',num2str(round(betaP(1,2,2),1)),'ms   ',num2str(round(gammaP(1,2,1),1)),'Hz; ',num2str(round(gammaP(1,2,2),1)),'ms'],...
            [num2str(round(betaP(1,3,1),1)),'Hz; ',num2str(round(betaP(1,3,2),1)),'ms   ',num2str(round(gammaP(1,3,1),1)),'Hz; ',num2str(round(gammaP(1,3,2),1)),'ms']}); hold off;
        
        set(0, 'CurrentFigure', Coh); subplot(4,2,2); yl = ylim; hold on;
        text(30,yl(2)*2/3,{[num2str(round(betaP(2,1,1),1)),'Hz; ',num2str(round(betaP(2,1,2),1)),'ms   ',num2str(round(gammaP(2,1,1),1)),'Hz; ',num2str(round(gammaP(2,1,2),1)),'ms'],...
            [num2str(round(betaP(2,2,1),1)),'Hz; ',num2str(round(betaP(2,2,2),1)),'ms   ',num2str(round(gammaP(2,2,1),1)),'Hz; ',num2str(round(gammaP(2,2,2),1)),'ms'],...
            [num2str(round(betaP(2,3,1),1)),'Hz; ',num2str(round(betaP(2,3,2),1)),'ms   ',num2str(round(gammaP(2,3,1),1)),'Hz; ',num2str(round(gammaP(2,3,2),1)),'ms']}); hold off;

        set(Coh,'visible','on')
        
        set(0, 'CurrentFigure', Coh); ax = axes; t1 = title('Epochs Coherence_ '); ax.Visible = 'off'; t1.Visible = 'on';
        print(Coh, '-dpsc2', fname, '-append'); close(Coh);
        
        
        %% Granger Causality - pairwise for all left to right channels and vice versa
        % define params
        chnsR = [5,6];%,11,12];
        chnsL = [lchn,lchn+1,[13,14,15,16]+17];
        params.beta      = beta;
        params.gamma     = gamma;
        
        params.regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default) tried LWR, takes much longer but gives same results
        params.icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
        
        params.morder    = fs*0.05; % 50ms, from PNAS paper
        params.momax     = 50;     % maximum model order for model order estimation
        
        params.acmaxlags = [];   % maximum autocovariance lags (empty for automatic calculation)
        
        params.statt     = 'F';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
        params.alpha     = 0.005;   % significance level for significance test
        params.mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
        
        params.fs        = fs;    % sample rate (Hz)
        params.fres      = [];     % frequency resolution (empty for automatic calculation)
        
        params.nvars     = 2;     % number of variables. only doing pairwise
        
        params.frange    = [9,100]; %range of frequencies to look at
        
        % Left trials
        power = a.doGC(params,trigPreL,trigCondL,trigPostL,inds,Filter,chnsL,chnsR,chnm,LTrial,line);
        power_all(1,win,:,:,:,:,:) = power;
        
        set(0,'CurrentFigure',LTrial);
        if win == 4
            ax = axes; t1 = title('Left Trial - Beta L2R & R2L, Gamma L2R & R2L_ '); ax.Visible = 'off'; t1.Visible = 'on';
            print(LTrial, '-dpsc2', fname, '-append'); close(LTrial);
        end
        
        % need to add power to B and G for overall figure
        if(win == 1)
            B(end+1,win,1,1,1) =  mean2(squeeze(power(2,1,:,:,1))) - mean2(squeeze(power(1,1,:,:,1)));  %cond vs pre, LTrial, l2r
        else
            B(end,win,1,1,1) =  mean2(squeeze(power(2,1,:,:,1))) - mean2(squeeze(power(1,1,:,:,1)));  %cond vs pre, LTrial, l2r
        end
        B(end,win,1,1,2) =  mean2(squeeze(power(2,1,:,:,2))) - mean2(squeeze(power(1,1,:,:,2)));
        B(end,win,2,1,1) = mean2(squeeze(power(3,1,:,:,1))) - mean2(squeeze(power(1,1,:,:,1))); 
        B(end,win,2,1,2) =  mean2(squeeze(power(3,1,:,:,2))) - mean2(squeeze(power(1,1,:,:,2)));
        if(win == 1)
            G(end+1,win,1,1,1) =  mean2(squeeze(power(2,2,:,:,1))) - mean2(squeeze(power(1,2,:,:,1)));  %cond vs pre, LTrial, l2r
        else
            G(end,win,1,1,1) =  mean2(squeeze(power(2,2,:,:,1))) - mean2(squeeze(power(1,2,:,:,1)));  %cond vs pre, LTrial, l2r
        end
        G(end,win,1,1,2) =  mean2(squeeze(power(2,2,:,:,2))) - mean2(squeeze(power(1,2,:,:,2)));
        G(end,win,2,1,1) = mean2(squeeze(power(3,2,:,:,1))) - mean2(squeeze(power(1,2,:,:,1)));
        G(end,win,2,1,2) =  mean2(squeeze(power(3,2,:,:,2))) - mean2(squeeze(power(1,2,:,:,2)));
        
        
        % Right trials
        power = a.doGC(params,trigPreR,trigCondR,trigPostR,inds,Filter,chnsL,chnsR,chnm,RTrial,line);
        power_all(2,win,:,:,:,:,:) = power;        
        
        if win == 4
            ax = axes; t1 = title('Right Trial - Beta L2R & R2L, Gamma L2R & R2L_ '); ax.Visible = 'off'; t1.Visible = 'on';
            print(RTrial, '-dpsc2', fname, '-append'); close(RTrial);
        end
        % need to add power to B and G for overall figure
        B(end,win,1,2,1) =  mean2(squeeze(power(2,1,:,:,1))) - mean2(squeeze(power(1,1,:,:,1)));  %cond vs pre, LTrial, l2r
        B(end,win,1,2,2) =  mean2(squeeze(power(2,1,:,:,2))) - mean2(squeeze(power(1,1,:,:,2)));
        B(end,win,2,2,1) = mean2(squeeze(power(3,1,:,:,1))) - mean2(squeeze(power(1,1,:,:,1))); 
        B(end,win,2,2,2) =  mean2(squeeze(power(3,1,:,:,2))) - mean2(squeeze(power(1,1,:,:,2)));
        G(end,win,1,2,1) =  mean2(squeeze(power(2,2,:,:,1))) - mean2(squeeze(power(1,2,:,:,1)));  %cond vs pre, LTrial, l2r
        G(end,win,1,2,2) =  mean2(squeeze(power(2,2,:,:,2))) - mean2(squeeze(power(1,2,:,:,2)));
        G(end,win,2,2,1) = mean2(squeeze(power(3,2,:,:,1))) - mean2(squeeze(power(1,2,:,:,1)));
        G(end,win,2,2,2) =  mean2(squeeze(power(3,2,:,:,2))) - mean2(squeeze(power(1,2,:,:,2)));
    end
%     SLnew(end+1) = SL;
end

% B = []; G = []; % Beta causality (day,window,epoch,trial type,l2r vs r2l), Gamma causality
% Sessions = cell(0); % Sessions
% RT = []; % reaction time (day,epoch,trial type)

%% Plotting comparison between days (try it with both RT and normalized delay?)
for j = 1:2
    if j == 1
        cond = Contra;
    else
        cond = ~Contra;
    end
    
    Beta = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(Beta,'visible','off');
    Gamma = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(Gamma,'visible','off');
    for i = 1:size(B,2)*size(B,3)*size(B,4) %number of windows * number of epochs * number of trial types
        win = floor((i-1)/(size(B,3)*size(B,4)))+1;
        epoch = mod(floor((i-1)/size(B,4)),size(B,3))+1;
        trial = mod(i-1,size(B,4))+1;
        
        set(0, 'CurrentFigure', Beta);
        subplot(size(B,2),size(B,3)*size(B,4),i);
        hold on; scatter(RT(cond,epoch,trial),B(cond,win,epoch,trial,1),'k'); scatter(RT(cond,epoch,trial),B(cond,win,epoch,trial,2),'r'); hold off;
        if(mod(i-1,size(B,3)*size(B,4)) == 0)
            if(win == 1)
                ylabel('Whole');
            elseif(win==2)
                ylabel('-500~0ms');
            elseif(win==3)
                ylabel('0~500ms');
            else
                ylabel('500~1000ms');
            end
        end
        
        if(win == 1)
            if(epoch == 1 && trial == 1)
                title('Cond Left Trial')
            elseif(epoch == 1 && trial == 2)
                title('Cond Right Trial')
            elseif(epoch == 2 && trial == 1)
                title('Post Left Trial')
            else
                title('Post Right Trial')
            end
        end
        
        set(0, 'CurrentFigure', Gamma);
        subplot(size(B,2),size(B,3)*size(B,4),i);
        hold on; scatter(RT(cond,epoch,trial),G(cond,win,epoch,trial,1),'k'); scatter(RT(cond,epoch,trial),G(cond,win,epoch,trial,2),'r'); hold off;
        if(mod(i-1,size(B,3)*size(B,4)) == 0)
            if(win == 1)
                ylabel('Whole');
            elseif(win==2)
                ylabel('-500~0ms');
            elseif(win==3)
                ylabel('0~500ms');
            else
                ylabel('500~1000ms');
            end
        end
        
        if(win == 1)
            if(epoch == 1 && trial == 1)
                title('Cond Left Trial')
            elseif(epoch == 1 && trial == 2)
                title('Cond Right Trial')
            elseif(epoch == 2 && trial == 1)
                title('Post Left Trial')
            else
                title('Post Right Trial')
            end
        end
    end
    
    set(0, 'CurrentFigure', Beta);
    ax = axes; t1 = title('Beta_ '); ax.Visible = 'off'; t1.Visible = 'on';
    print(Beta, '-dpsc2', fname, '-append'); close(Beta);
    
    set(0, 'CurrentFigure', Gamma);
    ax = axes; t1 = title('Gamma_ '); ax.Visible = 'off'; t1.Visible = 'on';
    print(Gamma, '-dpsc2', fname, '-append'); close(Gamma);
end

toc

