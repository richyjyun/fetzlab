load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat')

%% Choosing a good day
i = 9;
MRCA = 1;
temp = SL;
SL = SL(i);

%% Load Data
D = SL.Date;
S = SL.Session_Guger_Train;
Session = [char(D),'_',char(S(2))];
trainfile = [Session,'.f32'];
gugfile = Session;
if(~exist(trainfile) || ~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
    keyboard;
end

StimChns = strsplit(SL.Stim_Loc,'/');
n = char(StimChns(1));
StimChns(2) = cellstr([n(1:4),char(StimChns(2))]);

% down sampling rate
dwn = 10;
[accel, trig1, ~, lefttrials, righttrials, ~, lefttrialsuccess, righttrialsuccess] = u.LoadTrainUbi(trainfile,dwn);
[data, fs, chnm, ~] = u.LoadGug(gugfile, dwn);

fo = 60;  q = 35; bw = (fo/(fs/2))/q; % comb notch filter for 60Hz noise
[B,A] = iircomb(fs/fo,bw,'notch'); % Note type flag 'notch'
Filter = filtfilt(B,A,double(data)); 

% Account for different number of trials between synced and unsyned data
temp = lefttrials*1000/(fs*dwn);
check = zeros(length(lefttrials),1);
for i = 1:length(temp)
    norm = abs(temp(i,1)-SL.lefttrials(:,1));
    ind = find(norm == min(norm));
    check(ind) = 1;
end
lefttrials = lefttrials(find(check),:);
lefttrialsuccess = lefttrialsuccess(find(check));

temp = righttrials*1000/(fs*dwn);
check = zeros(1,length(righttrials));
for i = 1:length(temp)
    norm = abs(temp(i,1)-SL.righttrials(:,1));
    ind = find(norm == min(norm));
    check(ind) = 1;
end
righttrials = righttrials(find(check),:);
righttrialsuccess = righttrialsuccess(find(check));


%% Define triggers and window
% lefttrialRT = lefttrials(:,1) + SL.rts_l./1000.*(fs*dwn);   
% lefttrialsuccess(isnan(SL.rts_l)) = 0;
trigPreL = lefttrials((lefttrials(:,1)-15<trig1(1)) & lefttrialsuccess);
trigPostL = lefttrials((lefttrials(:,1)+15>trig1(end)) & lefttrialsuccess);
trigCondL = lefttrials(lefttrials(:,1)+15>=trig1(1) & lefttrials(:,1)<=trig1(end)-15 & lefttrialsuccess,:);
remove = [];
for i = 1:length(trig1)
    stim = trigCondL(:,1)-15 < trig1(i) & trigCondL(:,2) + 15 > trig1(i);
    if(any(stim))
        remove(end+1) = find(stim,1);
    end
end
trigCondL(remove,:) = [];
trigCondL = trigCondL(:,1);%+SL.rts_l(length(trigPreL)+1:length(trigPreL)+length(trigCondL));   
shift =0;% nanmedian(SL.rts_l);

trigPreR = righttrials(righttrials(:,1)-15<trig1(1) & righttrialsuccess,1); trigPreR(end) = [];
trigPostR = righttrials(righttrials(:,1)+15>trig1(end) & righttrialsuccess,1);
trigCondR = righttrials(righttrials(:,1)+15>=trig1(1) & righttrials(:,1)<=trig1(end)-15 & righttrialsuccess,:);
remove = [];
for i = 1:length(trig1)
    stim = trigCondR(:,1)-15 < trig1(i) & trigCondR(:,2) + 15 > trig1(i);
    if(any(stim))
        remove(end+1) = find(stim,1);
    end
end
trigCondR(remove,:) = [];
trigCondR = trigCondR(:,1);
   
% window - time in ms after trigger
% back - time before trigger as a fraction of window
window = 1500;
back = 3;

%% Get neural data snips and Moving window spectrum 
trig = trigPreR;
inds = floor(-window*(1/back)*fs/1000:1:window*fs/1000);
trialinds = repmat((trig')./dwn, length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds(:,floor(trialinds(1,:))<=0) = []; trialinds(:,floor(trialinds(end,:))>length(data)) = [];

rchn = 4;
stream = data(:,rchn);
NSnips = u.meanSubtract(stream(floor(trialinds)));

% Define parameters for Chronux
params.tapers = [3,5]; params.Fs = fs; params.fpass = [9,100]; params.trialave = 1;
movingwin = [0.25,0.01];

% Moving window spectrum
[S,t,f] = mtspecgramc(NSnips,movingwin,params);

% [S,f] = mtspectrumc(NSnips,params);

% Plotting for sanity
figure;
subplot(2,2,3)
imagesc(t,f,log(S'))
set(gca,'YDir','normal')
clim1 = caxis;
xl = xlim;

subplot(2,2,1)
plot((1:length(NSnips))./fs,NSnips)
xlim(xl);
title([char(chnm(rchn)),' Right Hemisphere'])

% Checking right hemisphere as well
% Filter = filtfilt(bbpf,abpf,double(d2));
lchn = 22;
stream = Filter(:,lchn);
NSnips = u.meanSubtract(stream(floor(trialinds)));
params.tapers = [3,5]; params.Fs = fs; params.fpass = [9,100]; params.trialave = 1;
movingwin = [0.25,0.01];

% Moving window spectrum
[S,t,f] = mtspecgramc(NSnips,movingwin,params);

% Plotting for sanity
subplot(2,2,4)
imagesc(t,f,log(S'))
clim2 = caxis;
set(gca,'YDir','normal')
xl = xlim;

subplot(2,2,2)
plot((1:length(NSnips))/fs,NSnips)
xlim(xl);
title([char(chnm(lchn)) ' Left Hemisphere'])

if(diff(clim1)>diff(clim2)), clim = clim1; else, clim=clim2; end;
subplot(2,2,3); caxis(clim); subplot(2,2,4); caxis(clim); 

%% Get snips from both hemispheres for both trial types across the three parts of the session
params.tapers = [3,5]; params.Fs = fs; params.fpass = [9,100]; params.trialave = 1;
movingwin = [0.25,0.01];

rchn = 5; lchn = 23;
FiltRHemi = Filter(:,rchn);;%filtfilt(bbpf,abpf,double(d));
FiltLHemi = Filter(:,lchn);;%filtfilt(bbpf,abpf,double(d2));

figure;
for i = 1:3
    switch i
        case 1 
            trigL = trigPreL;
            trigR = trigPreR;
            col = 'k';
        case 2
            trigL = trigCondL;
            trigR = trigCondR;
            col = 'b:';
        case 3
            trigL = trigPostL;
            trigR = trigPostR;
            col = 'r--';
    end
    
    inds = floor((-window*(1/back))*fs/1000:1:window*fs/1000);
    trialindsL = repmat((trigL')./dwn, length(inds), 1) + repmat(inds(:), 1, size(trigL,1));
    trialindsR = repmat((trigR')./dwn, length(inds), 1) + repmat(inds(:), 1, size(trigR,1));
    trialindsL(:,floor(trialindsL(1,:))<=0) = []; trialindsL(:,floor(trialindsL(end,:))>length(FiltLHemi)) = [];
    trialindsR(:,floor(trialindsR(1,:))<=0) = []; trialindsR(:,floor(trialindsR(end,:))>length(FiltLHemi)) = [];

    RSnips_LTrial = u.meanSubtract(FiltRHemi(floor(trialindsL)));
    LSnips_LTrial = u.meanSubtract(FiltLHemi(floor(trialindsL)));
    
    [C,phi,S12,S1,S2,f,]=coherencyc(RSnips_LTrial,LSnips_LTrial,params);
    
    subplot(4,2,1); hold on;
    plot(f,C,col,'LineWidth',1.5); xlim(params.fpass); hold off;
    title('Left Trials')
    
    subplot(4,2,i*2+1); hold on;
    plot(f,S1,'r','LineWidth',1.5);  plot(f,S2,'k','LineWidth',1.5); hold off;
    xlim(params.fpass); if(i==1), legend('R','L'); end
%     set(gca, 'YScale', 'log');    
        
    RSnips_RTrial = u.meanSubtract(FiltRHemi(floor(trialindsR)));
    LSnips_RTrial = u.meanSubtract(FiltLHemi(floor(trialindsR)));
   
    [C,phi,S12,S1,S2,f,]=coherencyc(RSnips_RTrial,LSnips_RTrial,params);
    
    subplot(4,2,2); hold on;
    plot(f,C,col,'LineWidth',1.5); xlim(params.fpass); hold off;
    title('Right Trials')
    subplot(4,2,(i+1)*2); hold on;
    plot(f,S1,'r','LineWidth',1.5); plot(f,S2,'k','LineWidth',1.5); xlim(params.fpass); hold off;
%     set(gca, 'YScale', 'log');
    
end


%% Coherence spectra over time between verious channels across hemispheres
% Define params
chnsR = [5,6]; 
chnsL = [7,8,13,14,15,16]+17; 
params.tapers = [3,5]; params.Fs = fs; params.fpass = [9,100]; params.trialave = 1;
movingwin = [0.25,0.01];

trig = trigPreL;

inds = floor((-window*(1/back)-shift)*fs/1000:1:(window-shift)*fs/1000);
trialinds = repmat((trig')./dwn, length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds(:,floor(trialinds(1,:))<=0) = []; trialinds(:,floor(trialinds(end,:))>length(Filter)) = [];

figure;
clim = [0,0];
for i = 1:length(chnsR)
    d = Filter(:,chnsR(i));%filtfilt(b,a,double(data(:,chnsR(i)))); %comb notch
    Snips = u.meanSubtract(d(floor(trialinds)));

    for j = 1:length(chnsL)
        d2 = Filter(:,chnsL(j));%filtfilt(b,a,double(data(:,chnsL(j)))); %comb notch
        Snips2 = u.meanSubtract(d2(floor(trialinds)));

        disp([i,j]) % just for checking where it's at
        
        [C,phi,S12,S1,S2,t,f] = cohgramc(Snips,Snips2,movingwin,params);
        
        % plot
        subplot(length(chnsR)*2,length(chnsL),(i-1)*(length(chnsL))+j+12);
        imagesc(t,f,log(squeeze(C')))
        if(diff(caxis)>diff(clim)), clim = caxis; end
        set(gca,'YDir','normal')
        if(j == 1), ylabel(char(chnm(chnsR(i))),'Fontsize',8); end
        if(i == 1), title(char(chnm(chnsL(j))),'Fontsize',8); end
        
    end
end

for i=1:length(chnsR)*2*length(chnsL), subplot(length(chnsR)*2,length(chnsL),i); caxis(clim); end
ax = axes; t1 = title('Top: LTrials, Bot: RTrials_ '); ax.Visible = 'off'; t1.Visible = 'on';


%% Plot coherence between multiple electrodes across hemispheres
chnsR = [5,6]; 
chnsL = [7,8,13,14,15,16]+17; 
params.tapers = [3,5]; params.Fs = fs; params.fpass = [9,100]; params.trialave = 1;
movingwin = [0.25,0.01];

trig = trigPreL;

inds = floor((-window*(1/back)-shift)*fs/1000:1:(window-shift)*fs/1000);
trialinds = repmat((trig')./dwn, length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds(:,floor(trialinds(1,:))<=0) = []; trialinds(:,floor(trialinds(end,:))>length(Filter)) = [];

div = 2;  % number of devisions over time
Coh = figure;
range = [];
for i = 1:length(chnsR)
    d = Filter(:,chnsR(i));
    Snips = u.meanSubtract(d(floor(trialinds)));
    for j = 1:length(chnsL)
        d2 = Filter(:,chnsL(j));
        Snips2 = u.meanSubtract(d2(floor(trialinds)));
        for k = 1:div
            Coh;
            Lbound = floor((k-1)/div*length(Snips)+1);
            Rbound = floor(k/div*length(Snips));
            [C,phi,~,~,~,f] = coherencyc(Snips(Lbound:Rbound,:),...
                Snips2(Lbound:Rbound,:),params);
%             [C,phi,~,~,~,f] = coherencyc(Snips,...
%                 Snips2,params);
            subplot(length(chnsR),length(chnsL),(i-1)*(length(chnsL))+j)
            switch k
                case 1
                    line = 'k';
                case 2
                    line = 'b:';
                case 3
                    line = 'r--';
                case 4
                    line = 'g-.';
            end
            hold on
            plot(f,C,line,'Linewidth',1.5); xlim([min(f),max(f)]);
            hold off
            range(end+1,:) = [min(min(C)),max(max(C))];
            if(j == 1)
                ylabel(char(chnm(chnsR(i))),'Fontsize',8)
            end
            if(i == 1)
                title(char(chnm(chnsL(j))),'Fontsize',8)
            end
        end
        if(i==1&&j==1), legend('1','2','3','4'); end
        yl = ylim; hold on; plot([20,20],yl,'k'); ylim(yl); hold off;
    end
end


%% Plot coherence points at maximum amplitude of beta/gamma across hemispheres
chnsR = [5,6]; 
chnsL = [7,8,13,14,15,16]+17; 
params.tapers = [3,5]; params.Fs = fs; params.fpass = [9,100]; params.trialave = 1;
movingwin = [0.25,0.01];

beta = [14,30];
gamma = [85,100];

div = 4;
Beta = figure; set(Beta,'visible','off'); 
Gamma = figure; set(Gamma,'visible','off'); 

for t = 1:3
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
    
    inds = floor(-window*(1/back)*fs/1000:1:window*fs/1000);
    trialinds = repmat((trig')./dwn, length(inds), 1) + repmat(inds(:), 1, size(trig,1));
    trialinds(:,floor(trialinds(1,:))<=0) = []; trialinds(:,floor(trialinds(end,:))>length(data)) = [];
    
    for i = 1:length(chnsR)
        d = data(:,chnsR(i));%filtfilt(b,a,double(data(:,chnsR(i)))); %comb notch
        Snips = u.meanSubtract(d(floor(trialinds)));
        for j = 1:length(chnsL)
            d2 = data(:,chnsL(j));%filtfilt(b,a,double(data(:,chnsL(j)))); %comb notch
            Snips2 = u.meanSubtract(d2(floor(trialinds)));
            B = []; G = []; BPhi = []; GPhi = []; Bf = []; Gf = [];
            for k = 1:div
                Lbound = floor((k-1)/div*length(Snips)+1);
                Rbound = floor(k/div*length(Snips));
                [C,phi,~,~,~,f] = coherencyc(Snips(Lbound:Rbound,:),...
                    Snips2(Lbound:Rbound,:),params);
                
                B(end+1) = max(C(f>=beta(1)&f<=beta(2)));
                G(end+1) = max(C(f>=gamma(1)&f<=gamma(2)));
                Bind = find(C==B(end),1); Gind = find(C==G(end),1);
                BPhi(end+1) = round(phi(Bind)/(2*pi*f(Bind))*1000);
                GPhi(end+1) = round(phi(Gind)/(2*pi*f(Gind))*1000);
                Bf(end+1) = f(Bind);
                Gf(end+1) = f(Gind);
            end
            set(0, 'CurrentFigure', Beta);
            subplot(length(chnsR)*2,length(chnsL),(i-1)*(length(chnsL))+j);
            %             set(gca,'XTick',1:div); set(gca,'XTickLabel',num2cell(Bf));
            hold on; plot(B,line); text(1:div,B,num2cell(BPhi),'Fontsize',7); hold off;
            if(j == 1), ylabel(char(chnm(chnsR(i))),'Fontsize',8); end
            if(i == 1), title(char(chnm(chnsL(j))),'Fontsize',8); end
        
%             set(0, 'CurrentFigure', Gamma)
            subplot(length(chnsR)*2,length(chnsL),(i-1)*(length(chnsL))+j+12)
            %             set(gca,'XTick',1:div); set(gca,'XTickLabel',num2cell(Gf));
            hold on; plot(G,line); text(1:div,G,num2cell(GPhi),'k','Fontsize',7); hold off;
            if(j == 1), ylabel(char(chnm(chnsR(i))),'Fontsize',8); end
            if(i == 1), title(char(chnm(chnsL(j))),'Fontsize',8); end
            
        end
    end
end

set(0, 'CurrentFigure', Beta); ax = axes; t1 = title('Beta = Top, Gamma = Bot_ '); ax.Visible = 'off'; t1.Visible = 'on';
set(0, 'CurrentFigure', Gamma); ax = axes; t1 = title('Gamma'); ax.Visible = 'off'; t1.Visible = 'on';
set(Beta,'visible','on'); set(Gamma,'visible','on');

%% Plot coherence points at maximum amplitude of beta/gamma across hemispheres over time
chnsR = [5,6]; 
chnsL = [7,8,13,14,15,16]+17; 
params.tapers = [3,5]; params.Fs = fs; params.fpass = [9,100]; params.trialave = 1;
movingwin = [0.25,0.01];

beta = [14,30];
gamma = [80,100];

Beta = figure; set(Beta,'visible','off'); 
% Gamma = figure; set(Gamma,'visible','off'); 

for t = 1:3
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
    
    inds = floor(-window*(1/back)*fs/1000:1:window*fs/1000);
    trialinds = repmat((trig')./dwn, length(inds), 1) + repmat(inds(:), 1, size(trig,1));
    trialinds(:,floor(trialinds(1,:))<=0) = []; trialinds(:,floor(trialinds(end,:))>length(data)) = [];
    
    for i = 1:length(chnsR)
        d = data(:,chnsR(i));%filtfilt(b,a,double(data(:,chnsR(i)))); %comb notch
        Snips = u.meanSubtract(d(floor(trialinds)));
        for j = 1:length(chnsL)
            d2 = data(:,chnsL(j));%filtfilt(b,a,double(data(:,chnsL(j)))); %comb notch
            Snips2 = u.meanSubtract(d2(floor(trialinds)));
            [C,phi,S12,S1,S2,t,f] = cohgramc(Snips,Snips2,movingwin,params);
            
            B = median(C(:,f>=beta(1)&f<=beta(2)),2); for k = 1:length(B), Bind(k) = find(C(k,:)==B(k));end; 
            G = median(C(:,f>=gamma(1)&f<=gamma(2)),2); for k = 1:length(G), Gind(k) = find(C(k,:)==G(k));end; 
%             [B,Bind] = max(C(:,f>=beta(1)&f<=beta(2)),[],2); 
%             [G,Gind] = max(C(:,f>=gamma(1)&f<=gamma(2)),[],2);  Gind = Gind+find(f>=gamma(1),1)-1;
            
            set(0, 'CurrentFigure', Beta);
            subplot(length(chnsR)*2,length(chnsL),(i-1)*(length(chnsL))+j);
%             hold on; yyaxis right; plot(t,f(Bind)); ylim([beta(1)-10,beta(2)+10]); hold off;
            hold on; yyaxis left; plot(t,B,line,'linewidth',1.5); hold off;% text(1:div,B,num2cell(BPhi),'Fontsize',7); hold off;
            if(j == 1), ylabel(char(chnm(chnsR(i))),'Fontsize',8); end
            if(i == 1), title(char(chnm(chnsL(j))),'Fontsize',8); end
            
%             set(0, 'CurrentFigure', Gamma)
            subplot(length(chnsR)*2,length(chnsL),(i-1)*(length(chnsL))+j+12)
%             hold on; yyaxis right; plot(t,f(Gind)); ylim([gamma(1)-10,gamma(2)+10]); hold off;
            hold on; yyaxis left; plot(t,G,line,'linewidth',1.5); hold off; %text(1:div,G,num2cell(GPhi),'k','Fontsize',7); hold off;
            if(j == 1), ylabel(char(chnm(chnsR(i))),'Fontsize',8); end
            if(i == 1), title(char(chnm(chnsL(j))),'Fontsize',8); end
            
        end
    end
end

set(0, 'CurrentFigure', Beta); ax = axes; t1 = title('Beta = Top, Gamma = Bottom_ '); ax.Visible = 'off'; t1.Visible = 'on';
% set(0, 'CurrentFigure', Gamma); ax = axes; t1 = title('Gamma'); ax.Visible = 'off'; t1.Visible = 'on';
set(Beta,'visible','on'); %set(Gamma,'visible','on');


%% Granger Causality - pairwise for all left to right channels and vice versa
% define params
chnsR = [11,12];%[5,6]; 
chnsL = [5,6,13,14,15,16]+17; 

params.ntrials   = length(trig);     % number of trials
params.nobs      = length(inds);   % number of observations per trial

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

R2L = figure; set(R2L,'visible','off'); 
L2R = figure; set(L2R,'visible','off'); 

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
    
    inds = floor(-window*(1/back)*fs/1000:1:window*fs/1000);
    trialinds = repmat((trig')./dwn, length(inds), 1) + repmat(inds(:), 1, size(trig,1));
    trialinds(:,floor(trialinds(1,:))<=0) = []; trialinds(:,floor(trialinds(end,:))>length(data)) = [];
    
    div = 1; % number of divisions over time
    for i = 1:length(chnsR)
        d = Filter(:,chnsR(i));%filtfilt(b,a,double(data(:,chnsR(i)))); %comb notch
        Snips = u.meanSubtract(d(floor(trialinds)));
        %     Snips = Snips(:,randperm(size(Snips,2))); % random permutation in trial order
        X(1,:,:) = Snips;
        for j = 1:length(chnsL)
            d2 = Filter(:,chnsL(j));%filtfilt(b,a,double(data(:,chnsL(j)))); %comb notch
            Snips2 = u.meanSubtract(d2(floor(trialinds)));
            %         Snips2 = Snips2(:,randperm(size(Snips2,2))); % random permutation in trial order
            X(2,:,:) = Snips2;
            
            for k = 1:div
                %             switch k
                %                 case 1
                %                     line = 'k';
                %                 case 2
                %                     line = 'b:';
                %                 case 3
                %                     line = 'r--';
                %                 case 4
                %                     line = 'g-.';
                %             end
                Lbound = floor((k-1)/div*length(X)+1);
                Rbound = floor(k/div*length(X));
                
                disp([num2str(i),num2str(j)]);
                
                [gc,f] = a.GC(X(:,Lbound:Rbound,:),params); % custom helper function for using MVGC toolbox
                
                % plot
                set(0, 'CurrentFigure', R2L)
                subplot(length(chnsR),length(chnsL),(i-1)*(length(chnsL))+j)
                hold on; plot(f,gc(1,:),line,'Linewidth',1.5); xlim(params.frange); hold off;
                if(j == 1), ylabel(char(chnm(chnsR(i))),'Fontsize',8); end
                if(i == 1), title(char(chnm(chnsL(j))),'Fontsize',8); end
                
                set(0, 'CurrentFigure', L2R)
                subplot(length(chnsR),length(chnsL),(i-1)*(length(chnsL))+j)
                hold on; plot(f,gc(2,:),line,'Linewidth',1.5); xlim(params.frange); hold off;
                if(j == 1), ylabel(char(chnm(chnsR(i))),'Fontsize',8); end
                if(i == 1), title(char(chnm(chnsL(j))),'Fontsize',8); end
            end
        end
    end
end
set(0, 'CurrentFigure', R2L); ax = axes; t1 = title('R2L'); ax.Visible = 'off'; t1.Visible = 'on';
set(0, 'CurrentFigure', L2R); ax = axes; t1 = title('L2R'); ax.Visible = 'off'; t1.Visible = 'on';
set(R2L,'visible','on'); set(L2R,'visible','on');



%% Granger Causality - with permutations removed?
% Ahh, do permutations for the time domain GC to see if the p-values given
% (whether or not there is causality essentially) is significant. Basically
% getting a p-value of the p-values.



%% Granger Causality Over TIME. Basically stealing the method used for spectra in chronux for better comparisons
movingwin = [0.25,0.01];
win = round(movingwin(1)*fs); wstep = round(movingwin(2)*fs);

chnsR = [5,6]; 
chnsL = [7,8,13,14,15,16]+17; 

params.regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
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

R2L = figure; set(R2L,'visible','off'); 
L2R = figure; set(L2R,'visible','off'); 

beta = [14,30]; gamma = [80,100];

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

    inds = floor(-window*(1/back)*fs/1000:1:window*fs/1000);
    trialinds = repmat((trig')./dwn, length(inds), 1) + repmat(inds(:), 1, size(trig,1));
    trialinds(:,floor(trialinds(1,:))<=0) = []; trialinds(:,floor(trialinds(end,:))>length(data)) = [];
    
    params.ntrials   = length(trig);     % number of trials
    params.nobs      = length(inds);   % number of observations per trial

    for i = 1:length(chnsR)
        d = Filter(:,chnsR(i));%filtfilt(b,a,double(data(:,chnsR(i)))); %comb notch
        Snips = u.meanSubtract(d(floor(trialinds)));
        
        X(1,:,:) = Snips;
        for j = 1:length(chnsL)
            d2 = Filter(:,chnsL(j));%filtfilt(b,a,double(data(:,chnsL(j)))); %comb notch
            Snips2 = u.meanSubtract(d2(floor(trialinds)));
            
            X(2,:,:) = Snips2;
            disp([num2str(i),num2str(j)]);
            winstart = 1:wstep:size(X,2)-win+1; winmid=winstart+round(win/2);
            
            gcb = zeros(length(winstart),2); fb = zeros(length(winstart),2);
            gcg = zeros(length(winstart),2); fg = zeros(length(winstart),2);
            
            for k = 1:length(winstart)
                disp([num2str(round(((i-1)/length(chnsR)+ (1/length(chnsR))*(j-1)/length(chnsL) + (1/(length(chnsR)*length(chnsL))) * k/length(winstart))*100,1)) ,'%'])
                Lbound = winstart(k); Rbound = Lbound+win-1;
                [gc,f] = a.GC(X(:,Lbound:Rbound,:),params);
                B = median(gc(:,f>=beta(1)&f<=beta(2))'); 
                G = median(gc(:,f>=gamma(1)&f<=gamma(2))'); 
%                 [B,Bind] = max(gc(:,f>=beta(1)&f<=beta(2))'); F1 = f(Bind+find(f>=beta(1),1)-1);
%                 [G,Gind] = max(gc(:,f>=gamma(1)&f<=gamma(2))'); F2 = f(Gind+find(f>=gamma(1),1)-1);
                gcb(k,:) = B; fb(k,:) = F1;  gcg(k,:) = G; fg(k,:) = F2;
            end
            
            set(0, 'CurrentFigure', R2L)
            subplot(length(chnsR)*2,length(chnsL),(i-1)*(length(chnsL))+j)
            hold on; plot(winmid/fs,gcb(:,1),line,'Linewidth',1.5);   hold off;
%             if(i == 1 && j == 1), legend('Beta','Gamma'); end;
            if(j == 1), ylabel(char(chnm(chnsR(i))),'Fontsize',8); end
            if(i == 1), title(char(chnm(chnsL(j))),'Fontsize',8); end
            subplot(length(chnsR)*2,length(chnsL),(i-1)*(length(chnsL))+j+length(chnsR)*length(chnsL))
%             hold on; yyaxis left; plot(winmid/fs,fb(:,1),line);
            %         yyaxis right; plot(winmid/fs,fg(:,1),'r');  hold off;
            hold on; plot(winmid/fs,gcg(:,1),line,'Linewidth',1.5); hold off;
            if(j == 1), ylabel(char(chnm(chnsR(i))),'Fontsize',8); end
            if(i == 1), title(char(chnm(chnsL(j))),'Fontsize',8); end
            
            set(0, 'CurrentFigure', L2R)
            subplot(length(chnsR)*2,length(chnsL),(i-1)*(length(chnsL))+j)
            hold on; plot(winmid/fs,gcb(:,2),line,'Linewidth',1.5); hold off;
%             if(i == 1 && j == 1), legend('Beta','Gamma'); end;
            if(j == 1), ylabel(char(chnm(chnsR(i))),'Fontsize',8); end
            if(i == 1), title(char(chnm(chnsL(j))),'Fontsize',8); end
            subplot(length(chnsR)*2,length(chnsL),(i-1)*(length(chnsL))+j+length(chnsR)*length(chnsL))
%             hold on; yyaxis left; plot(winmid/fs,fb(:,2),line);
            %         yyaxis right; plot(winmid/fs,fg(:,2),'r');  hold off;
            hold on; plot(winmid/fs,gcg(:,2),line,'Linewidth',1.5);  hold off;
            if(j == 1), ylabel(char(chnm(chnsR(i))),'Fontsize',8); end
            if(i == 1), title(char(chnm(chnsL(j))),'Fontsize',8); end
            
        end
    end
end

set(0, 'CurrentFigure', R2L); ax = axes; t1 = title('R2L, Top:Beta Bot:Gamma_ '); ax.Visible = 'off'; t1.Visible = 'on';
set(0, 'CurrentFigure', L2R); ax = axes; t1 = title('L2R, Top:Beta Bot:Gamma_ '); ax.Visible = 'off'; t1.Visible = 'on';
set(R2L,'visible','on'); set(L2R,'visible','on');



%% Granger Causality Spectra
movingwin = [0.25,0.01];
win = round(movingwin(1)*fs); wstep = round(movingwin(2)*fs);

chnsR = 5;%[5,6]; 
chnsL = [24,30,32]%[7,8,13,14,15,16]+17; 

params.regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
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

R2L = figure; set(R2L,'visible','off'); 
L2R = figure; set(L2R,'visible','off'); 

beta = [14,30]; gamma = [80,100];
clim = [0,0];
for t = 1:1
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

    inds = floor(-window*(1/back)*fs/1000:1:window*fs/1000);
    trialinds = repmat((trig')./dwn, length(inds), 1) + repmat(inds(:), 1, size(trig,1));
    trialinds(:,floor(trialinds(1,:))<=0) = []; trialinds(:,floor(trialinds(end,:))>length(data)) = [];
    
    params.ntrials   = length(trig);     % number of trials
    params.nobs      = length(inds);   % number of observations per trial

    for i = 1:length(chnsR)
        d = Filter(:,chnsR(i));%filtfilt(b,a,double(data(:,chnsR(i)))); %comb notch
        Snips = u.meanSubtract(d(floor(trialinds)));
        
        X(1,:,:) = Snips;
        for j = 1:length(chnsL)
            d2 = Filter(:,chnsL(j));%filtfilt(b,a,double(data(:,chnsL(j)))); %comb notch
            Snips2 = u.meanSubtract(d2(floor(trialinds)));
            
            X(2,:,:) = Snips2;
            disp([num2str(i),num2str(j)]);
            winstart = 1:wstep:size(X,2)-win+1; winmid=winstart+round(win/2);
            
            len = 500;
            gcb = zeros(length(winstart),len); fb = zeros(length(winstart),2);
            gcg = zeros(length(winstart),len); fg = zeros(length(winstart),2);
            
            for k = 1:length(winstart)
                disp([num2str(round(((i-1)/length(chnsR)+ (1/length(chnsR))*(j-1)/length(chnsL) + (1/(length(chnsR)*length(chnsL))) * k/length(winstart))*100,1)) ,'%'])
                Lbound = winstart(k); Rbound = Lbound+win-1;
                [gc,f] = a.GC(X(:,Lbound:Rbound,:),params);
                gcb(k,:) = u.Interp2Length(gc(1,:),len);
                gcg(k,:) = u.Interp2Length(gc(2,:),len);
            end
            
            set(0, 'CurrentFigure', R2L)
            subplot(length(chnsR),length(chnsL),(i-1)*(length(chnsL))+j)
            imagesc(winmid/fs,linspace(params.frange(1),params.frange(2),len),gcb'); set(gca,'YDir','normal');  
            if(diff(caxis)>diff(clim)), clim = caxis; end
            if(j == 1), ylabel(char(chnm(chnsR(i))),'Fontsize',8); end
            if(i == 1), title(char(chnm(chnsL(j))),'Fontsize',8); end
            
            set(0, 'CurrentFigure', L2R)
            subplot(length(chnsR),length(chnsL),(i-1)*(length(chnsL))+j)
            imagesc(winmid/fs,linspace(params.frange(1),params.frange(2),len),gcg'); set(gca,'YDir','normal');
            if(diff(caxis)>diff(clim)), clim = caxis; end
            if(j == 1), ylabel(char(chnm(chnsR(i))),'Fontsize',8); end
            if(i == 1), title(char(chnm(chnsL(j))),'Fontsize',8); end
            
        end
    end
end
for i = 1:(length(chnsR)*length(chnsL))
    set(0, 'CurrentFigure', R2L); subplot(length(chnsR),length(chnsL),i); caxis(clim);
    set(0, 'CurrentFigure', L2R); subplot(length(chnsR),length(chnsL),i); caxis(clim);
end
set(0, 'CurrentFigure', R2L); ax = axes; t1 = title('R2L, Top:Beta Bot:Gamma_ '); ax.Visible = 'off'; t1.Visible = 'on';
set(0, 'CurrentFigure', L2R); ax = axes; t1 = title('L2R, Top:Beta Bot:Gamma_ '); ax.Visible = 'off'; t1.Visible = 'on';
set(R2L,'visible','on'); set(L2R,'visible','on');


%% Granger Causality Trial by Trial (every 10 trials since GC calculation fails with too few trials)
% doesn't seem to show anything significant? 
movingwin = [0.25,0.01];
win = round(movingwin(1)*fs); wstep = round(movingwin(2)*fs);

chnsR = 5;%[5,6]; 
chnsL = [24,30,32];%[7,8,13,14,15,16]+17; 

params.regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
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

R2L = figure; set(R2L,'visible','off'); 
L2R = figure; set(L2R,'visible','off'); 

beta = [14,30]; gamma = [80,100];

trig = trigPreL;
rtsl = SL.rts_l((lefttrials(:,1)-15<trig1(1)) & lefttrialsuccess);

inds = floor(-window*(1/back)*fs/1000:1:window*fs/1000);
trialinds = repmat((trig')./dwn, length(inds), 1) + repmat(inds(:), 1, size(trig,1));
bad = floor(trialinds(1,:))<=0 | floor(trialinds(end,:))>length(data);
trialinds(:,bad) = []; rtsl(bad) = []; [rtsl,rti] = sort(rtsl);

params.ntrials   = length(trig);     % number of trials
params.nobs      = length(inds);   % number of observations per trial

X = [];
for i = 1:length(chnsR)
    d = Filter(:,chnsR(i));%filtfilt(b,a,double(data(:,chnsR(i)))); %comb notch
    Snips = u.meanSubtract(d(floor(trialinds)));
    
    X(1,:,:) = Snips;
    for j = 1:length(chnsL)
        d2 = Filter(:,chnsL(j));%filtfilt(b,a,double(data(:,chnsL(j)))); %comb notch
        Snips2 = u.meanSubtract(d2(floor(trialinds)));
        
        X(2,:,:) = Snips2;
        X = X(:,:,rti);
        
        RT = []; r2l = []; l2r = [];
        for k = 1:floor(size(X,3)/10)
            [gc,f] = a.GC(X(:,:,(k-1)*10+1:k*10),params);
            RT(k) = nanmedian(rtsl((k-1)*10+1:k*10));
            B = median(gc(:,f>=beta(1)&f<=beta(2))'); 
            G = median(gc(:,f>=gamma(1)&f<=gamma(2))'); 
            r2l(k,:) = [B(1),G(1)]; l2r(k,:) = [B(2),G(2)];
        end
                
        [RT,I] = sort(RT); r2l = r2l(I,:); l2r = l2r(I,:);
        
        set(0, 'CurrentFigure', R2L)
        subplot(length(chnsR)*2,length(chnsL),(i-1)*(length(chnsL))+j)
        plot(RT,r2l(:,1))
        subplot(length(chnsR)*2,length(chnsL),(i-1)*(length(chnsL))+j+length(chnsR)*length(chnsL))
        plot(RT,r2l(:,2))
        if(j == 1), ylabel(char(chnm(chnsR(i))),'Fontsize',8); end
        if(i == 1), title(char(chnm(chnsL(j))),'Fontsize',8); end
        
        set(0, 'CurrentFigure', L2R)
        subplot(length(chnsR)*2,length(chnsL),(i-1)*(length(chnsL))+j)
        plot(RT,l2r(:,1))
        subplot(length(chnsR)*2,length(chnsL),(i-1)*(length(chnsL))+j+length(chnsR)*length(chnsL))
        plot(RT,l2r(:,2))
        if(j == 1), ylabel(char(chnm(chnsR(i))),'Fontsize',8); end
        if(i == 1), title(char(chnm(chnsL(j))),'Fontsize',8); end
    end
end

set(0, 'CurrentFigure', R2L); ax = axes; t1 = title('R2L, Top:Beta Bot:Gamma_ '); ax.Visible = 'off'; t1.Visible = 'on';
set(0, 'CurrentFigure', L2R); ax = axes; t1 = title('L2R, Top:Beta Bot:Gamma_ '); ax.Visible = 'off'; t1.Visible = 'on';
set(R2L,'visible','on'); set(L2R,'visible','on');



%% Test
test = [];
for i = 1:floor(length(SL.rts_l)/10)
    test(i) = nanmean(SL.rts_l((i-1)*10+1:i*10));
end


