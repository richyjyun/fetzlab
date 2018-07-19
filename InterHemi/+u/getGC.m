function SLGC = getGC(SL)

SLGC = struct([]);

for i = 1:length(SL)
    if (~isempty(SL(i).Bad) || strcmp(SL(i).Condition, 'Control')...
            || strcmp(SL(i).Condition, 'NaN') || strcmp(SL(i).Condition, 'nostim')...
            || isempty(SL(i).trig1) || strcmp(SL(i).Condition(end),'R'))
        continue;
    end
    
    D = str2num(SL(i).Date);
    Bad = 20120515;%[20120508,20120515,20120612,20120814];
    if( any(Bad==D))
        continue;
    end
    
    disp(SL(i).Date);
    
    %% load neural data
    [data, fs, chnm] = loadData(SL(i));
    if(isempty(data)), continue; end
    
    %% Left trials
    Trials = SL(i).lefttrials; RT = SL(i).rts_l; Epochs = []; Stim = [];
    bounds = [find(Trials(:,2)<SL(i).trig1(1),1,'last'),find(Trials(:,1)>SL(i).trig1(end),1)];
    % If triggered by lefttrials for kato/igor, no conditioning epoch
    if (~strcmp(SL(i).Animal,'Ubi') || (D>=20170226 && D<=20170306))&& ...
            ((strcmp(SL(i).StimHemi,'L') && strcmp(SL(i).Condition(1:4),'Ipsi')) ||...
            (strcmp(SL(i).StimHemi,'R') && strcmp(SL(i).Condition(1:6),'Contra')))
        Epochs = [1,bounds(1);bounds(2),length(Trials)]; 
    else
        %%%%%%%%% NEED TO REMOVE TRIALS WHEN STIM OCCURED 
        if(strcmp(SL(i).Animal,'Ubi') && strcmp(SL(i).Condition(1:4),'Ipsi'))
            Stim = getStimTrials(SL(i).trig1,Trials);
        end
        Epochs = [1,bounds(1);bounds(1)+1,bounds(2)-1;bounds(2),length(Trials)];
    end     
    SLtemp = calculateGC(Trials,RT,Epochs,Stim,SL(i),data,fs,chnm,'L');
    
    %% Right Trials
    Trials = SL(i).righttrials; RT = SL(i).rts_r; Epochs = []; Stim = [];
    bounds = [find(Trials(:,2)<SL(i).trig1(1),1,'last'),find(Trials(:,1)>SL(i).trig1(end),1)];
    % If triggered by righttrials for kato/igor, no conditioning epoch
    if (~strcmp(SL(i).Animal,'Ubi') || (D>=20170226 && D<=20170306))&& ...
            ((strcmp(SL(i).StimHemi,'R') && strcmp(SL(i).Condition(1:4),'Ipsi')) ||...
            (strcmp(SL(i).StimHemi,'L') && strcmp(SL(i).Condition(1:6),'Contra')))
        Epochs = [1,bounds(1);bounds(2),length(Trials)]; 
    else
        %%%%%%%%% NEED TO REMOVE TRIALS WHEN STIM OCCURED 
        if(strcmp(SL(i).Animal,'Ubi') && strcmp(SL(i).Condition(1:6),'Contra'))
             Stim = getStimTrials(SL(i).trig1,Trials);
        end
        Epochs = [1,bounds(1);bounds(1)+1,bounds(2)-1;bounds(2),length(Trials)];
    end     
    SLtemp = calculateGC(Trials,RT,Epochs,Stim,SLtemp,data,fs,chnm,'R');
    
    %% Set up new struct
    SLnew = struct;
    SLnew.Date = SL(i).Date;
    SLnew.fs = fs;
    SLnew.tGC = SLtemp.tGC;
    SLnew.chnm =  SLtemp.chnm;
    SLnew.chID = SLtemp.chID;
    SLnew.LGC = SLtemp.LGC;
    SLnew.RGC = SLtemp.RGC;
    
    if(isempty(SLGC))
        SLGC = SLnew;
    else
        SLGC(end+1) = SLnew;
    end
end
end

function SL = calculateGC(Trials,RT,Epochs,Stim,SL,data,fs,chnm,side)
window = 1; chID = []; inds = (-window*fs:1:window*fs);
Beta = [];
Gamma = [];

% get viable channels
for chn = 1:length(chnm)
    nm = chnm{chn};
    if(strcmp(SL.Animal,'Ubi') && ~(strcmp(nm(1),'R') || strcmp(nm(1),'L')))
        continue;
    end
    chID(end+1) = chn;
end

for j = 1:size(Epochs,1)
    
    trial = Epochs(j,1):Epochs(j,2); trial = trial(~ismember(trial,Stim));
    trig = Trials(trial,1) + RT(trial);
    trig(isnan(trig)) = []; trig = trig*fs/1000;
    trialinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
    trialinds(:,trialinds(1,:)<=0) = []; trialinds(:,trialinds(end,:) > length(data)) = [];
           
    % parameters for GC
    params.tapers = [3,5]; params.Fs = fs;
    params.regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default) tried LWR, takes much longer but gives same results
    params.icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
    
    params.morder    = round(fs*0.05); % 50ms, from PNAS paper
    params.momax     = params.morder*3;% maximum model order for model order estimation. Default takes way too long, gives same results. Apparently often the same as model order, but 3* to give it more
    
    params.acmaxlags = [];   % maximum autocovariance lags (empty for automatic calculation)
    
    params.statt     = 'F';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
    params.alpha_sig     = 0.005;   % significance level for significance test
    params.mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
    
    params.fs        = fs;    % sample rate (Hz)
    params.fres      = [];     % frequency resolution (empty for automatic calculation)
    
    params.nvars     = 2;     % number of variables. only doing pairwise
    
    params.frange    = [9,100]; %range of frequencies to look at
    
    params.ntrials   = length(trig);     % number of trials
    params.nobs      = length(inds);   % number of observations per trial
   
    movingwin = [0.25,0.01];
    win = round(movingwin(1)*fs); wstep = round(movingwin(2)*fs);
    
    % CHANNELS FOR UBI
    goodChannels = [3,4,12,24,25,31];
    
    NSnips = zeros(length(goodChannels),length(inds),size(trialinds,2)); badtrials = [];
    for chn = 1:length(goodChannels)
        d = data(:,goodChannels(chn));
        NSnips(chn,:,:) = u.meanSubtract(d(floor(trialinds)),params);
        badtrials = [badtrials;find(any(NSnips(chn,:,:) > mean(std(NSnips(chn,:,:)))*10))]; %remove ones with outliers (noise, artifact, etc)
    end
    NSnips(:,:,badtrials) = [];
    
    % Do GC in windows
    winstart = 1:wstep:size(NSnips,2)-win+1; winmid=winstart+round(win/2);
    startadd = nan; GC = [];
    for k = 1:length(winstart)
        if(isnan(startadd))
            Lbound = winstart(k); Rbound = Lbound+win-1;
        else
            Lbound = winstart(startadd); Rbound = winstart(k)+win-1;
        end
        x = NSnips(:,Lbound:Rbound,:);
        [gc,p,s] = a.GC_time(x,params);
        GC(j,k,:,:) = gc;
        % if causality diverges, add more to the time span
        if(isnan(gc))
            startadd = k;
        else
            startadd = nan;
        end
    end 
end

SL.tGC = winmid;
SL.([side,'GC']) = GC;
SL.chnm = chnm;
SL.chID = goodChannels;

end


function Stim = getStimTrials(trig,Trials)
Stim = nan(1,length(trig));
for j = 1:length(trig)
    norm = Trials(:,1)-((trig(j))+50);
    ind = find(norm>0,1)-1;
    
    if isempty(ind)
        continue;
    end
    
    if(abs(Trials(ind,2)-trig(j)) > abs(Trials(ind+1,1)-trig(j)))
        ind = ind+1;
    end
    Stim(j) = ind;
end
Stim(isnan(Stim)) = []; 
end



function [data, fs, chnm] = loadData(SL)
% Loading in the neural data. No need to load in .f32 since everything
    % should already be set to the aligned data
    if(strcmp(SL.Animal,'Ubi'))
        D = SL.Date;
        S = SL.Session_Guger_Train;
        Session = [char(D),'_',char(S(2))];
        gugfile = Session;
        if(~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
            data = []; fs = []; chnm = []; return;
        end
        % down sampling rate
        dwn = 20;
        [data, fs, chnm, ~] = u.LoadGug(gugfile, dwn); % modify so I'm only loading in one channel to make it go faster
        
    elseif(strcmp(SL.Animal,'Kato'))
        D = SL.Date;
        Sessions = strsplit(char(SL.Sessions),'_');
        dwn = 10;
        data = [];
        for j = 1:3
            gugfile = [char(D),'_',char(Sessions(j))];
            if(~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
                 data = []; fs = []; chnm = []; return;
            end
            [d, fs, chnm, ~] = u.LoadGug(gugfile, dwn);
            data = [data;d]; 
        end
        
    elseif(strcmp(SL.Animal,'Igor'))
        D = SL.Date;
        Sessions = strsplit(char(SL.Sessions),'_');
        data = [];
        for j = 1:3
            gugfile = [char(D),'_',char(Sessions(j))];
            if( ~exist([gugfile,'.bin']) || ~exist([gugfile,'.cfg']))
                 data = []; fs = []; chnm = []; return;
            end
            [d, fs, chnm] = u.ImportIgorNeuralData(gugfile);
            
            data = [data;d];
        end

    end
end




