function SLNeuro = getLFPTime(SL)

SLNeuro = struct([]);

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
        Epochs = [1,bounds(1);bounds(1)+1,bounds(2)-1;bounds(2),length(Trials)];
    else
        if(str2num(SL(i).Date) >= 20170226  || str2num(SL(i).Date) <= 20170306)
            Stim = [];
        elseif(strcmp(SL(i).Animal,'Ubi') && strcmp(SL(i).Condition(1:4),'Ipsi'))
            Stim = getStimTrials(SL(i).trig1,Trials);
        end
        Epochs = [1,bounds(1);bounds(1)+1,bounds(2)-1;bounds(2),length(Trials)];
    end
    SLtemp = getOscillations(Trials,RT,Epochs,Stim,SL(i),data,fs,chnm,'L');
    
    %% Right Trials
    Trials = SL(i).righttrials; RT = SL(i).rts_r; Epochs = []; Stim = [];
    bounds = [find(Trials(:,2)<SL(i).trig1(1),1,'last'),find(Trials(:,1)>SL(i).trig1(end),1)];
    % If triggered by righttrials for kato/igor, no conditioning epoch
    if (~strcmp(SL(i).Animal,'Ubi') || (D>=20170226 && D<=20170306))&& ...
            ((strcmp(SL(i).StimHemi,'R') && strcmp(SL(i).Condition(1:4),'Ipsi')) ||...
            (strcmp(SL(i).StimHemi,'L') && strcmp(SL(i).Condition(1:6),'Contra')))
        Epochs = [1,bounds(1);bounds(1)+1,bounds(2)-1;bounds(2),length(Trials)];
    else
        if(str2num(SL(i).Date) >= 20170226  || str2num(SL(i).Date) <= 20170306)
            Stim = [];
        elseif(strcmp(SL(i).Animal,'Ubi') && strcmp(SL(i).Condition(1:6),'Contra'))
            Stim = getStimTrials(SL(i).trig1,Trials);
        end
        Epochs = [1,bounds(1);bounds(1)+1,bounds(2)-1;bounds(2),length(Trials)];
    end
    SLtemp = getOscillations(Trials,RT,Epochs,Stim,SLtemp,data,fs,chnm,'R');
    
    %% Set up new struct
    SLnew = struct;
    SLnew.Date = SL(i).Date;
    SLnew.fs = fs;
    SLnew.tneuro = SLtemp.tneuro;
    SLnew.chnm =  SLtemp.chnm;
    SLnew.chID = SLtemp.chID;
    SLnew.Lbeta = SLtemp.Lbeta;
    SLnew.Lgamma = SLtemp.Lgamma;
    SLnew.Rbeta = SLtemp.Rbeta;
    SLnew.Rgamma = SLtemp.Rgamma;
%     SLnew.LRP = SLtemp.LRP;
%     SLnew.RRP = SLtemp.RRP;
    
    if(isempty(SLNeuro))
        SLNeuro = SLnew;
    else
        SLNeuro(end+1) = SLnew;
    end
end
end

function SL = getOscillations(Trials,RT,Epochs,Stim,SL,data,fs,chnm,side)
window = 2; chID = []; inds = (-window*fs:1:window*fs);
Bind = {0,0.5,1}; Bind = cellfun(@(x) find(inds>=(x*fs),1),Bind); %[find(inds>(0.3*fs),1),find(inds<(0.7*fs),1,'last')];
Gind = {-0.5,0,0.5}; Gind = cellfun(@(x) find(inds>=(x*fs),1),Gind); %find(inds>(-0.1*fs),1),find(inds<(0.1*fs),1,'last')];
Beta = {};
Gamma = {};
% RP = [];

% get viable channels
for chn = 1:length(chnm)
    nm = chnm{chn};
    if(strcmp(SL.Animal,'Ubi') && ~(strcmp(nm(1),'R') || strcmp(nm(1),'L')))
        continue;
    elseif(strcmp(SL.Animal,'Kato') && ~(strcmp(nm(1),'E') || strcmp(nm(1),'I')))
        continue;
    end
    chID(end+1) = chn;
end

for j = 1:size(Epochs,1)
    
    trial = Epochs(j,1):Epochs(j,2); trial = trial(~ismember(trial,Stim));
    trig = Trials(trial,1) + RT(trial);
    triglength = length(trig);
    trig = trig*fs/1000;
    if(length(trig) == 0)
        continue;
    end
    trialinds = repmat(trig', length(inds), 1) + repmat(inds(:), 1, size(trig,1));
    naninds = isnan(trig) | (trialinds(1,:)<=0)' | (trialinds(end,:) > length(data))';
    trialinds(:,trialinds(1,:)<=0) = []; trialinds(:,trialinds(end,:) > length(data)) = [];
    trialinds(:,isnan(trialinds(1,:))) = [];

    
    params.tapers = [3,5]; params.Fs = fs; params.trialave = 1;
    %     rmv = [];
    
    % loop through viable channels
    for chn = 1:length(chID)
        
        d = data(:,chID(chn));
        NSnips = u.meanSubtract(d(floor(trialinds)),params);
%         if(strcmp(SL.Animal,'Ubi'))
%             NSnips(:,any(NSnips > mean(std(NSnips))*10)) = []; %remove ones with outliers (noise, artifact, etc)
%         end
        
        %         Testing
        %         chnm(chID(chn))
        %         figure; subplot(2,1,1); plot(NSnips);
        %         subplot(2,1,2); [pxx,fr] = periodogram(NSnips,[],[],fs); plot(fr(fr>10&fr<150),mean(pxx(fr>10&fr<150,:),2));
        
        beta = [15,25]; gamma = [45,120];
        if(strcmp(SL.Animal,'Igor'))
            beta = [15,21]; gamma = [80,110];
        elseif strcmp(SL.Animal,'Ubi')
            beta = [19,25]; gamma = [85,100];
        elseif strcmp(SL.Animal,'Kato')
            beta = [15,25]; gamma = [80,105];
        end
        
        %         % Testing
        %         bp = gamma;
        %         movingwin = [0.25,0.01]; params.fpass = bp;
        %         [S,t,freq] = mtspecgramc(NSnips,movingwin,params);
        %         figure; subplot(2,2,1); plot(t-window,S);
        %         subplot(2,2,2); plot(t-window,mean(S')); yyaxis right; plot(inds/1000,mean(NSnips,2))
        %
        %         bandpass = bp;
        %         Wn_theta = [bandpass(1)/(fs/2) bandpass(2)/(fs/2)]; % normalized by the nyquist frequency
        %         [btheta,atheta] = butter(3,Wn_theta);
        %         sig = filtfilt(btheta,atheta,NSnips);
        %         power = abs(hilbert(sig));
        %         subplot(2,2,3); plot(inds,smooth(median(power,2),20)); yyaxis right; plot(inds,mean(NSnips,2))
        %         rms = 5; % low pass cutoff hz
        %         [brms,arms] = butter(1,rms/(params.Fs/2),'low'); % lowpass (see daq_sapi_*)
        %         power = max(filtfilt(brms,arms,abs(sig)),zeros(size(sig))); % filter with lowpass
        %         power = median(power,2);
        %         subplot(2,2,4); plot(inds,smooth(power,50)); yyaxis right; plot(inds,mean(NSnips,2))
        %
        windows = [-0.3,0; 0,0.3; 0.3,0.6];
        windows = windows*fs + window*fs;

        Power = getPower(NSnips,beta,params,windows);
        if(isempty(Power))
            Beta{j,chn} = nan;
        else
            temp = nan(length(Bind),triglength);
            temp(:,~naninds) = Power;
            Beta{j,chn} = Power;
        end
        
        Power = getPower(NSnips,gamma,params,windows);
        if(isempty(Power))
            Gamma{j,chn} = nan;
        else
            temp = nan(length(Gind),triglength);
            temp(:,~naninds) = Power;
            Gamma{j,chn} = Power;
        end
        
%         RP(j,chn,:) = mean(NSnips,2);
        
    end
    
end

% Beta(:,rmv,:) = [];
% Gamma(:,rmv,:) = [];

SL.tneuro = inds/fs;
SL.([side,'beta']) = Beta;
SL.([side,'gamma']) = Gamma;
% SL.([side,'RP']) = RP;
SL.chnm = chnm;
SL.chID = chID;

end

function power = getPower(NSnips,bp,params,windows)

% Richardson method
Wn_theta = [bp(1)/(params.Fs/2) bp(2)/(params.Fs/2)]; % normalized by the nyquist frequency
[btheta,atheta] = butter(3,Wn_theta);
sig = filtfilt(btheta,atheta,NSnips);

% Richardson method
rms = 5; % low pass cutoff hz
[brms,arms] = butter(1,rms/(params.Fs/2),'low'); % lowpass (see daq_sapi_*)
power = max(filtfilt(brms,arms,abs(sig)),zeros(size(sig))); % filter with lowpass

if(isempty(power))
    return;
end

power = [nanmean(power(windows(1,1):windows(1,2),:));...
    nanmean(power(windows(2,1):windows(2,2),:));...
    nanmean(power(windows(3,1):windows(3,2),:))];

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
    dwn = 10;
    [data, fs, chnm, ~] = u.LoadGug(gugfile, dwn); % modify so I'm only loading in one channel to make it go faster
    
elseif(strcmp(SL.Animal,'Kato'))
    D = SL.Date;
    Sessions = strsplit(char(SL.Sessions),'_');
    dwn = 5;
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







