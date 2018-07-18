clear; close all;

%% Load
% 2. Epochs - Thd1 (threshold + eventgen1 values), Dscm (discrimination)
% 3. Snips - Beta (snippets, 24kHz)
% 4. Streams - Mani (manipulandum, 3kHz), LFPs (3kHz), SUAs (24kHz), Filt (3kHz)
% 5. Scalars (all at 1Hz) - Trig (trigger channel), SUAc (SUA channels), Stim (Stim params and times) 
tankpath = 'Y:\~NeuroWest\Spanky\CycleTriggered-170710-143939\';
blockname = 'Spanky-170920-134709';
T1 = 60; T2 = 70;  % in seconds. 0 to denote start or end of entire recording
TT = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2);
% Mani = TDT2mat([tankpath,blockname],'TYPE',4,'STORE','Mani');
trigChn = TT.scalars.Trig.data(1);


%% Filter data exactly as the TDT system filters it
% Biquad filter with 2nd order Butterworth. Need to do low and high pass separately instead of a bandpass. 
beta = [12,25];
data = TT.streams.LFPs.data(trigChn,:);
fs = TT.streams.LFPs.fs;
[z,p,k] = butter(2,beta(2)/(fs/2),'low');  % low pass filter
sos = zp2sos(z,p,k);
Filter = sosfilt(sos,data);
[z,p,k] = butter(2,beta(1)/(fs/2),'high'); % high pass filter
sos = zp2sos(z,p,k);
Filter = double(sosfilt(sos,Filter));


%% Plot comparison between online and offline filtered data
if(isfield(TT.streams,'Filt'))
    figure; subplot(2,1,1)
    plot(TT.streams.Filt.data); hold on; plot(Filter);
    
    params.tapers = [3,5]; params.Fs = fs; params.fpass = [5,30]; params.trialave = 1;
    subplot(2,1,2)
    [S,f] = mtspectrumc(TT.streams.Filt.data,params); plot(f,S);
    [S,f] = mtspectrumc(Filter,params); hold on; plot(f,S);
    
    figure; plot(f,S);
    [S,f] = mtspectrumc(data,params); hold on; plot(f,S);
     
    [r,lags] = xcorr(Filter,TT.streams.Filt.data,floor(fs/2),'coeff');
    figure; plot(lags/fs,r);
    offset = lags(r==max(r))/fs;
end;

%% Characterizing the beta. Have resting period to determine threshold? Use boxID? Manipulandum data?
Hilb = abs(hilbert(Filter)); figure; plot((1:length(Hilb))/fs,Hilb); % instantaneous power
window = 0.05; % 50 ms intervals to envelope over
window = round(window*fs);
% Env = envelope(Hilb,window,'peak'); hold on; plot((1:length(Hilb))/fs,Env);
M = mean(Hilb); S = std(Hilb); 
hold on; plot([1,length(Hilb)]/fs,[M,M],'r'); 
hold on; plot([1,length(Hilb)]/fs,[M+S,M+S],'r--');
S = S;
figure; plot((1:length(Hilb))/fs,Filter); hold on; plot([1,length(Hilb)]/fs,[3*S,3*S],'r--'); 


%% Find Offset of generated sine waves. Basically delay of just the filter
% the delay seems to be very frequency dependent...
if(isfield(TT.streams,'SinI'))
    in = TT.streams.SinI.data; out = TT.streams.SinO.data; Filt = sosfilt(sos,in);
    fs = TT.streams.SinI.fs;
    figure; plot(in); hold on; plot(out); hold on; plot(Filt);
    [r,lags] = xcorr(out,in,floor(0.02*fs),'coeff');
    figure; plot(lags/fs,r); inds = lags>=0; lags = lags(inds);
    offset = lags(r(inds)==max(r(inds)))/fs;
    
    params.tapers = [3,5]; params.Fs = fs; params.fpass = [5,30]; params.trialave = 1;
    [S,f] = mtspectrumc(in,params); figure; plot(f,S);
    [S,f] = mtspectrumc(out,params); hold on; plot(f,S);
    [S,f] = mtspectrumc(Filt,params); hold on; plot(f,S);
end

%% Snippet data
% Snippets used in the win discrim script. Don't really need this for analysis,
% as we can extract snippets easily. 
snips = TT.snips.Beta.data;
figure;
for i = 1:size(snips,1)
    hold on; plot(snips(i,:))
end

%% Find Offset to original raw data compared to online filtered data
% This is the offset that will be used for phase calculations
data = TT.streams.LFPs.data(trigChn,:);
fs = TT.streams.LFPs.fs;
figure; plot(Filter); hold on; plot(data);
[r,lags] = xcorr(Filter,data,floor(fs/2),'coeff');
figure; plot(lags/fs,r); inds = lags>=0; lags = lags(inds);
offset = lags( r(inds)==max(r(inds)))/fs;

%% Plot filter with when stimulation occured
thr = TT.epocs.Thdb.data(1); tic = [0,1];
data = Filter;
figure; plot(data); hold on; plot([0,length(data)],[thr,thr],'k');
stim = TT.epocs.Thd1.onset-T1;
stim = stim(logical(TT.epocs.Thd1.data));
window2 = TT.epocs.Thdf;

yyaxis right; ofs = 0;%0.0375*fs; %76 % ah, epoch data is triggering at the second window limit!
for i = 1:length(stim)
    hold on; plot([stim(i)*fs-ofs,stim(i)*fs-ofs],tic,'r')
end

%% Fitting a sine curve to find the frequency at each stim
ofs = offset; %s NEED TO FIGURE OUT OFFSET TO USE FOR PHASE CALCULATION
stim = TT.epocs.Epo1.onset-T1;
stim = stim(logical(TT.epocs.Epo1.data));
freq = nan(length(stim),1);
phase = nan(length(stim),1);
options.Algorithm = 'levenberg-marquardt';
options.MaxIterations = 1000; options.Display = 'off';
window = 0.05; %window to look before stimulus
for i = 1:length(freq)
    start = round((stim(i)-window)*fs); finish = round(stim(i)*fs); 
    if(~(start > 0 && finish <= length(Filter)))
        continue;
    end
    snip = double(Filter(start:finish)); snip = snip-mean(snip);
    t = (1:length(snip))/fs;
    sinusoid = @(x) x(1)*sin(((2*pi)/x(2))*t+x(3)) +x(4) - snip;
    x = [max(snip),.05,0,0]; 
    X = lsqnonlin(sinusoid,double(x),[],[],options);
    %     figure; plot(t,snip); hold on; plot(t,X(1)*sin(((2*pi)/X(2))*t+X(3)));
    freq(i) = 1/X(2);
    fitted = X(1)*sin(((2*pi)/X(2))*t+X(3));
    [~,p0] = findpeaks(diff(fitted)); 
    if(~isempty(p0))
        phase(i) = ((length(fitted)-p0(end)+ofs*fs)/fs)/abs(X(2))*360;
    else
        [~,p0] = findpeaks(-diff(fitted));
        if(isempty(p0))
            continue;
        end
        phase(i) = ((length(fitted)-p0(end))/fs)/abs(X(2))*360+180;
    end
end
figure; subplot(2,2,1); histogram(freq); xlabel('Frequency of Stim Wave');
hold on; plot([nanmedian(freq),nanmedian(freq)],ylim,'r','linewidth',1.5);
title([num2str(nanmedian(freq)),'Hz']);
subplot(2,2,3); boxplot(freq);
freq2 = freq(freq<25 & freq>10);
subplot(2,2,2); histogram(freq2); xlabel('Frequency (10-25Hz)');
hold on; plot([nanmedian(freq2),nanmedian(freq2)],ylim,'r','linewidth',1.5);
title([num2str(nanmedian(freq2)),'Hz']);
subplot(2,2,4); boxplot(freq2);
figure; histogram(phase); 
yyaxis right; x = linspace(1,360,1000); plot(x,sin(2*pi/360.*x),'linewidth',1.5); 
hold on; plot([nanmedian(phase),nanmedian(phase)],ylim,'r','linewidth',1.5);
title([num2str(nanmedian(phase)),' degrees']);
xlabel('Phase of Stim (degrees)');


%% Manipulanum data
Manip = TT.streams.Mani.data(1:2,:)'; plot((1:length(Manip))/fs,Manip(:,1));
movingavg = movmean(Manip(:,1),150);
hold on; plot((1:length(Manip))/fs,movingavg);
yyaxis right; plot((1:length(Manip)-1)/fs,diff(movingavg));
figure; plot(diff(Manip));

%% Check for CCEPs
%% Load data
Scalars = TDT2mat([tankpath,blockname],'TYPE',5);
LFPs = TDT2mat([tankpath,blockname],'TYPE',4,'STORE','LFPs');
Stim = Scalars.scalars.Stim; LFPs = LFPs.streams.LFPs; 
%% Run
fs = LFPs.fs; T1 = 0;
current = Stim.data(5,:); inds = [1,find(diff(current)~=0)+1];
t = Stim.ts - T1; amplitudes = nan(length(inds)-1,1);
window = 0.05; range = round(0:1:window*fs); %50 ms window for capturing CCEPs
params.tapers = [3,5]; params.Fs = fs;
for i = 1:length(inds)-1
    figure;
    trig = round(t(inds(i):inds(i+1)-1)*fs);
    trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,1));
    trialinds(:,floor(trialinds(1,:))<=0) = []; 
    trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];
    for j = 1:size(LFPs.data,1)
        [r,c,e] = GetWadeChannelPosition(j);
        subplot(10,10,(r-1)*10+c);
        d = LFPs.data(j,:);
%         if(size(trialinds,2) == 1)
%             d = u.meanSubtract(d(trialinds),params);
%         else
        d = mean(d(trialinds),2);
%         end
%         subplot(12,8,j); hold on; 
        plot(range/fs,d); ylim([-40*1e-6,40*1e-6]);
        axis off
        sub_pos = get(gca,'position'); % get subplot axis position
        set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
        title(num2str(j),'fontsize',7)
    end
    amplitudes(i) = current(inds(i));
end
amplitudes(end+1) = current(inds(end));








