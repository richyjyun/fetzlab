clear;

tankpath = 'Y:\~NeuroWest\Spanky\IFNN\';
% tankpath = 'Y:\~NeuroWest\Spanky\FncGenTest\';
% tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';
blockname = 'Spanky-180713-113402';

% chn = 78;
chn = 78;
T1 = 0; T2 = 5*60;
LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','Channel',chn);
LFPs = LFPs.streams.LFPs; fs = LFPs.fs;

% h = hilbert(LFPs.data);

%% define noise
thresh = 1.5e-4;

noise = abs(h)>thresh; d = diff(noise);
nInds = [find(d==1)',find(d==-1)'];

i = 1;
while i < length(nInds)
    if(nInds(i+1,1)-nInds(i,2) < fs)
        nInds(i,2) = nInds(i+1,2);
        nInds(i+1,:) = [];
    else
        i = i+1;
    end
end

% events = diff(nInds');

%% plot
t = (1:length(LFPs.data))/fs;
figure; plot(t,LFPs.data); hold on; plot(t,abs(h));

for i = 1:length(nInds)
    plot(nInds(i,:)/fs,[thresh,thresh]*2,'k','linewidth',2);
end

%% get power spectrum of noise
F = {};
P = {};

maxval = 0; maxind = 0;
% figure; subplot(2,1,1);
for i = 1:length(nInds)
    if(nInds(i,2)-nInds(i,1) < 8)
        continue;
    end
    [P{i},F{i}] = pwelch(LFPs.data(nInds(i,1):nInds(i,2)),[],[],[],fs);
%     plot(F{i},P{i}); hold on;
    if(length(F{i}) > maxval)
        maxind = i;
        maxval = length(F{i});
    end
end

power = [];
for i = 1:length(F)
    if(isempty(F{i}))
        continue;
    end
    power(:,end+1) = interp1(F{i},P{i},F{maxind});
end

figure; yyaxis left; plot(F{maxind},mean(power,2));

%% define good
gInds = [1,nInds(1,1)-1];
for i = 2:length(nInds)-1
    gInds(i,:) = [nInds(i-1,2)+1,nInds(i,1)-1];
end
gInds(end+1,:) = [nInds(end,2)+1,length(LFPs.data)];

F = {};
P = {};

maxval = 0; maxind = 0;
% subplot(2,1,2);
for i = 1:length(gInds)
    if(gInds(i,2)-gInds(i,1) < 8)
        continue;
    end
    [P{i},F{i}] = pwelch(LFPs.data(gInds(i,1):gInds(i,2)),[],[],[],fs);
%     plot(F{i},P{i}); hold on;
    if(length(F{i}) > maxval)
        maxind = i;
        maxval = length(F{i});
    end
end

power = [];
for i = 1:length(F)
    if(isempty(F{i}))
        continue;
    end
    power(:,end+1) = interp1(F{i},P{i},F{maxind});
end

yyaxis right; plot(F{maxind},mean(power,2));

%% how much of the recording is noise?
noiselength = sum(diff(nInds'));
noiseratio = noiselength/length(LFPs.data);

%% using chronux for moving window spectra

params.tapers = [3,5]; params.Fs = fs; params.fpass = [9,100]; params.trialave = 1;
movingwin = [5,1];

% Moving window spectrum
figure; subplot(3,1,1); plot(LFPs.data); xlim([1,length(LFPs.data)]); title('Raw LFP'); axis off;
[S,t,f] = mtspecgramc(LFPs.data,movingwin,params);
subplot(3,1,[2,3]); imagesc(t,f,log(S')); title('Moving Window Spectra')
ylabel('Freq (Hz)'); xlabel('Time (s)');

% [S,f] = mtspectrumc(temp,params);
% figure; plot(f,S);

%% coherence
params.tapers = [3,5]; params.Fs = fs; params.fpass = [5,100]; params.trialave = 1;
movingwin = [0.25,0.01];

% chn = 81;
chn = 81;
LFPs2 = TDT2mat([tankpath,blockname],'TYPE',4,'STORE','LFPs','Channel',chn);
LFPs2 = LFPs2.streams.LFPs; fs = LFPs2.fs;

start = 1; window = 10000;

C = []; S1 = []; S2 = [];
while start+window < length(LFPs.data)
    [C(end+1,:),phi,S12,S1(end+1,:),S2(end+1,:),f]=coherencyc(LFPs.data(start:start+window-1),LFPs2.data(start:start+window-1),params);
    start = start+window;
end

figure;
subplot(3,1,1); plot(f,mean(S1)); title('Chn1 Spectrum');
subplot(3,1,2); plot(f,mean(S2)); title('Chn2 Spectrum');
subplot(3,1,3); plot(f,mean(C)); title('Coherence');
xlabel('Frequency (Hz)');

% figure; plot(f,C);




