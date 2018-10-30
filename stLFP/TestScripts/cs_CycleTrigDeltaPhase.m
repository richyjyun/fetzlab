clear; close all;

%% Load in proper times
% tankpath = 'Y:\~NeuroWest\Spanky\RandomStim-180314-124242\';
tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';
% tankpath = 'Y:\~NeuroWest\Spanky\IFNN\';
% tankpath = 'Y:\~NeuroWest\Spanky\CycleTriggered-170710-143939\';
% tankpath = 'Y:\~NeuroWest\Spanky\Connectivity-180207-131758\';
% tankpath = 'Y:\~NeuroWest\Spanky\SpankyUpstairs\';
blockname = 'Spanky-180122-130528';
TT = TDT2mat([tankpath,blockname],'TYPE',2);
Dscm = TT.epocs.Dscm;
[val,ind] = findpeaks(Dscm.data); 
val(end+1) = Dscm.data(end); ind(end+1) = length(Dscm.data);
ind = ind(val>1000); val = val(val>1000); 
tests = 1:2:length(val); 
times = [ind(tests)-val(tests),ind(tests)];%-val(tests)+val(2)];  % times = [ind(1) - val(1),ind(1)-val(1)+15000;ind(3)-val(3),ind(3)-val(3)+10000];
times = Dscm.onset(times);

%% Load in the two channels
T1 = times(1,1); T2 = times(1,2);

LFP32 = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','Channel',32); LFP32 = LFP32.streams.LFPs;
LFP64 = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','Channel',64); LFP64 = LFP64.streams.LFPs;

fs = LFP32.fs;

% get cycles 
bp = [15,25];
stim = bpfilt(LFP32.data',bp,fs,3)'; %filter one channel
testtime = round(2*160*fs); % 2 minutes of test
thresh = std(stim(1:testtime)); % set threshold to be a standard deviation

% find threshold crossings
trig = stim<=-thresh; trig = diff(trig); trig(trig<0) = 0;
trig = find(trig);

% window discrimination
period = 1/mean(bp);
delay1 = round((period/5)*fs); % 10ms delay
delay2 = round((period/5+period/2)*fs); % 25ms delay after delay1
trig((trig+delay2)>length(stim)) = [];
trig = trig((stim(trig+delay1) <= -thresh) & (stim(trig+delay2) >= thresh));

% determine trough time
range = 0:round(0.05*fs);
trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(LFP32.data)) = [];

temp = stim(trialinds);

h = hilbert(stim);
h = angle(h);
amp = abs(h);
h = h(trialinds);

d = diff(h); 
[~,ind] = max(abs(d));

trough = trig+ind;
mag = amp(trough);

% filter and hilbert 64 for comparison
filt = bpfilt(LFP64.data',bp,fs,3)'; %filter one channel

temp = filt(trialinds);

h = hilbert(filt);
h = angle(h);
% h = h(trialinds);
% 
% d = diff(h); 
% [~,ind2] = max(abs(d));
% 
% delta = ind2-ind;

phases = h(trough);

phases = phases-pi;
phases = wrapToPi(phases);

bins = -pi:pi/32:pi;
figure; histogram(phases,bins,'normalization','probability')


%% after stim
%% Load in the two channels
T1 = times(2,1); T2 = times(2,2);

LFP32 = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','Channel',32); LFP32 = LFP32.streams.LFPs;
LFP64 = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','Channel',64); LFP64 = LFP64.streams.LFPs;

fs = LFP32.fs;

% get cycles 
bp = [15,25];
stim = bpfilt(LFP32.data',bp,fs,3)'; %filter one channel
testtime = round(2*160*fs); % 2 minutes of test
thresh = std(stim(1:testtime)); % set threshold to be a standard deviation

% find threshold crossings
trig = stim<=-thresh; trig = diff(trig); trig(trig<0) = 0;
trig = find(trig);

% window discrimination
period = 1/mean(bp);
delay1 = round((period/5)*fs); % 10ms delay
delay2 = round((period/5+period/2)*fs); % 25ms delay after delay1
trig((trig+delay2)>length(stim)) = [];
trig = trig((stim(trig+delay1) <= -thresh) & (stim(trig+delay2) >= thresh));

% determine trough time
range = 0:round(0.05*fs);
trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(LFP32.data)) = [];

temp = stim(trialinds);

h = hilbert(stim);
h = angle(h);
amp = abs(h);
h = h(trialinds);

d = diff(h); 
[~,ind] = max(abs(d));

trough = trig+ind;
mag2 = amp(trough);

d = diff(h); 
[~,ind] = max(abs(d));

trough = trig+ind;

% filter and hilbert 64 for comparison
filt = bpfilt(LFP64.data',bp,fs,3)'; %filter one channel
h = hilbert(filt);
h = angle(h);
phases2 = h(trough);

phases2 = phases2-pi;
phases2 = wrapToPi(phases2);

hold on; histogram(phases2,bins,'normalization','probability')


xbins = -pi:pi/16:pi;
ybins = linspace(min([mag,mag2]),max([mag,mag2]),32);

figure; h1 = histogram2(phases,mag,xbins,ybins,'FaceColor','flat','normalization','probability'); 
% view(2); colorbar; h1.BinCounts(h1.BinCounts==0) = 1;

figure; h2 = histogram2(phases2,mag2,xbins,ybins,'FaceColor','flat','normalization','probability'); 
% view(2); colorbar; h2.BinCounts(h2.BinCounts==0) = 1;
% cl(2,:) = caxis;

c1 = h1.Values; x1 = h1.XBinEdges; y1 = h1.YBinEdges;
figure; imagesc(c1'); set(gca,'YDir','normal')
xticks(linspace(1,size(c1,1),5));
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
yticks(linspace(1,size(c1,2),3))
yticklabels(ybins([1,floor(length(ybins)/2),length(ybins)]));
colorbar; cl(1,:) = caxis;

c2 = h2.Values; x2 = h2.XBinEdges; y2 = h2.YBinEdges;
figure; imagesc(c2'); set(gca,'YDir','normal')
xticks(linspace(1,size(c2,1),5));
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
yticks(linspace(1,size(c2,2),3))
yticklabels(ybins([1,floor(length(ybins)/2),length(ybins)]));
colorbar; cl(2,:) = caxis;

caxis([0,max(cl(:,2))]);





