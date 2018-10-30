% clear; close all;

%% Load in data
chn = 32; sc = 1;

tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';
blknm = 'Spanky-180122-130528';

T1 = times(1,1);
T2 = times(1,2);
LFPs1 = TDT2mat([tankpath,blknm],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','verbose',false);
LFPs1 = LFPs1.streams.LFPs;

% get all snippets from spike sorting
Snips1 = TDT2mat([tankpath,blknm],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu','Channel',chn,'verbose',false);
Snips1 = Snips1.snips.eNe1;


T1 = times(2,1);
T2 = times(2,2);
LFPs2 = TDT2mat([tankpath,blknm],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','verbose',false);
LFPs2 = LFPs2.streams.LFPs;

% get all snippets from spike sorting
Snips2 = TDT2mat([tankpath,blknm],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu','Channel',chn,'verbose',false);
Snips2 = Snips2.snips.eNe1;

fs = LFPs1.fs;

%% phase preference
stimChn = 64;

noise = findNoise(LFPs1.data(stimChn,:),fs);
filt1 = LFPs1.data(stimChn,:); filt1(noise) = 0;
filt1 = bpfilt(filt1,[15,25],fs,3); 

h1 = hilbert(filt1); h1(noise) =nan;
amp1 = abs(h1); ang1 = angle(h1);

spk1 = Snips1.ts(Snips1.sortcode==sc)-times(1,1); spk1 = round(spk1*fs);

phase1 = ang1(spk1);

noise = findNoise(LFPs2.data(stimChn,:),fs);
filt2 = LFPs2.data(stimChn,:); filt2(noise) = 0;
filt2 = bpfilt(filt2,[15,25],fs,3); 

h2 = hilbert(filt2); h2(noise) = nan;
amp2 = abs(h2); ang2 = angle(h2);

spk2 = Snips2.ts(Snips2.sortcode==sc)-(times(2,1)); spk2 = round(spk2*fs);

phase2 = ang2(spk2);

bins = -pi:pi/16:pi;

figure; histogram(phase1,bins); hold on; histogram(phase2,bins);

%% difference in the beta oscillations
noise = findNoise(LFPs1.data(chn,:),fs);
filt1 = LFPs1.data(chn,:); filt1(noise) = 0;
filt1 = bpfilt(filt1,[15,25],fs,3);

h = hilbert(filt1); h(noise) = nan; ang = angle(h); amp = abs(h);

d1 = circDiff(ang,ang1); 

amp = amp(spk1);
b = discretize(ang(spk1),bins);
a = [];
for i = 1:(length(bins)-1)
    a(i) = nanmean(amp(b==i));
end

noise = findNoise(LFPs2.data(chn,:),fs);
filt2 = LFPs2.data(chn,:); filt2(noise) = 0;
filt2 = bpfilt(filt2,[15,25],fs,3);

h = hilbert(filt2); h(noise) = nan; ang = angle(h);

d2 = circDiff(ang,ang2);

bins = 0:pi/32:pi;
figure; histogram(d1,bins,'Normalization','probability'); hold on; histogram(d2,bins,'Normalization','probability');

%% difference in phase distribution during spiking

d1spk = d1(spk1);
d2spk = d2(spk2);

figure; histogram(d1spk,bins,'Normalization','probability'); hold on; histogram(d2spk,bins,'Normalization','probability');

figure; histogram(d1(trough1),bins,'Normalization','probability'); hold on; histogram(d2(trough2),bins,'Normalization','probability');

%% phase distribution in stim channel at time of trough in trig channel
noise = findNoise(LFPs1.data(chn,:),fs);
filt1 = LFPs1.data(chn,:); filt1(noise) = 0;
filt1 = bpfilt(filt1,[15,25],fs,3);

h = hilbert(filt1); ang = angle(h); amp = abs(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% This is where i messed up!!!!
[~,trough1] = findpeaks(ang);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trough1(amp(trough1)<std(filt1)) = [];

noise = findNoise(LFPs2.data(chn,:),fs);
filt2 = LFPs2.data(chn,:); filt2(noise) = 0;
filt2 = bpfilt(filt2,[15,25],fs,3);

h = hilbert(filt2); ang = angle(h); amp = abs(h);

[~,trough2] = findpeaks(ang);
trough2(amp(trough2)<std(filt2)) = [];

p1 = ang1(trough1);

p2 = ang2(trough2);

bins = -pi:pi/16:pi;
figure; histogram(p1,bins,'Normalization','probability'); hold on; histogram(p2,bins,'Normalization','probability');

%% 2 window discrim
% get cycle trigger
bp = [15,25];
stim = LFPs2.data(stimChn,:);
noise = findNoise(stim,fs);
stim(noise) = 0; stim = bpfilt(stim,bp,fs,3);
stim(noise) = nan;

thresh = nanstd(stim); % set threshold to be a standard deviation

% find threshold crossings
trig = stim<=-thresh; trig = diff(trig); trig(trig<0) = 0;
trig = find(trig);

% window discrimination
period = 1/mean(bp);
delay1 = round((period/5)*fs); % 10ms delay
delay2 = round((period/5+period/2)*fs); % 25ms delay after delay1
trig((trig+delay2)>length(stim)) = [];
trig = trig((stim(trig+delay1) <= -thresh) & (stim(trig+delay2) >= thresh));

noise = findNoise(LFPs1.data(stimChn,:),fs);
filt1 = LFPs1.data(stimChn,:); filt1(noise) = 0;
filt1 = bpfilt(filt1,[15,25],fs,3); 

h1 = hilbert(filt1); h1(noise) =nan;
amp1 = abs(h1); ang1 = angle(h1);

noise = findNoise(LFPs2.data(stimChn,:),fs);
filt1 = LFPs2.data(stimChn,:); filt1(noise) = 0;
filt1 = bpfilt(filt1,[15,25],fs,3); 

h1 = hilbert(filt1); h1(noise) =nan;
amp1 = abs(h1); ang2 = angle(h1);

figure; histogram(ang1(trig1),bins,'Normalization','probability'); hold on; histogram(ang2(trig2),bins,'Normalization','probability');


%% plot triggered LFPs
window = [-0.05,0.05];
range = round(window(1)*fs:1:window(2)*fs);

trig = trough1;

trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(LFPs1.data)) = [];

noise = findNoise(LFPs1.data(chn,:),fs);
ctLFP1 = LFPs1.data(stimChn,:); 
ctLFP1(noise) = 0;
ctLFP1 = bpfilt(ctLFP1,[15,25],fs,3);

h1 = hilbert(ctLFP1); amp1 = abs(h1); phase1 = angle(h1);

ctLFP1(noise) = nan;

ctLFP1 = ctLFP1(trialinds);
ctAmp1 = amp1(trialinds); ctAng1 = phase1(trialinds);
% figure; hist3([amp1(trig)',phase1(trig)'],[100,32])

trig = trough2;

trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(LFPs2.data)) = [];

noise = findNoise(LFPs2.data(chn,:),fs);
ctLFP2 = LFPs2.data(stimChn,:); 
ctLFP2(noise) = 0;
ctLFP2 = bpfilt(ctLFP2,[15,25],fs,3);

h2 = hilbert(ctLFP2); amp2 = abs(h2); phase2 = angle(h2);
% figure; hist3([amp2(trig)',phase2(trig)'],[100,32])

ctLFP2(noise) = nan;

ctLFP2 = ctLFP2(trialinds); ctLFP2 = zscore(ctLFP2);
ctAmp2 = amp2(trialinds); ctAng2 = phase2(trialinds);

figure; plot(range/fs,zscore(nanmean(ctLFP1,2))); hold on; plot(range/fs,zscore(nanmean(ctLFP2,2)));
figure; plot(range/fs,(nanmean(ctAmp1,2))); hold on; plot(range/fs,(nanmean(ctAmp2,2)));
figure; plot(range/fs,(nanmean(ctAng1,2))); hold on; plot(range/fs,(nanmean(ctAng2,2)));



% spike trig
trig = spk1';

trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(LFPs1.data)) = [];

noise = findNoise(LFPs1.data(chn,:),fs);
stLFP1 = LFPs1.data(stimChn,:); stLFP1(noise) = 0;
stLFP1 = bpfilt(stLFP1,[15,25],fs,3);

h1 = hilbert(stLFP1); amp1 = abs(h1); phase1 = angle(h1);

stLFP1(noise) = nan;

stLFP1 = stLFP1(trialinds); stLFP1 = stLFP1-nanmean(stLFP1);
stAmp1 = amp1(trialinds); stAng1 = phase1(trialinds);
% figure; hist3([amp1(spk1)',phase1(spk1)'],[100,32])

x = 1:size(stAmp1,1)'; y = nanmean(stAmp1,2)'; dy = nanstd(stAmp1,1,2)';
figure;
% fill([x;flipud(x)],[y-dy;flipud(y+dy)],[.5 .5 .5],'linestyle','none');
% plot(x,y); hold on; plot(x,y-dy,'k--'); plot(x,y+dy,'k--');


trig = spk2';

trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(LFPs2.data)) = [];

noise = findNoise(LFPs2.data(chn,:),fs);
stLFP2 = LFPs2.data(stimChn,:); stLFP2(noise) = 0;
stLFP2 = bpfilt(stLFP2,[15,25],fs,3);

h2 = hilbert(stLFP2); amp2 = abs(h2); phase2 = angle(h2);

stLFP2(noise) = nan;

stLFP2 = stLFP2(trialinds); stLFP2 = stLFP2-nanmean(stLFP2);
stAmp2 = amp2(trialinds); stAng2 = phase2(trialinds);
% figure; hist3([amp2(spk2)',phase2(spk2)'],[100,32])


figure; plot(range/fs,zscore(nanmean(stLFP1,2))); hold on; plot(range/fs,zscore(nanmean(stLFP2,2)));
figure; plot(range/fs,(nanmean(stAmp1,2))); hold on; plot(range/fs,(nanmean(stAmp2,2)));
figure; plot(range/fs,(nanmean(stAng1,2))); hold on; plot(range/fs,(nanmean(stAng2,2)));



%% cross correlation
figure; CrossCorr(spk1/fs, 'ts2',spk1/fs,'binsize', 0.001,'lag',[window(1),window(2)]*2,'suppress_plot',0); 
figure; CrossCorr(spk2/fs, 'ts2',spk2/fs,'binsize', 0.001,'lag',[window(1),window(2)]*2,'suppress_plot',0); 

figure; CrossCorr(spk1/fs,'ts2',spk2/fs,'binsize', 0.001,'lag',[window(1),window(2)]*2,'suppress_plot',0); 


%% find closest differences
closestVal = zeros(1,length(spk1));
closestInd = zeros(1,length(spk1));
for i = 1:length(spk1)
    [closestVal(i),closestInd(i)] = min(abs(spk1(i)-spk2));
end

temp = closestInd(closestVal<10);

inds = find(Snips1.sortcode==1);
figure; plot(mean(Snips2.data(inds(closestVal<10),:)))

inds = find(Snips2.sortcode == 2);
figure; plot(mean(Snips2.data(inds(temp),:)))









