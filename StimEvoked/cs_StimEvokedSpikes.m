%% Set variables
clear; close all;
path = 'Y:\~NeuroWest\Spanky\Neurochip';
day = 'S_20181029_14';

%% Load meta data
[fpath,fname,Channels,fs,session_time] = getNCData(path,day);

%% Load relevant data
trigChn = 7; event = 1; bw = [300,3000];
[data, names, ~] = nc3data(trigChn, 0, session_time, fs, bw, fname);

events = nc3events(fname);
stim = events.stim{1};
load(fname);

win1min = -p.event_window1_max(event);
win1max = -p.event_window1_min(event);
win2min = -p.event_window2_max(event);
win2max = -p.event_window2_min(event);

%% Get when stim evoked potential
window = [0.0015,0.003];
range = round(window(1)*fs):1:round(window(2)*fs);

trig = round(stim*fs)';
trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(data)) = [];

trials = data(trialinds);

norm = median(trials,2);

trials = trials-norm;

evoked = any(trials>win2min & trials<win2max) & any(trials<win1max & trials>win1min); %may have to adjust on a per day basis
evokedts = stim(evoked);
[~,conddelay] = min(trials(:,evoked)); 

%% Bin over time
binwidth = 60;
bins = 0:binwidth:session_time;
[stimn,~] = histcounts(stim,bins);
[evokedn,~] = histcounts(evokedts,bins);

figure; subplot(2,1,1); bar(bins(2:end)-binwidth/2,evokedn./stimn,1,'k');
ylabel('% Evoked Spikes');
yyaxis right; ax = gca; ax.YColor = 'r';
plot(bins(2:end)-binwidth/2,stimn/binwidth,'r');
ylabel('Firing Rate (Hz)');
xlabel(['Seconds; Bins: ',num2str(binwidth),'s']);
xlim([0,bins(end)]); title('Evoked Spikes & Firing Rate');
xt = xticks;

idx = discretize(evokedts,bins);
delaydist = [];
bw = 2;
b = 1:bw:(max(conddelay)+1);
for i = 1:(length(bins)-1)
    delays = conddelay(idx==i);
    [delaydist(i,:),~] = histcounts(delays,b);
end

subplot(2,1,2); imagesc((delaydist./sum(delaydist,2))'); set(gca,'Ydir','normal');
yt = [10,20,30]/bw; ytl = ((yt*bw)/fs+window(1))*1000;
yticks(yt); yticklabels(ytl);
xticks(xt/binwidth); xticklabels(xt);
ylabel(['Evoked Delay (ms), Bins: ',num2str(bw),'ms']);
xlabel(['Seconds; Bins: ',num2str(binwidth),'s']);
title('Normalized Evoked Delay Distribution');
c = colorbar('southoutside');
ylabel(c,'% Distribution');

% printFile = ['F:\S\Packets\NC\',day,'_EvokedSpikes.ps'];
% print('-painters','-fillpage',gcf, '-dpsc2', printFile, '-append');
% callps2pdf(printFile);




%% Bin by number of stimuli
binwidth = 500; 
evokedshaped = [evoked,zeros(1,binwidth*ceil(length(stim)/binwidth)-length(stim))];
evokedshaped = reshape(evokedshaped,binwidth,[]);

figure; bar((0:binwidth:(binwidth*(size(evokedshaped,2)-1)))+binwidth/2,sum(evokedshaped)/binwidth,1,'k');
ylabel('% Evoked Spikes'); xlabel(['N Stimuli, Bins: ',num2str(binwidth),' Stim']);
xlim([0,binwidth*size(evokedshaped,2)])

%% Cross correlation of cell
window = 0.1;
idx = discretize(stim,bins);
cor = [];
for i = 1:(length(bins)-1)
   [cor(:,i),lags] = CrossCorr(stim(idx==i), 'ts2',stim(idx==i),'binsize', 0.001,'lag',[-window,window]); 
end

cor = cor./sum(cor);

figure; bar(lags,cor(:,1),1,'facealpha',0.5); 
hold on; bar(lags,cor(:,end),1,'r','facealpha',0.5);

figure;
initc = [0.8,0.8,0.8];
for i = 1:size(cor,2)
    plot(lags,cor(:,i),'Color',initc*(1-i/size(cor,2)));
    hold on;
end










