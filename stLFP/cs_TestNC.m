clear; close all;

ncPath = 'Y:\~NeuroWest\Spanky\Neurochip\';
session = 'S_20180613_07\';
file = 'S_20180613_07.mat';
file_name = [ncPath,session,file];

% events = nc3events(file_name);
% 
% [data, names, session_time] = nc3data(17:32, 0, 60, 20000, [], file_name);

load(file_name);

bw = [p.event_lower_bandwidth(7),p.event_upper_bandwidth(7)];
thresh = p.event_threshold(7);
win1del = p.event_window1_delay(7); %in ms
win1max = p.event_window1_max(7); %in uV
win1min = p.event_window1_min(7);
win2del = p.event_window2_delay(7);
win2max = p.event_window2_max(7);
win2min = p.event_window2_min(7);

session = 'S_20180605_01\';
file = 'S_20180605_01.mat';
file_name = [ncPath,session,file];

sr = 20000;

[data, names, session_time] = nc3data(17:32, 0, 60, sr, [], file_name);

trigChn = 11; stimChn = 9;

trigData = bpfilt(data(:,trigChn),bw,sr,3);
trigData = -trigData;

crossing = trigData>thresh;
crossing = find(diff(crossing));

win1del = round(win1del/1000*sr);
win2del = round(win2del/1000*sr);

win1 = trigData(crossing+win1del) >= win1min & trigData(crossing+win1del) <= win1max; 
win2 = trigData(crossing+win2del) >= win2min & trigData(crossing+win2del) <= win2max; 

spks = crossing(win1&win2)';

window = 0.001;
range = round(-window*sr:1:1.5*window*sr);
trialinds = repmat(spks, length(range), 1) + repmat(range(:), 1, size(spks,2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(trigData)) = [];

trigData = -trigData;

snips = trigData(trialinds);

figure; plot(mean(snips,2))




