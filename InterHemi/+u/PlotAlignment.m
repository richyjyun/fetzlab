%% Alignment Figures
% Showing that aligning data off of (re)calculated reaction times gives
% tighter alignment of both accelerometer(movement) and guger(neural) data

% After loading in data (20170221) without downsampling
% 
% SL.accel_raw_r = accel(:,2);
% SL.accel_raw_l = accel(:,1);
% SL.fs = 9600;
% SL.lefttrials = lefttrials;
% SL.righttrials = righttrials;
% SL.lefttrialsuccess = lefttrialsuccess;
% SL.righttrialsuccess = righttrialsuccess; 
% SL.trig1 = trig1; 
% 
% SL = a.AppendReactionTimes2(SL); % for calculating reaction time with upsampled data

% Define neural data and accelerometer data
d = Filter(:,rchn+1);

%% Aligning off of trial onset
fs = 960; dwn = 10; trig1 = SL.trig1;
trig = trig1;
inds = floor(0:1:window*fs*dwn/1000);
trialinds = repmat((trig'), length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds(:,floor(trialinds(1,:))<=0) = []; trialinds(:,floor(trialinds(end,:))>length(SL.accel_raw_r)) = [];

figure; plot(SL.accel_raw_r(floor(trialinds)))

% accel data
figure;
Snip = accel(:,2); Snip = u.meanSubtract(Snip(floor(trialinds(:,170:250)))); Snip = abs(Snip); hold on; plot(Snip); plot(mean(Snip'),'r','LineWidth',2);
plot(median(Snip'),'k','LineWidth',2); hold off;

% neural data
figure;
trialinds(:,floor(trialinds(1,:)/dwn)<=0) = []; trialinds(:,floor(trialinds(end,:)/dwn)>length(data)) = [];
Snip = u.meanSubtract(d(floor(trialinds(:,170:250)/dwn))); hold on; plot((1:length(Snip))/(fs*dwn),Snip); plot((1:length(Snip))/(fs*dwn),mean(Snip'),'r','LineWidth',2);
plot((1:length(Snip))/(fs*dwn),median(Snip'),'k','LineWidth',2); hold off;

%% Aligning off of reaction time
trig = trigCondL; 
inds = floor(0:1:window*fs/1000);
trialinds = repmat((trig')/dwn, length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds(:,floor(trialinds(1,:))<=0) = []; trialinds(:,floor(trialinds(end,:))>length(accel)) = [];

% accel data
figure;
Snip = accel(:,2); Snip = u.meanSubtract(Snip(floor(trialinds(:,170:250)))); Snip = abs(Snip); hold on; plot(Snip); plot(mean(Snip'),'r','LineWidth',2);
plot(median(Snip'),'k','LineWidth',2); hold off;

% neural data
figure;
Snip = u.meanSubtract(d(floor(trialinds))); hold on; plot((1:size(Snip,1))/(fs),Snip); plot((1:size(Snip,1))/(fs),mean(Snip'),'r','LineWidth',2);
plot((1:size(Snip,1))/(fs),median(Snip'),'k','LineWidth',2); hold off;

