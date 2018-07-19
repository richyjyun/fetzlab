% cs_plotTwitchRemoval

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat')

i = 6;

window = 600;
window = floor(window*SL(i).fs/1000);
figure
inds = -window/3:1:window;
trialinds = repmat(SL(i).trig1'.*SL(i).fs./1000, length(inds), 1) + repmat(inds(:), 1, size(SL(i).trig1,1));

RFilter = u.FilterAccForMovementTimes(SL(i).accel_raw_r, SL(i).fs, 'richardson');
LFilter = u.FilterAccForMovementTimes(SL(i).accel_raw_l, SL(i).fs, 'richardson');

Snips = RFilter(trialinds);
subplot(2,2,1); plot(inds,Snips); title('Twitch Right Trials');
Snips = LFilter(trialinds);
subplot(2,2,2); plot(inds,Snips); title('Twitch Left Trials');

SL = a.RemoveTwitch(SL);

RFilter = u.FilterAccForMovementTimes(SL(i).accel_raw_r, SL(i).fs, 'richardson');
LFilter = u.FilterAccForMovementTimes(SL(i).accel_raw_l, SL(i).fs, 'richardson');

Snips = RFilter(trialinds);
subplot(2,2,3); plot(inds,Snips); title('No Twitch Right Trials');
Snips = LFilter(trialinds);
subplot(2,2,4); plot(inds,Snips); title('No Twitch Left Trials');
