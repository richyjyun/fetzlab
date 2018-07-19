%%
fname = '20170331_02';

dwn = 1;
[acc, trig1, ~, lefttrials, righttrials, fs, lefttrialsuccess, righttrialsuccess] = u.LoadTrainUbi([fname,'.f32'],dwn);

%% convert to ts
lefttrials = lefttrials/(fs);
righttrials = righttrials/(fs);
%%
tic
[data, fs, chnm, ~] = u.LoadGug(fname, 3);
toc
% RMC 3.5/6.5
% LMC 3.5/8.5

bpf = [1 200]; % chL/R transform: 1. bandpass filter (Hz) (around beta wave range)
[bbpf,abpf] = butter(1,bpf/(fs/2)); % 1st order bandpass (see daq_sapi_*)
Filter = filtfilt(bbpf,abpf,double(data));

%%
ts = 0:1/fs:(size(acc,1)-1)/fs;

ichs = [6 12 23 33];

window = [30 40];

tf = ts>window(1) &ts<window(2);

Filter(:,1) = u.FilterAccForMovementTimes(acc(:,1), fs, 'richardson');
Filter(:,2) = u.FilterAccForMovementTimes(acc(:,2), fs, 'richardson');

tmplt = lefttrials(lefttrials(:,1)>window(1)&lefttrials(:,2)<window(2),:);
tmprt = righttrials(righttrials(:,1)>window(1)&righttrials(:,2)<window(2),:);

figure,
for i = 2%:2%length(ichs)
%     subplot(2, 1, i),
    %     plot(ts(tf), Filter(tf, ichs(i))), hold on
    plot(ts(tf), acc(tf, i)), hold on
    plot(ts(tf), Filter(tf,i) + median(acc(tf, i)),'k','LineWidth',1.5), hold on
    plot(xlim,[SL(15).Max(i)/3 + median(acc(tf, i)),SL(15).Max(i)/3 + median(acc(tf, i))],'r','LineWidth',1.5), hold on
%     ylims = ylim;
%     plot(repmat(tmplt(:)', 2, 1),repmat(ylims(:), 1, size(tmplt,1)*2), 'b');
%     plot(repmat(tmprt(:)', 2, 1),repmat(ylims(:), 1, size(tmprt,1)*2), 'r');
end