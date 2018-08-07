% % % Find the spike times
% % % convert the spike times into samples (onset)
% % % using the first spike time and the firs value of onset onset(1)
% % % 	- determine (if they exist) the 153 preceeding values the LFPs  and the 152 following values of the LFPs and average them for all 96 channels
% % % 	- return a 96 cell array containing the LFP averages
% % % 	- this array will be value 1 in a new array
% % % 	- perform this analysis for every spike of a given channel and cluster creating an array the LFPs for 96 channels at every spike
% 	- result should be a 96 x onset length matrix ie. 96x 2649
close all
tankpath = 'S:\~NeuroWest\Spanky\CycleTriggered-170710-143939\';
blockname = 'Spanky-180613-135031';
T1 = 0;
T2 = 0;

Cat = TDT2mat(tankpath, blockname, 'T1', T1, 'T2', T2, 'TYPE', [3,4]);
dog = Cat.streams.LFPs.data;
Fs = Cat.streams.LFPs.fs;
GOOD_CHANNELS = goodChan(blockname);
trig = GOOD_CHANNELS(:);

%###~~~~~~~~~~~~~~~~###
%###    SPIKERATIO  ###
%###~~~~~~~~~~~~~~~~###
% To determine what days are good and what days are bad by finding the
% the ratio of spikes that occur during noise events vs. total spikes.
% for j = 1:size(dog,1)
%     d = dog(j,:);
%     timeout = 1.5*fs;
%     fb = find(abs(d)>5e-4);
%     for i = 1:length(fb)
%         temp1 = round(fb(i)-timeout);
%         temp2 = round(fb(i)+timeout);
%         temp1(temp1<1)=1;
%         temp2(temp2>length(d))=length(d);
%         d(temp1:temp2) = nan;
%     end
% end
% d(isnan(d)) = [];
%###~~~~~~~~~~~~~~~~###



[pxx, freq] = pwelch(dog',[],[],[],Fs);
% pxx = (20*log10(pxx))';
pxx = pxx';

keep = freq > 5 & freq < 50;
pxx = pxx(:,keep);
freq = freq(keep);

magCC = corrcoef(pxx');
magTrigCC = magCC(trig,:);
% figure; plot(magTrigCC);

rawCC = corrcoef(dog');
rawTrigCC = rawCC(trig,:);
% figure; hold on; plot(rawTrigCC); 



magCombCC = mean(magCC(trig,:));

rawCombCC = mean(rawCC(trig,:));


for i = 1:96
    comb(i) = magCombCC(i)*rawCombCC(i);
end    

figure; plot(dog(81,:));
figure; plot(comb);
S = std(comb);







% Svec = ones([1,96]); S = std(rawTrigCC); Svec = Svec*S*1.5; 
% plot(Svec, 'r'); title('raw'); hold off
% bad(bad >(S*1.5)) = 0;
% bad = find(bad);

% figure; plot(comb);
% Svec = ones([1,96]);S = std(comb); Svec = Svec*S*1.5;
% hold on; plot(Svec, 'r'); title('combined'); hold off;
