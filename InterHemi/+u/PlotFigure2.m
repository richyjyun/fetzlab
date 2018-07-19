%% Figure 2. Anatomy of the movements
%     A.    Show unaligned raw accel, then aligned by RT detection
%     B.    Example average rectified acceleration aligned by RT detection
%     C.    average RP from contralateral hemisphere aligned with above
%     D.    Aligned beta power contra/ipsi
%     E.     Aligned gamma power contra/ipsi
%     F.     Indicate potential timing of stimulation with above
% % Just with Ubi's data
% clear
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbiFinal.mat');
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbiNeuro.mat');
% 
% d = 8; n = 7; % SL and SLNeuro indices


h = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); %set(h,'visible','off');

%% A. Show unaligned raw accel, then aligned by RT detection
window = [-0.5,1];
inds = window(1)*SL(d).fs:1:window(2)*SL(d).fs;

% unaligned
trig = SL(d).fs/1000*SL(d).lefttrials(:,1)';
trialinds = repmat(trig, length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds = floor(trialinds);
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(SL(d).accel_raw_l)) = [];
% figure; plot(inds/SL(d).fs,mean(SL(d).accel_raw_l(trialinds),2));

% aligned
trig = SL(d).fs/1000*SL(d).lefttrials(:,1)' + SL(d).rts_l';
trig(isnan(trig)) = [];
trialinds = repmat(trig, length(inds), 1) + repmat(inds(:), 1, size(trig,1));
trialinds = floor(trialinds);
trialinds(:,trialinds(1,:)<=0) = [];
trialinds(:,trialinds(end,:)>length(SL(d).accel_raw_l)) = [];
subplot(4,1,1); 
sig = SL(d).accel_raw_l(trialinds); sig = abs(sig-mean(sig)); sig = mean(sig,2);
plot(inds/SL(d).fs,sig - mean(sig),'linewidth',1.5);

%% B. Example average rectified acceleration aligned by RT detection
LFilter = u.FilterAccForMovementTimes(SL(d).accel_raw_l, SL(d).fs, 'richardson');
hold on; plot(inds/SL(d).fs,mean(LFilter(trialinds),2),'linewidth',1.5);
xlim(window); title('2A-2, 2B. Accelerometer, unfiltered & filtered')

%% C. average RP from contralateral hemisphere aligned with above\
smth = 50; % smoothing factor

rchn = 4; %RMC 2.5
lchn = 24; %LMC 4
subplot(4,1,2);
plot(SLNeuro(n).tneuro,smooth(squeeze(SLNeuro(n).LRP(1,rchn,:)),smth),'linewidth',1.5);
hold on; 
plot(SLNeuro(n).tneuro,smooth(squeeze(SLNeuro(n).RRP(1,rchn,:)),smth),'linewidth',1.5);
title('2C. RP'); xlim(window)

%% D. Aligned beta power contra/ipsi
subplot(4,1,3);
ninds = [find(SLNeuro(n).tneuro<window(1),1,'last'),find(SLNeuro(n).tneuro>window(2),1)];
sig = squeeze(SLNeuro(n).Lbeta(1,rchn,:));
%sig = sig-mean(sig(ninds(1):ninds(2)));
% sig = sig-mean(sig);
sig = smooth(sig,smth);
plot(SLNeuro(n).tneuro,sig,'linewidth',1.5);
sig = squeeze(SLNeuro(n).Rbeta(1,rchn,:));
%sig = sig-mean(sig(ninds(1):ninds(2)));
% sig = sig-mean(sig);
sig = smooth(sig,smth);
hold on; plot(SLNeuro(n).tneuro,sig,'linewidth',1.5);
xlim(window)
title('2D. Beta (blue contra, orange ipsi)');

%% E. Aligned gamma power contra/ipsi
subplot(4,1,4);
sig = squeeze(SLNeuro(n).Lgamma(1,rchn,:)); %sig = sig-mean(sig);
sig = smooth(sig,smth);
plot(SLNeuro(n).tneuro,sig,'linewidth',1.5);
sig = squeeze(SLNeuro(n).Rgamma(1,rchn,:)); %sig = sig-mean(sig);
sig = smooth(sig,smth);
hold on; plot(SLNeuro(n).tneuro,sig,'linewidth',1.5);
xlim(window)
title('2E. Gamma (blue contra, orange ipsi)');
xlabel('Time (s)');

%% F. Indicate potential timing of stimulation with above



