% cs_CheckRandDelays
close all
fname = ['20170406_02';'20170412_03'];
bins = -250:100:450; % ms

figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 10])

for i = 1:size(fname,1)
    file = fname(i,:);
    [accdata, trig1, trig2, lefttrials, righttrials, fs, lefttrialsuccess, righttrialsuccess] = u.LoadTrainUbi([file,'.i16']); % load data
    
    S.Date = file;
    S.Animal = 'Ubi';
    S.Condition = 'Contra_R';
    S.Stim_Delay = '0';
    S.StimHemi = 'L';
    S.lefttrials = lefttrials;
    S.righttrials = righttrials;
    S.fs = fs;
    S.accel_raw_r = double(accdata(:,2));
    S.accel_raw_l = double(accdata(:,1));
    S.lefttrialsuccess = lefttrialsuccess;
    S.righttrialsuccess = righttrialsuccess;
    S.trig1 = trig1;
    S.trig2 = trig2;
    
    S = u.AppendThreshold(S);
    S = a.RemoveTwitch(S);
    S = a.AppendReactionTimes(S);
    S = a.AppendNormalizedDelay(S,1); 
    
    rts = nan(1,length(trig1));
    lts = rts;
    for j = 1:length(trig1)
        norm = abs(righttrials(:,1)-trig1(j));
        ind = find(norm == min(norm),1);
        rts(j) = S.rts_r(ind);
        
        norm = abs(righttrials(:,1)-trig1(j));
        ind = find(norm == min(norm),1);
        lts(j) = S.rts_l(ind);
    end
    
    % Plotting
%     figure
%     plot(S.Delays)
%     figure
%     hist(S.Delays)
    
%     figure;
%     scatter(S.Delays,rts)
%     xlabel('Normalized Delay (ms)')
%     ylabel('Reaction Time (ms)')
%     ylim([100 460])
%     title(S.Date);

% figure of change against delay
[c, ~, w] = histcounts(S.Delays, bins);
mean_rts = nan(length(unique(w))-1,1);
std_rts = mean_rts;
for ii = unique(w(:)')
    if ii==0, continue; end
    mean_rts(ii) = mean(rts(w==ii));
    std_rts(ii) = std(rts(w==ii));%/sqrt(c(ii));
end
subplot(2,2,1)
errorbar(bins(1:end-1)+diff(bins(1:2))/2, mean_rts, std_rts, 'linewidth', 2), hold on
ylim([100 475])
xlim([-250 450])
ylabel('RT Right (ms)')
xlabel('LHEM STIM DELAY (ms)')

% left trials
[c, ~, w] = histcounts(S.Delays, bins);
mean_rts = nan(length(unique(w))-1,1);
std_rts = mean_rts;
for ii = unique(w(:)')
    if ii==0, continue; end
    mean_rts(ii) = mean(lts(w==ii));
    std_rts(ii) = std(lts(w==ii));%/sqrt(c(ii));
end
subplot(2,2,2)
errorbar(bins(1:end-1)+diff(bins(1:2))/2, mean_rts, std_rts, 'linewidth', 2), hold on
ylim([100 475])
xlim([-250 450])
ylabel('RT Left (ms)')
xlabel('LHEM STIM DELAY (ms)')

% change against number of stims
windowsize = 20;
bounds = 1:windowsize:length(S.Delays);
meanL = nan(length(bounds),1);
stdL = meanL;
meanR = meanL;
stdR = meanL;
for ii = 1:length(bounds)-1
   tmp = S.rts_l(S.lefttrials(:,1)>S.trig1(bounds(ii)) & S.lefttrials(:,1)<=S.trig1(bounds(ii+1)));
   meanL(ii) = nanmean(tmp);
   stdL(ii) = nanstd(tmp);
   tmp = S.rts_r(S.righttrials(:,1)>S.trig1(bounds(ii)) & S.righttrials(:,1)<=S.trig1(bounds(ii+1)));
   meanR(ii) = nanmean(tmp);
   stdR(ii) = nanstd(tmp);
end

meanL(end) = [];
stdL(end) = [];
meanR(end) = [];
stdR(end)= [];

subplot(2,2,3), title(['Right RT against stim number (mean of ' num2str(windowsize) ')'])
errorbar(bounds(1:end-1)+diff(bounds(1:2))/2, meanR, stdR, 'linewidth', 2), hold on

ylim([100 475])
xlim([-10 270])
ylabel('RT RIGHT (ms)')
xlabel('N STIM')

subplot(2,2,4), title(['Left RT against stim number (mean of ' num2str(windowsize) ')'])
errorbar(bounds(1:end-1)+diff(bounds(1:2))/2, meanL, stdL, 'linewidth', 2), hold on

ylim([100 475])
xlim([-10 270])
ylabel('RT LEFT (ms)')
xlabel('N STIM')

end

print(gcf, '-dpsc2', '../packets/controlrandstim.ps')