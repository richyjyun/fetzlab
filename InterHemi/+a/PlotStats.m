function PlotStats(SL, ch_hem, fname)
%
% for an SL, calls the getstats function to collect median rts and then plots it like the old
% packets
%

% if strcmp(SL(1).monkey, 'Kato'),
%     ch_hem = ['L','L','L','L','R','R','R','R'];
% elseif strcmp(SL(1).monkey, 'Igor')
%     iL = strfind(chnm, 'ML');
%     iL = ~cellfun(@isempty, iL, 'unif', 1);
%     iR = strfind(chnm, 'MR');
%     iR = ~cellfun(@isempty, iR, 'unif', 1);
%     ch_hem(iL) = 'L';
%     ch_hem(iR) = 'R';
% end

ts = linspace(SL(1).TA_ts(1), SL(1).TA_ts(2), size(SL(1).TA_Lbeta,1));

[ipsiCSipsiRT, ipsiCScontRT, contCSipsiRT, contCScontRT, ipsiCScontRTLFPbeta, ipsiCSipsiRTLFPbeta, contCScontRTLFPbeta, contCSipsiRTLFPbeta, controlRTL, ipsiCScontbeta, ipsiCSipsibeta, contCScontbeta, contCSipsibeta, controlRTR, ipsiCScontRT_rtvt, ipsiCSipsiRT_rtvt, contCScontRT_rtvt, contCSipsiRT_rtvt] = ana.GetStats(SL, ch_hem);

figure('position', [0 0 850 1100], 'paperposition', [.01 .01 8.4 10.98])
 
%%%%%% 1
subplot(6, 4, [1 2 5 6])
line(repmat([1;2;3], 1, size(ipsiCScontRT,1)), ipsiCScontRT(:,[1 2 3])', 'color', [.6 .6 .6], 'linewidth', 1, 'marker', 'o'), hold on
bar(nanmedian(ipsiCScontRT(:,[1 2 3]), 1), .5, 'edgecolor', 'k', 'facecolor', 'none', 'linewidth', 1.5)
title('Ipsilateral CS, Contralateral RT'), set(gca, 'xtick', [1 2 3], 'xticklabel', {'Pre', 'Con', 'Pos'}), ylim([150, 425])

ax1 = subplot(6, 4, 9);
set(gca, 'colororder', [.7 .7 1; .4 .4 1; .1 .1 1])
line(repmat(ts(:), 1, 3), median(ipsiCScontRTLFPbeta(:,[1 3 5],:),3), 'linewidth', 1.5)
xlim([ts(1) ts(end)])
%ylim([.8*min(median(ipsiCScontRTLFPbeta(:,1,:),3)) 1.2*max(median(ipsiCScontRTLFPbeta(:,1,:),3))])
title('Ipsi \beta')

ax2 = subplot(6, 4, 10);
set(gca, 'colororder', [1 .7 .7; 1 .4 .4; 1 .1 .1])
line(repmat(ts(:), 1, 3), median(ipsiCScontRTLFPbeta(:,[2 4 6],:),3), 'linewidth', 1.5)
xlim([ts(1) ts(end)])
%ylim([.8*min(median(ipsiCScontRTLFPbeta(:,1,:),3)) 1.2*max(median(ipsiCScontRTLFPbeta(:,1,:),3))])
title('Cont \beta')

c = get([ax1 ax2], 'ylim');
set([ax1 ax2], 'ylim', [min([c{1}(1) c{2}(1)]) max([c{1}(2) c{2}(2)])]);

%%%%%% 2
subplot(6, 4, [3 4 7 8])
line(repmat([1;2;3], 1, size(ipsiCSipsiRT,1)), ipsiCSipsiRT(:,[1 2 3])', 'color', [.6 .6 .6], 'linewidth', 1, 'marker', 'o'), hold on
bar(nanmedian(ipsiCSipsiRT(:,[1 2 3]), 1), .5, 'edgecolor', 'k', 'facecolor', 'none', 'linewidth', 1.5)
title('Ipsilateral CS, Ipsilateral RT'), set(gca, 'xtick', [1 2 3], 'xticklabel', {'Pre', 'Con', 'Pos'}), ylim([150, 425])

ax1 = subplot(6, 4, 11);
set(gca, 'colororder', [.7 .7 1; .4 .4 1; .1 .1 1])
line(repmat(ts(:), 1, 3), median(ipsiCSipsiRTLFPbeta(:,[1 3 5],:),3), 'linewidth', 1.5)
xlim([ts(1) ts(end)])
%ylim([.8*min(median(ipsiCSipsiRTLFPbeta(:,1,:),3)) 1.2*max(median(ipsiCSipsiRTLFPbeta(:,1,:),3))])
title('Ipsi \beta')

ax2 = subplot(6, 4, 12);
set(gca, 'colororder', [1 .7 .7; 1 .4 .4; 1 .1 .1])
line(repmat(ts(:), 1, 3), median(ipsiCSipsiRTLFPbeta(:,[2 4 6],:),3), 'linewidth', 1.5)
xlim([ts(1) ts(end)])
%ylim([.8*min(median(ipsiCSipsiRTLFPbeta(:,1,:),3)) 1.2*max(median(ipsiCSipsiRTLFPbeta(:,1,:),3))])
title('Cont \beta')

c = get([ax1 ax2], 'ylim');
set([ax1 ax2], 'ylim', [min([c{1}(1) c{2}(1)]) max([c{1}(2) c{2}(2)])]);

%%%%%% 3
subplot(6, 4, [13 14 17 18])
line(repmat([1;2], 1, size(ipsiCScontRT,1)), ipsiCScontRT(:,[1 3])', 'color', [.6 .6 .6], 'linewidth', 1, 'marker', 'o'), hold on
bar(nanmedian(ipsiCScontRT(:, [1,3]), 1), .5, 'edgecolor', 'k', 'facecolor', 'none', 'linewidth', 1.5)
title('Ipsilateral CS, Contralateral RT'), set(gca, 'xtick', [1 2], 'xticklabel', {'Pre', 'Pos'}), ylim([150, 425])

ax1 = subplot(6, 4, 21);
set(gca, 'colororder', [.7 .7 1; .1 .1 1])
line(repmat(ts(:), 1, 2), median(ipsiCScontRTLFPbeta(:,[1 5],:),3), 'linewidth', 1.5)
xlim([ts(1) ts(end)])
title('Ipsi \beta')

ax2 = subplot(6, 4, 22);
set(gca, 'colororder', [1 .7 .7; 1 .1 .1])
line(repmat(ts(:), 1, 2), median(ipsiCScontRTLFPbeta(:,[2 6],:),3), 'linewidth', 1.5)
xlim([ts(1) ts(end)])
title('Cont \beta')

c = get([ax1 ax2], 'ylim');
set([ax1 ax2], 'ylim', [min([c{1}(1) c{2}(1)]) max([c{1}(2) c{2}(2)])]);

%%%%%% 4
subplot(6, 4, [15 16 19 20])
line(repmat([1;2], 1, size(ipsiCSipsiRT,1)), ipsiCSipsiRT(:, [1,3])', 'color', [.6 .6 .6], 'linewidth', 1, 'marker', 'o'), hold on
bar(nanmedian(ipsiCSipsiRT(:,[1,3]), 1), .5, 'edgecolor', 'k', 'facecolor', 'none', 'linewidth', 1.5)
title('Ipsilateral CS, Ipsilateral RT'), set(gca, 'xtick', [1 2], 'xticklabel', {'Pre', 'Pos'}), ylim([150, 425])

ax1 = subplot(6, 4, 23);
set(gca, 'colororder', [.7 .7 1; .1 .1 1])
line(repmat(ts(:), 1, 2), median(ipsiCSipsiRTLFPbeta(:,[1 5],:),3), 'linewidth', 1.5)
xlim([ts(1) ts(end)])
title('Ipsi \beta')

ax2 = subplot(6, 4, 24);
set(gca, 'colororder', [1 .7 .7; 1 .1 .1])
line(repmat(ts(:), 1, 2), median(ipsiCSipsiRTLFPbeta(:,[2 6],:),3), 'linewidth', 1.5)
xlim([ts(1) ts(end)])
title('Cont \beta')

c = get([ax1 ax2], 'ylim');
set([ax1 ax2], 'ylim', [min([c{1}(1) c{2}(1)]) max([c{1}(2) c{2}(2)])]);

print(gcf, '-dpsc2', fname, '-append')

clf

%% 1

signrank(contCScontRT(:,1),contCScontRT(:,2))
signrank(contCScontRT(:,1),contCScontRT(:,3))

signrank(contCSipsiRT(:,1),contCScontRT(:,2)) 
signrank(contCSipsiRT(:,1),contCScontRT(:,3))

subplot(6, 4, [1 2 5 6])
line(repmat([1;2;3], 1, size(contCScontRT,1)), contCScontRT(:,[1 2 3])', 'color', [.6 .6 .6], 'linewidth', 1, 'marker', 'o'), hold on
bar(nanmedian(contCScontRT(:,[1 2 3]), 1), .5, 'edgecolor', 'k', 'facecolor', 'none', 'linewidth', 1.5)
title('Contralateral CS, Contralateral RT'), set(gca, 'xtick', [1 2 3], 'xticklabel', {'Pre', 'Con', 'Pos'}), ylim([150, 425])

ax1 = subplot(6, 4, 9);
set(gca, 'colororder', [.7 .7 1; .4 .4 1; .1 .1 1])
line(repmat(ts(:), 1, 3), median(contCScontRTLFPbeta(:,[1 3 5],:),3), 'linewidth', 1.5)
xlim([ts(1) ts(end)])
%ylim([.8*min(median(contCScontRTLFPbeta(:,1,:),3)) 1.2*max(median(contCScontRTLFPbeta(:,1,:),3))])
title('Ipsi \beta')

ax2 = subplot(6, 4, 10);
set(gca, 'colororder', [1 .7 .7; 1 .4 .4; 1 .1 .1])
line(repmat(ts(:), 1, 3), median(contCScontRTLFPbeta(:,[2 4 6],:),3), 'linewidth', 1.5)
xlim([ts(1) ts(end)])
%ylim([.8*min(median(contCScontRTLFPbeta(:,1,:),3)) 1.2*max(median(contCScontRTLFPbeta(:,1,:),3))])
title('Cont \beta')

c = get([ax1 ax2], 'ylim');
set([ax1 ax2], 'ylim', [min([c{1}(1) c{2}(1)]) max([c{1}(2) c{2}(2)])]);

%% 2
subplot(6, 4, [3 4 7 8])
line(repmat([1;2;3], 1, size(contCSipsiRT,1)), contCSipsiRT(:,[1 2 3])', 'color', [.6 .6 .6], 'linewidth', 1, 'marker', 'o'), hold on
bar(nanmedian(contCSipsiRT(:,[1 2 3]), 1), .5, 'edgecolor', 'k', 'facecolor', 'none', 'linewidth', 1.5)
title('Contralateral CS, Ipsilateral RT'), set(gca, 'xtick', [1 2 3], 'xticklabel', {'Pre', 'Con', 'Pos'}), ylim([150, 425])

ax1 = subplot(6, 4, 11);
set(gca, 'colororder', [.7 .7 1; .4 .4 1; .1 .1 1])
line(repmat(ts(:), 1, 3), median(contCSipsiRTLFPbeta(:,[1 3 5],:),3), 'linewidth', 1.5)
xlim([ts(1) ts(end)])
%ylim([.8*min(median(contCSipsiRTLFPbeta(:,1,:),3)) 1.2*max(median(contCSipsiRTLFPbeta(:,1,:),3))])
title('Ipsi \beta')

ax2 = subplot(6, 4, 12);
set(gca, 'colororder', [1 .7 .7; 1 .4 .4; 1 .1 .1])
line(repmat(ts(:), 1, 3), median(contCSipsiRTLFPbeta(:,[2 4 6],:),3), 'linewidth', 1.5)
xlim([ts(1) ts(end)])
%ylim([.8*min(median(contCSipsiRTLFPbeta(:,1,:),3)) 1.2*max(median(contCSipsiRTLFPbeta(:,1,:),3))])
title('Cont \beta')

c = get([ax1 ax2], 'ylim');
set([ax1 ax2], 'ylim', [min([c{1}(1) c{2}(1)]) max([c{1}(2) c{2}(2)])]);

%%3
subplot(6, 4, [13 14 17 18])
line(repmat([1;2], 1, size(contCScontRT,1)), contCScontRT(:,[1 3])', 'color', [.6 .6 .6], 'linewidth', 1, 'marker', 'o'), hold on
bar(nanmedian(contCScontRT(:, [1,3]), 1), .5, 'edgecolor', 'k', 'facecolor', 'none', 'linewidth', 1.5)
title('Contralateral CS, Contralateral RT'), set(gca, 'xtick', [1 2], 'xticklabel', {'Pre', 'Pos'}), ylim([150, 425])

ax1 = subplot(6, 4, 21);
set(gca, 'colororder', [.7 .7 1; .1 .1 1])
line(repmat(ts(:), 1, 2), median(contCScontRTLFPbeta(:,[1 5],:),3), 'linewidth', 1.5)
xlim([ts(1) ts(end)])
title('Ipsi \beta')

ax2 = subplot(6, 4, 22);
set(gca, 'colororder', [1 .7 .7; 1 .1 .1])
line(repmat(ts(:), 1, 2), median(contCScontRTLFPbeta(:,[2 6],:),3), 'linewidth', 1.5)
xlim([ts(1) ts(end)])
title('Cont \beta')

c = get([ax1 ax2], 'ylim');
set([ax1 ax2], 'ylim', [min([c{1}(1) c{2}(1)]) max([c{1}(2) c{2}(2)])]);

%%4
subplot(6, 4, [15 16 19 20])
line(repmat([1;2], 1, size(contCSipsiRT,1)), contCSipsiRT(:, [1,3])', 'color', [.6 .6 .6], 'linewidth', 1, 'marker', 'o'), hold on
bar(nanmedian(contCSipsiRT(:,[1,3]), 1), .5, 'edgecolor', 'k', 'facecolor', 'none', 'linewidth', 1.5)
title('Contralateral CS, Ipsilateral RT'), set(gca, 'xtick', [1 2], 'xticklabel', {'Pre', 'Pos'}), ylim([150, 425])

ax1 = subplot(6, 4, 23);
set(gca, 'colororder', [.7 .7 1; .1 .1 1])
line(repmat(ts(:), 1, 2), median(contCSipsiRTLFPbeta(:,[1 5],:),3), 'linewidth', 1.5)
xlim([ts(1) ts(end)])
title('Ipsi \beta')

ax2 = subplot(6, 4, 24);
set(gca, 'colororder', [1 .7 .7; 1 .1 .1])
line(repmat(ts(:), 1, 2), median(contCSipsiRTLFPbeta(:,[2 6],:),3), 'linewidth', 1.5)
xlim([ts(1) ts(end)])
title('Cont \beta')

c = get([ax1 ax2], 'ylim');
set([ax1 ax2], 'ylim', [min([c{1}(1) c{2}(1)]) max([c{1}(2) c{2}(2)])]);

print(gcf, '-dpsc2', fname, '-append')

clf

%% plot control results
subplot(6, 4, [1 2 5 6])
line(repmat([1;2;3], 1, size(controlRTL,1)), controlRTL', 'color', [.6 .6 .6], 'linewidth', 1, 'marker', 'o'), hold on
bar(nanmedian(controlRTL, 1), .5, 'edgecolor', 'k', 'facecolor', 'none', 'linewidth', 1.5)
title('Control RT L'), set(gca, 'xtick', [1 2 3], 'xticklabel', {'Pre', 'Con', 'Pos'}), ylim([150, 425])

print(gcf, '-dpsc2', fname, '-append')

subplot(6, 4, [3 4 7 8])
line(repmat([1;2;3], 1, size(controlRTR,1)), controlRTR', 'color', [.6 .6 .6], 'linewidth', 1, 'marker', 'o'), hold on
bar(nanmedian(controlRTR, 1), .5, 'edgecolor', 'k', 'facecolor', 'none', 'linewidth', 1.5)
title('Control RT R'), set(gca, 'xtick', [1 2 3], 'xticklabel', {'Pre', 'Con', 'Pos'}), ylim([150, 425])

print(gcf, '-dpsc2', fname, '-append')

clf

%% plot RT as a function of time during conditioning session

subplot(2, 2, 1), ana.PlotRTvTime(ipsiCScontRT_rtvt), title('RT v T Ipsi CS Cont RT')
subplot(2, 2, 2), ana.PlotRTvTime(ipsiCSipsiRT_rtvt), title('RT v T Ipsi CS Ipsi RT')
subplot(2, 2, 3), ana.PlotRTvTime(contCScontRT_rtvt), title('RT v T Cont CS Cont RT')
subplot(2, 2, 4), ana.PlotRTvTime(contCSipsiRT_rtvt), title('RT v T Cont CS Ipsi RT')

print(gcf, '-dpsc2', fname, '-append')

clf

%% plot beta stats

subplot(2, 2, 1)
line(repmat([1;2;3], 1, size(ipsiCScontbeta,1)), ipsiCScontbeta', 'color', [.6 .6 .6], 'linewidth', 1, 'marker', 'o'), hold on
bar(nanmedian(ipsiCScontbeta, 1), .5, 'edgecolor', 'k', 'facecolor', 'none', 'linewidth', 1.5)
title('Ipsilateral Condition, Beta Cont to Stim'), set(gca, 'xtick', [1 2 3], 'xticklabel', {'Pre', 'Con', 'Pos'}), ylim([0 .2])

subplot(2, 2, 2)
line(repmat([1;2;3], 1, size(ipsiCSipsibeta,1)), ipsiCSipsibeta', 'color', [.6 .6 .6], 'linewidth', 1, 'marker', 'o'), hold on
bar(nanmedian(ipsiCSipsibeta, 1), .5, 'edgecolor', 'k', 'facecolor', 'none', 'linewidth', 1.5)
title('Ipsilateral Condition, Beta Ipsi to Stim'), set(gca, 'xtick', [1 2 3], 'xticklabel', {'Pre', 'Con', 'Pos'}), ylim([0 .2])

subplot(2, 2, 3)
line(repmat([1;2;3], 1, size(contCScontbeta,1)), contCScontbeta', 'color', [.6 .6 .6], 'linewidth', 1, 'marker', 'o'), hold on
bar(nanmedian(contCScontbeta, 1), .5, 'edgecolor', 'k', 'facecolor', 'none', 'linewidth', 1.5)
title('Contralateral Condition, Beta Cont to Stim'), set(gca, 'xtick', [1 2 3], 'xticklabel', {'Pre', 'Con', 'Pos'}), ylim([0 .2]) 

subplot(2, 2, 4)
line(repmat([1;2;3], 1, size(contCSipsibeta,1)), contCSipsibeta', 'color', [.6 .6 .6], 'linewidth', 1, 'marker', 'o'), hold on
bar(nanmedian(contCSipsibeta, 1), .5, 'edgecolor', 'k', 'facecolor', 'none', 'linewidth', 1.5)
title('Contralateral Condition, Beta Ipsi to Stim'), set(gca, 'xtick', [1 2 3], 'xticklabel', {'Pre', 'Con', 'Pos'}), ylim([0 .2]) 

print(gcf, '-dpsc2', fname, '-append')

clf

subplot(4, 2, 1) % change in RT from baseline during conditioning as a function of stim delay from target presentation
scatter(contCScontRT(:,4), -diff(contCScontRT(:,[1 2]),1,2), '*')
title('Contralateral CS, Contralateral RT'), xlabel('stim lag from target onset (ms)'), ylabel('change in RT from baseline in conditioning session')

subplot(4, 2, 2)
scatter(contCSipsiRT(:,4), -diff(contCSipsiRT(:,[1 2]),1,2), '*')
title('Contralateral CS, Ipsilateral RT'), xlabel('stim lag from target onset (ms)'), ylabel('change in RT from baseline in conditioning session')

subplot(4, 2, 3) % change in RP from baseline during postconditioning as a function fo stim delay from target presentation
scatter(contCScontRT(:,4), -diff(contCScontRT(:,[1 3]),1,2), '*')
title('Contralateral CS, Contralateral RT'), xlabel('stim lag from target onset (ms)'), ylabel('change in RT from baseline in conditioning session')

subplot(4, 2, 4)
scatter(contCSipsiRT(:,4), -diff(contCSipsiRT(:,[1 3]),1,2), '*')
title('Contralateral CS, Ipsilateral RT'), xlabel('stim lag from target onset (ms)'), ylabel('change in RT from baseline in conditioning session')

subplot(4, 2, 5) % change in RT from baseline during conditioning as a function of stim delay from target presentation
scatter(ipsiCScontRT(:,4), -diff(ipsiCScontRT(:,[1 2]),1,2), '*')
title('Ipsilateral CS, Contralateral RT'), xlabel('stim lag from target onset (ms)'), ylabel('change in RT from baseline in conditioning session')

subplot(4, 2, 6)
scatter(ipsiCSipsiRT(:,4), -diff(ipsiCSipsiRT(:,[1 2]),1,2), '*')
title('Ipsilateral CS, Ipsilateral RT'), xlabel('stim lag from target onset (ms)'), ylabel('change in RT from baseline in conditioning session')

subplot(4, 2, 7) % change in RP from baseline during postconditioning as a function fo stim delay from target presentation
scatter(ipsiCScontRT(:,4), -diff(ipsiCScontRT(:,[1 3]),1,2), '*')
title('Ipsilateral CS, Contralateral RT'), xlabel('stim lag from target onset (ms)'), ylabel('change in RT from baseline in conditioning session')

subplot(4, 2, 8)
scatter(ipsiCSipsiRT(:,4), -diff(ipsiCSipsiRT(:,[1 3]),1,2), '*')
title('Ipsilateral CS, Ipsilateral RT'), xlabel('stim lag from target onset (ms)'), ylabel('change in RT from baseline in conditioning session')

print(gcf, '-dpsc2', fname, '-append')

clf