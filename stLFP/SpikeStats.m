clear; close all;

%% Load
% 2. Epochs - Thd1 (threshold) and all corresponding values
% 3. Snips - Beta (snippets, 24kHz)
% 4. Streams - Mani (manipulandum, 3kHz), LFPs (3kHz), SUAs (24kHz), Filt
% (3kHz)
% 5. Scalars (all at 1Hz) - Trig (trigger channel), SUAc (SUA channels), Stim (Stim
% params and times) 
tankpath = 'Y:\~NeuroWest\Spanky\CycleTriggered-170710-143939\';
blockname = 'Spanky-171206-162516';
Dscm = TDT2mat([tankpath,blockname],'TYPE',2); Dscm = Dscm.epocs.Dscm; 
T = find(Dscm.data == 9980); T = T([1,3,5]); T = Dscm.onset([T-9950,T]); 
Snips = TDT2mat([tankpath,blockname],'T1',T(1,1),'T2',T(1,2),'TYPE',3,'STORE','eNe1'); Snips = Snips.snips.eNe1;

Chn = 83;
ind = Snips.chan == Chn & Snips.sortcode == 1;

% Plot snippets
snips = Snips.data(ind,:); figure; plot(snips','Color',[0.7,0.7,0.7]); hold on; plot(mean(snips),'k','LineWidth',2);

% Plot autocorrelation
[cor,lags] = u.CrossCorr(Snips.ts(ind), Snips.ts(ind),'binsize', 0.002,'lag',[-0.2,0.2],'suppress_plot',0);
figure; plot(lags,r);

% Plot ISI
intervals = diff(Snips.ts(ind)); intervals = intervals(intervals < 0.5); 
figure; histogram(intervals);

% Raster plot
figure;
for i = 1:96
    ind = Snips.chan == i & Snips.sortcode == 1;
    times = Snips.ts(ind); hold on;
    scatter(times,ones(length(times),1)*i,3,'k','filled');
end
