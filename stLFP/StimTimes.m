clear;

%% Load in proper times
tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';
% tankpath = 'Y:\~NeuroWest\Spanky\CycleTriggered-170710-143939\';
% tankpath = 'Y:\~NeuroWest\Spanky\Connectivity-180207-131758\';
blockname = 'Spanky-180122-130528';
TT = TDT2mat([tankpath,blockname],'TYPE',2);
Dscm = TT.epocs.Dscm;
TT = TDT2mat([tankpath,blockname],'TYPE',5);
Stim = TT.scalars.Stim;

% just subtract what the delay should be (from sheet) and find closes
% spike, and get error from that, easiest way. 
delay = 0.005;

err = [];
for i = 1:length(Stim.ts)
   norm = Stim.ts(i)-delay-Dscm.onset;
   err(i) = norm(find(abs(norm)== min(abs(norm)),1));
end

outliers = err>(mean(err)+3*std(err)) | err<(mean(err)-3*std(err));
err(outliers) = [];
histogram(err);
delays = err+delay;

