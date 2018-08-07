%% For plotting the phase gradient
% Loads in data, filters to beta, calculates instantaneous phase, loops to 
% find gradient at each time step. 

clear; close all;

%% Load in proper times
tankpath = 'Y:\~NeuroWest\Spanky\RandomStim-180314-124242\';
% tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';
% tankpath = 'Y:\~NeuroWest\Spanky\IFNN\';
% tankpath = 'Y:\~NeuroWest\Spanky\CycleTriggered-170710-143939\';
% tankpath = 'Y:\~NeuroWest\Spanky\Connectivity-180207-131758\';
blockname = 'Spanky-180717-135403';

T1 = times(2,1); T2 = times(2,2);
LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs;
Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu'); Snips = Snips.snips.eNe1;

fs = LFPs.fs;

trigChns = 81; codes = 1;
trig = (Snips.ts(Snips.chan == trigChns & Snips.sortcode == codes)' - T1)*fs;

LFPs.data = bpfilt(LFPs.data',[15,25],LFPs.fs,3)';

h = hilbert(LFPs.data')';

phase = angle(h(:,round(trig)));

bad = [1,18,20,51,52,95];

[avg,plv] = meanPhase(phase,2);
plv(bad) = nan;

grid = [];
for chn = 1:size(avg,1)
    [c,r,~] = GetWadeChannelPosition(chn);
    if(ismember(chn,bad))
        grid(r,c) = nan(1,1);
    else
        grid(r,c) = avg(chn);
    end
end

grid([1,1,10,10],[1,10,1,10]) = nan;


% f = figure; 
subplot(1,2,2); colormap hsv; clims = [-pi,pi];

im = imagesc(grid,clims); axis off;
set(im,'AlphaData',~isnan(grid)); colorbar

[fx,fy] = phaseGradient(grid,2);
hold on; quiver(1:10,1:10,fx,fy,'color','k','linewidth',1,'autoscalefactor',0.5);%,'autoscale','off');
%     hold on; contour(1:10,1:10,grid(:,:,i));
x = [1,1,6,6]; y = [1,6,1,6];

% calculating quadrants
% generalize matrix widths?
quadx = mat2cell(fx,[5,5],[5,5]);
quadx = cellfun(@(x) nanmean(x(:)),quadx);
quady = mat2cell(fy,[5,5],[5,5]);
quady = cellfun(@(x) nanmean(x(:)),quady);

quiver([3,8],[3,8],quadx,quady,'color','w','linewidth',5,'autoscalefactor',0.5);

c = colorbar; set(c,'ticks',[-pi,-pi/2,0,pi/2,pi],'ticklabels',[{'-\pi'},{'-\pi/2'},{'0'},{'\pi/2'},{'\pi'}])



