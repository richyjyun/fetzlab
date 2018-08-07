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
times = [0,20*60];

T1 = times(1,1); T2 = times(1,2);
LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs;
fs = LFPs.fs;

LFPs.data = bpfilt(LFPs.data',[15,25],LFPs.fs,3)';

h = hilbert(LFPs.data')';

phase = angle(h);

bad = [1,3,18,67,20,22,15,29,31,95];

grid = [];
for chn = 1:size(phase,1)
    [c,r,~] = GetWadeChannelPosition(chn);
    if(ismember(chn,bad))
        grid(r,c,:) = nan(1,1,size(phase,2));
    else
        grid(r,c,:) = phase(chn,:);
    end
end
grid([1,1,10,10],[1,10,1,10],:) = nan;

vidfile = VideoWriter(fullfile('F:\S\Packets\Video',[blockname,'_PhaseGrad2.mp4']),'MPEG-4');
open(vidfile);

f = figure; colormap hsv; clims = [-pi,pi]; 
Fx = []; Fy = [];
for i = 1:size(grid,3)

%     im = imagesc(grid(:,:,i),clims); axis off;
%     set(im,'AlphaData',~isnan(grid(:,:,i)))
%     
    [fx,fy] = phaseGradient(grid(:,:,i),2);
    Fx(:,:,end+1) = fx; Fy(:,:,end+1) = fy;
%     hold on; quiver(1:10,1:10,fx,fy,'color','k','linewidth',1,'autoscalefactor',0.5);%,'autoscale','off');
% %     hold on; contour(1:10,1:10,grid(:,:,i));
%     x = [1,1,6,6]; y = [1,6,1,6];
%     
%     % calculating quadrants
%     % generalize matrix widths? 
%     quadx = mat2cell(fx,[5,5],[5,5]);
%     quadx = cellfun(@(x) nanmean(x(:)),quadx);
%     quady = mat2cell(fy,[5,5],[5,5]);
%     quady = cellfun(@(x) nanmean(x(:)),quady);
% 
%     quiver([3,8],[3,8],quadx,quady,'color','w','linewidth',5,'autoscalefactor',0.5);
%     
%     c = colorbar; set(c,'ticks',[-pi,-pi/2,0,pi/2,pi],'ticklabels',[{'-\pi'},{'-\pi/2'},{'0'},{'\pi/2'},{'\pi'}])
% 
%     hold off;
%     title(i);%/fs);
%     
%     if(~ishandle(f))
%         break;
%     end
%     
%     drawnow
% %     pause(0.5);
% %     
% %     F = getframe(gcf);
% %     writeVideo(vidfile,F);
% %     
end


f = figure; colormap hsv; clims = [-pi,pi]; 

highBeta = find(abs(h(81,:)) > 2*std(abs(h(81,:))));

highPh = grid(:,:,highBeta); highPh = meanPhase(highPh,3);
highFx = nansum(Fx(:,:,highBeta),3);
highFy = nansum(Fy(:,:,highBeta),3);

im = imagesc(highPh,clims); axis off;
set(im,'AlphaData',~isnan(grid(:,:,i)))

hold on; quiver(1:10,1:10,highFx,highFy,'color','k','linewidth',1,'autoscalefactor',0.5);%,'autoscale','off');

% close(vidfile);



difference = exp(1j*phase(78,:))-exp(1j*phase(81,:));
avg = angle(sum(difference));






