clear; %close all;

Theta = cell(1,2);
PLV = cell(1,2);

path = 'F:\S\Spectra\32_64\';
files = dir([path,'*.mat']);

Amp = []; Phase = [];
for i = 1:length(files)
    load([path,files(i).name]);
    Amp = [Amp,amp];
    Phase = [Phase,phase];
end

[Theta{1},PLV{1}] = getPercentiles(Amp,Phase,0:10:100);

path = 'F:\S\Spectra\Post\32_64\';
files = dir([path,'*.mat']);

Amp = []; Phase = [];
for i = 1:length(files)
    load([path,files(i).name]);
    Amp = [Amp,amp];
    Phase = [Phase,phase];
end

[Theta{2},PLV{2}] = getPercentiles(Amp,Phase,0:10:100);


T = Theta{2}-Theta{1};
bad = abs(T)>pi;
T(bad) = 2*pi-abs(T(bad));

P = abs(PLV{2}-PLV{1});

map1 = hot; map1 = map1(24:end,:);
map2 = fliplr(map1);
map = [map1(1:end-1,:);flipud(map2)];

% thresh = prctile(P(:),70);
% inds = P>=thresh;
% M = max(max(abs(T(inds))));
mint = -pi; maxt = pi;

ncol = size(map,1);
scaled = round(1+(ncol-1)*(T'-mint)/(maxt-mint));
RGB = ind2rgb(scaled,map);
HSV = rgb2hsv(RGB);
minp = min(P(:)); maxp = max(P(:));
newVal = (P'-minp)/(maxp-minp);
HSV(:,:,2) = newVal;

freq = 2:2:90;
colormap(map); imagesc(hsv2rgb(HSV));
set(gca,'YDir','normal'); c = colorbar;
xticks(1.5:1:10.5); xticklabels(10:10:100); xlabel('Amplitude Envelope (Percentile)');
yticks(freq-0.5); yticklabels(freq*2); ylabel('Frequency (Hz)');
% set(c,'Ticks',0:0.5:1,'TickLabels',{'-\pi/4',0,'\pi/4'});
set(c,'Ticks',0:1/8:1,'TickLabels',{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'});


