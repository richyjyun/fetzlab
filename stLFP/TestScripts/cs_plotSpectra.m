clear; %close all;

path = 'F:\S\Spectra\Post\32_32\';
files = dir([path,'*.mat']);

LFP = []; Amp = []; Phase = []; meanLFP = []; Ratio = [];
for i = 1%:length(files)
    load([path,files(i).name]);
    Amp = [Amp,amp];
    Phase = [Phase,phase];
    LFP = [LFP,sweeps-mean(sweeps)];
    meanLFP = [meanLFP,mean(sweeps-mean(sweeps),2)];
%     Ratio(:,:,i) = ratio;
end

window = 0.2; x = linspace(-window,window,size(LFP,1));
% % stLFP = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); 
% figure;
% stderr = std(meanLFP,0,2);
% fill([x,fliplr(x)],...
%     [mean(LFP,2)+stderr;flipud(mean(LFP,2)-stderr)],...
%     [.7,.7,.7],'LineStyle','none'); hold on;
% plot(x,mean(LFP,2),'k','linewidth',1.5); 
% xlim([-window,window]); xlabel('Time (s)');
% ylabel('Amplitude (V)');
% title(['stLFP, ',num2str(size(LFP,2)),' spikes, ',num2str(length(files)),' days'])

% Split them off into averages according to percentile of amplitude
collected = []; spikes = []; prcnt = 0:10:100; strength=[];
for p = 1:length(prcnt)-1
    for ang = 1:size(Phase,1)
        p1 = prctile(Amp(ang,:),prcnt(p)); p2 = prctile(Amp(ang,:),prcnt(p+1));
        inds = Amp(ang,:) > p1 & Amp(ang,:) < p2;
        spikes(p,ang) = sum(inds);
        % mean of circular values. kind of undoing the "angle"
        % function when saving phase
        r = sum(exp(1j*Phase(ang,inds)));
        collected(p,ang) = angle(r);
        strength(p,ang) = abs(r);
    end
end

% 0 is peak, pi/-pi are trough
map = hsv;

% Make regions lighter depending on number of spikes
% spikeProb = mean(Ratio,3);
% minv = -pi; maxv = pi;
% ncol = size(map,1);
% scaled = round(1+(ncol-1)*(collected'-minv)/(maxv-minv));
% RGB = ind2rgb(scaled,map);
% HSV = rgb2hsv(RGB);
% minv = min(spikeProb(:)); maxv = max(spikeProb(:));
% newVal = (spikeProb-minv)/(maxv-minv);
% HSV(:,:,2) = newVal;

% strength = sqrt(strength);

minv = -pi; maxv = pi;
ncol = size(map,1);
scaled = round(1+(ncol-1)*(collected'-minv)/(maxv-minv));
RGB = ind2rgb(scaled,map);
HSV = rgb2hsv(RGB);
minv = min(strength(:)); maxv = max(strength(:));
newVal = (strength'-minv)/(maxv-minv);
HSV(:,:,2) = newVal;

% Plot
% spectra = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); %set(gcf,'visible','off');
figure;
freq = 2:2:90;

subplot(2,1,1); imagesc(collected');
colormap hsv; set(gca,'YDir','normal'); c = colorbar; caxis([-pi,pi])
xticks(1.5:1:10.5); xticklabels(10:10:100); xlabel('Amplitude Envelope (Percentile)');
yticks(freq-0.5); yticklabels(freq*2); ylabel('Frequency (Hz)');
set(c,'Ticks',-pi:2*pi/8:pi,'TickLabels',{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'});
title('Average Phase')

subplot(2,1,2); imagesc(strength'); colorbar
set(gca,'YDir','normal');

subplot(2,1,2); imagesc(hsv2rgb(HSV));
colormap hsv; set(gca,'YDir','normal'); c = colorbar;
xticks(1.5:1:10.5); xticklabels(10:10:100); xlabel('Amplitude Envelope (Percentile)');
yticks(freq-0.5); yticklabels(freq*2); ylabel('Frequency (Hz)');
set(c,'Ticks',0:1/8:1,'TickLabels',{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'});
title('With Spiking Probability')
% 
% % Plot histogram of phases per frequency
% for ang = 30:40
%     figure;
%     for p = 1:length(prcnt)-1
%         color = [1,1,1] - p/length(prcnt);
%         p1 = prctile(Amp(ang,:),prcnt(p)); p2 = prctile(Amp(ang,:),prcnt(p+1));
%         inds = Amp(ang,:) > p1 & Amp(ang,:) < p2;
%         data = Phase(ang,inds);
%         [N,edges] = histcounts(data,-pi:pi/18:pi);
%         plot(edges(1:length(N))+pi/36,N,'color',color); hold on;
%     end
%     title([num2str(freq(ang)),'-',num2str(freq(ang+1)),'Hz'])
% end
% 
% Collected = [];
% for i = 1:length(files)
%     load([path,files(i).name]);
%     collected = [];
%     for p = 1:length(prcnt)-1
%         p1 = prctile(amp(:),prcnt(p)); p2 = prctile(amp(:),prcnt(p+1));
%         for ang = 1:size(phase,1)
%             inds = amp(ang,:) > p1 & amp(ang,:) < p2;
%             spikes(p,ang) = sum(inds);
%             % mean of circular values. kind of undoing the "angle"
%             % function when saving phase
%             r = sum(exp(1j*phase(ang,inds)));
%             collected(p,ang) = angle(r);
%         end
%     end
%     Collected(:,:,i) = collected;
% end
% 
% r = sum(exp(1j*Collected),3);
% avgPhase = angle(r)';
% imagesc(avgPhase);
% colormap hsv; set(gca,'YDir','normal'); c = colorbar; caxis([-pi,pi])
% xticks(1.5:1:10.5); xticklabels(10:10:100); xlabel('Amplitude Envelope (Percentile)');
% yticks(freq-0.5); yticklabels(freq*2); ylabel('Frequency (Hz)');
% set(c,'Ticks',-pi:2*pi/8:pi,'TickLabels',{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'});
% title('Average Phase')
% 
