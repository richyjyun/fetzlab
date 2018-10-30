%% Set variables
clear; close all;
path = 'Y:\~NeuroWest\Spanky\Neurochip';
day = 'S_20180605_01';

%% Load meta data
[fpath,fname,Channels,fs,session_time] = getNCData(path,day);

start = 80*60; finish = 100*60; 
dwnfs = 500;
width = 8;
step = 10;
f = 0.5:0.1:50;

n = floor((finish-width)/step);
t = start; spectra = zeros(n,length(f));
ts = zeros(1,length(n));

load('Y:\~NeuroWest\Spanky\Neurochip\S_20180605_01\spk_18.mat');

lfpfs = 1000;
window = 0.05;
range = round(-window*lfpfs):1:round(window*lfpfs);

stlfp = zeros(n,length(range));

chns = 1:16; bad = [3,12,14,16]; chns(bad) = [];
tic 
parfor i = 1:n
    t = start+(i-1)*step;
    
    %     [data, ~, ~] = nc3data(channel, t, width, lfpfs, [], fname);
    %     data = bpfilt(data,[15,25],lfpfs,3);
    %
    %     trig = round((spkTime(spkTime>t & spkTime<t+width) - t)*lfpfs);
    %     trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
    %     trialinds(:,floor(trialinds(1,:))<=0) = [];
    %     trialinds(:,floor(trialinds(end,:))>length(data)) = [];
    %
    %     stlfp(i,:) = mean(data(trialinds),2)';
    
    [data, ~, ~] = nc3data(chns+16, t, width, dwnfs, [0,60], fname);
    
    [pxx,~] = pwelch(detrend(data),[],[],f,dwnfs);
    
%     F = fft(data(:,1));
%     
%     L = length(data);
%     P2 = abs(F/L);
%     P1 = P2(1:L/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);
%     
%     f = fs*(0:(L/2))/L;
%     plot(f,P1)
%     
    
%     temp = zeros(length(f),1);
%     for ch = 1:size(data,2)
%         s = spectrogram(data(:,ch),[],[],f,fs);
%         s = mag2db(abs(s));
%         temp = temp+mean(abs(s),2);
%     end
%     
%     temp = temp./size(data,2);
    
%     spectra(i,:) = temp';
    
    
    %
    spectra(i,:) = mean(pxx./sum(pxx),2);
    
    ts(i) = t+width/2;
    
%     t = t+step;
    
end
toc

spectra_raw = spectra;
save('Y:\~NeuroWest\Spanky\Neurochip\S_20180605_01\Spectra_8s.mat','spectra_raw','ts','start','dwnfs','width','step','finish','f','-v7.3');

spectra = spectra_raw./(sum(spectra_raw,2));

% delta = [1,4]; theta = [4,8]; alpha = [7.5,12.5]; beta = [13,25];
% delta = [findClosestInd(f,delta(1)),findClosestInd(f,delta(2))];
% theta = [findClosestInd(f,theta(1)),findClosestInd(f,theta(2))];
% alpha = [findClosestInd(f,alpha(1)),findClosestInd(f,alpha(2))];
% beta = [findClosestInd(f,beta(1)),findClosestInd(f,beta(2))];
% 
% 
% spectra = [mean(spectra(:,delta(1):delta(2)),2),...
%     mean(spectra(:,theta(1):theta(2)),2),...
%     mean(spectra(:,alpha(1):alpha(2)),2),...
%     mean(spectra(:,beta(1):beta(2)),2)];

for i = 1:size(spectra,1)
    spectra(i,:) = smooth(spectra(i,:),20);
end

Accel = nc3data(38, 0, step*n, 100, [], fname);
 
%% clustering

[U,S,V] = svd(spectra','econ');

lambda = (S.^2)/(length(S)-1); % variance
r = find(cumsum(diag(lambda)/sum(diag(lambda))) > 0.99,1); % get the rank that accounts for more than 95% of the variance
if(r<20)
    r = 20;
end

% Low dimensional space
Ur = U(:,1:r);
Sr = S(1:r,1:r);
Vr = V(:,1:r);

% figure; scatter3(Vr(:,1),Vr(:,2),Vr(:,3));

idx = kmeans(Vr,4,'Replicates',50,'MaxIter',1000);

% % Plot in svd space
% figure;
% for i = 1:max(idx)
%     scatter3(Vr(idx==i,1),Vr(idx==i,2),Vr(idx==i,3));
%     hold on;
% end
% 
% idx(idx==1) = 5;
% idx(idx==2) = 8;
% idx(idx==3) = 6;
% idx(idx==4) = 7;
% idx(idx==5) = 1;
% idx(idx==6) = 2;
% idx(idx==7) = 3;
% idx(idx==8) = 4;

% Individual spectra
figure; 
T = {'Movement','Awake at rest','REM-like','NREM-like'};
for i = 1:max(idx)
    subplot(2,2,i);
    x = f';
    y = mean(spectra(idx==i,:))';
    dy = std(spectra(idx==i,:))';
    fill([x;flipud(x)],[y-dy;flipud(y+dy)],[.8 .8 .8],'edgealpha',0);
    
    hold on;
    plot(x,y,'k','linewidth',2);
    xlim([0,50]); xlabel('Frequency (Hz)');
    ylabel('Normalized Spectral Density');
    title(T{i});
end
mtit('Individual Spectra');
print('-painters','-fillpage',gcf, '-dpsc2', ['F:\S\Packets\NC\',day,'_Classification.ps'], '-append');
close(gcf);

% Spectra stacked on top
figure; colors = get(gca,'colororder'); p = [];
for i = 1:max(idx)
    x = f';
    y = mean(spectra(idx==i,:))';
    dy = std(spectra(idx==i,:))';
    h = fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(i,:),'edgealpha',0);
    set(h,'facealpha',0.15);
    hold on;
    p(i) = plot(f,mean(spectra(idx==i,:)),'Color',colors(i,:),'linewidth',2);
end
xlabel('Frequency (Hz)');
ylabel('Normalized Spectral Density');
legend(p,T);
print('-opengl','-fillpage',gcf, '-dpsc2', ['F:\S\Packets\NC\',day,'_Classification.ps'], '-append');
close(gcf);

% Accel vs spectra
figure;
subplot(2,1,1)
for i = 1:max(idx)
    ind = idx==i;
    scatter(ts(ind),1500+(i*100)*ones(1,sum(ind)),'.');
    hold on;
end
plot((1:length(Accel))/100,Accel); xlabel('S');
yticks([]);
xlim([0,length(Accel)/200])
title('Classification and z-axis Acceleration over time')
subplot(2,1,2)
for i = 1:max(idx)
    ind = idx==i;
    scatter(ts(ind),1500+(i*100)*ones(1,sum(ind)),'.');
    hold on;
end
plot((1:length(Accel))/100,Accel); xlabel('S');
yticks([]);
xlim([length(Accel)/200,length(Accel)/100])
legend(T);
print('-painters','-fillpage',gcf, '-dpsc2', ['F:\S\Packets\NC\',day,'_Classification.ps'], '-append');
close(gcf);

% % Plot stLFP
% figure;
% trough = [];
% for i = 1:max(idx)
%     ind = idx==i;
%     plot(range/lfpfs,(mean(stlfp_18(ind,:))),'linewidth',2)
%     hold on;
%     trough(i) = find(zscore(mean(stlfp_24(ind,:))) == min(zscore(mean(stlfp_24(ind,:)))));
% end

%% Spiking dynamics
% ISI histogram
spkTime = spkTime_18;
ISI = diff(spkTime); ISI = [nan,ISI];

% without overlapping bins
edges = ts-width/2; edges = [edges,edges(end)+width];
Y = discretize(spkTime,edges);
figure; yl =0;
for i = 1:max(idx)
    bins = find(idx==i);
    binspk = find(any(Y==bins));
    subplot(2,2,i);
    histogram(ISI(binspk),0:0.001:0.1,'normalization','probability');
    title(T{i},'fontsize',8);
    xlabel('ms');
    yl = max([yl,ylim]);
end
for i = 1:max(idx)
    subplot(2,2,i);
    ylim([0,yl]);
end
mtit('ISI Histogram');
print('-painters','-fillpage',gcf, '-dpsc2', ['F:\S\Packets\NC\',day,'_Classification.ps'], '-append');
close(gcf);

% % with overlapping bins
% edges = ts-width/2; edges = [edges;edges+width];
% Y = discretizeOverlap(spkTime,edges);
% figure; yl =0;
% for i = 1:max(idx)
%     bins = find(idx==i);
%     binspk = unique([Y{bins}]);
%     subplot(2,2,i);
%     histogram(ISI(binspk),0:0.001:0.1,'normalization','probability');
%     title(T{i},'fontsize',8);
%     xlabel('ms');
%     yl = max([yl,ylim]);
% end
% for i = 1:max(idx)
%     subplot(2,2,i);
%     ylim([0,yl]);
% end
% mtit('ISI Histogram');

% autocorrelation
figure; yl=0;
for i = 1:max(idx)
    bins = find(idx==i);
    binspk = find(any(Y==bins));
    subplot(2,2,i);
    CrossCorr(spkTime(binspk), 'ts2',spkTime(binspk),'binsize', 0.001,'lag',[-0.1,0.1],'suppress_plot',0); %axis off;
    title(T{i});     xlabel('ms');
    yl = max([yl,ylim]);
end
for i = 1:max(idx)
    subplot(2,2,i);
    ylim([0,yl]);
end
mtit('Autocorrelation');
print('-painters','-fillpage',gcf, '-dpsc2', ['F:\S\Packets\NC\',day,'_Classification.ps'], '-append');
close(gcf);

% spike rate over time (with overlap)
edges = ts-width/2; edges = [edges;edges+width];
spkRate = histOverlap(spkTime,edges);%histcounts(spkTime,edges);
figure;
subplot(2,1,1)
for i = 1:max(idx)
    ind = idx==i;
    scatter(ts(ind),max(spkRate)+(i*200)*ones(1,sum(ind)),'.');
    hold on;
end
plot(ts,spkRate); xlabel('S');
yticks([]);
xlim([0,ts(end)/2]); ylim([0,max(spkRate)+(max(idx)+1)*200]);
title('Classification and Spike rate over time')
subplot(2,1,2)
for i = 1:max(idx)
    ind = idx==i;
    scatter(ts(ind),max(spkRate)+(i*200)*ones(1,sum(ind)),'.');
    hold on;
end
plot(ts,spkRate); xlabel('S');
yticks([]);
xlim([ts(end)/2,ts(end)]); ylim([0,max(spkRate)+(max(idx)+1)*200]);
legend(T);
print('-painters','-fillpage',gcf, '-dpsc2', ['F:\S\Packets\NC\',day,'_Classification.ps'], '-append');
close(gcf);

% spike rate histogram
figure; yl = 0;
for i = 1:max(idx)
    subplot(2,2,i);
    ind = idx==i;
    histogram(spkRate(ind)/60,0:1:80,'normalization','probability');
    xlabel('Hz');
    yl = max([yl,ylim]); xlim([0,80]);
    title(T{i},'fontsize',8);
end
for i = 1:max(idx)
    subplot(2,2,i);
    ylim([0,yl]);    
    ind = idx==i;
    med = median(spkRate(ind));
    hold on; plot([med,med],[0,yl],'k--','linewidth',1.5);
end
mtit('Spike rate Histogram');
print('-painters','-fillpage',gcf, '-dpsc2', ['F:\S\Packets\NC\',day,'_Classification.ps'], '-append');
close(gcf);

% cross correlation
% without overlap
spk1 = spkTime_18; spk2 = spkTime_24;
edges = ts-width/2; edges = [edges,edges(end)+width];
Y1 = discretize(spk1,edges);
Y2 = discretize(spk2,edges);
figure;
for i = 1:max(idx)
    bins = find(idx==i);
    binspk1 = find(any(Y1==bins));
    binspk2 = find(any(Y2==bins));
    
    subplot(2,2,i);
    CrossCorr(spk1(binspk1), 'ts2',spk2(binspk2),'binsize', 0.001,'lag',[-0.2,0.2],'suppress_plot',0); %axis off;
    title(T{i});
end
mtit('Cross Correlation');
print('-painters','-fillpage',gcf, '-dpsc2', ['F:\S\Packets\NC\',day,'_Classification.ps'], '-append');
close(gcf);

% % with overlap
% edges = ts-width/2; edges = [edges;edges+width];
% Y1 = discretizeOverlap(spk1,edges);
% Y2 = discretizeOverlap(spk2,edges);
% figure;
% for i = 1:max(idx)
%     bins = find(idx==i);
%     binspk1 = unique([Y1{bins}]);
%     binspk2 = unique([Y2{bins}]);
%     
%     subplot(2,2,i);
%     CrossCorr(spk1(binspk1), 'ts2',spk2(binspk2),'binsize', 0.001,'lag',[-0.2,0.2],'suppress_plot',0); %axis off;
%     title(T{i});
% end

% closest spike
spk2 = sort(spk2);
edges = [-Inf,mean([spk2(1:end-1); spk2(2:end)]),+Inf];
I = discretize(spk1,edges);
delta = spk1-spk2(I);

edges = ts-width/2; edges = [edges,edges(end)+width];
ISI = diff(spk1); ISI = [nan,ISI];
Y = discretize(spk1,edges);

edges = linspace(-0.1,0.1,1000);
figure; yl = 0;
for i = 1:max(idx)
    bins = find(idx==i);
    binspk = find(any(Y==bins));
    subplot(2,2,i);
    histogram(delta(binspk),edges,'normalization','probability');
    title(T{i},'fontsize',8); xlabel('ms');
    yl = max([yl,ylim]);
end
for i = 1:max(idx)
    subplot(2,2,i);
    ylim([0,yl]);
end
mtit('Closest spike between two cells');
print('-painters','-fillpage',gcf, '-dpsc2', ['F:\S\Packets\NC\',day,'_Classification.ps'], '-append');
close(gcf);


