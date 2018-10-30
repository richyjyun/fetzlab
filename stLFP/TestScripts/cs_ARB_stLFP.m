tankpath = 'R:\Fetz Lab\neurowest\ARB_spankybackup\OP_DT1_052915\';
blockname = 'S20161014';

T1 = 0; T2 = 0;


% % get LFPs
LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs;
% define variables
window = 0.5;
fs = LFPs.fs;
range = round(-window*fs:1:window*fs);

range = round(-0.01*fs:1:0.04*fs);

% % get all snippets from spike sorting
% Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu','Channel',trigChn);
% Snips = Snips.snips.eNeu;

% stLFPs = zeros(96,length(range));
% window = 0.05;
% for s = 1:length(Snips.ts)
%     t1 = Snips.ts(s) - window; t2 = Snips.ts(s) + window;
%     if(t1 < T1), continue; end
%     
%     LFPs = TDT2mat([tankpath,blockname],'T1',t1,'T2',t2,'TYPE',4,'STORE','LFPs','verbose',false); 
% 
%     if(isempty(stLFPs))
%         stLFPs = LFPs.streams.LFPs.data;
%     else
%         stLFPs = stLFPs + LFPs.streams.LFPs.data(:,1:length(range));
%     end
% end

% stLFPs = stLFPs./length(Snips.ts);

% % % filter
% % if(~isempty(filt))
% %     LFPs.data = bpfilt(LFPs.data',filt,LFPs.fs,3)';
% % end

% data = Snips.data'; ts = Snips.ts;
% bad = any(data>1e-4 | data<-1e-4);
% data(:,bad) = [];
% ts(bad) = [];
% 
% [p,spk,ts,noise] = getSortingParams(double(data),fs*8,[],ts);
% [spkTime,spk,noise,spkInd] = sortChannel(double(Snips.data'),p);
% 
% 
% 
% 
% data = Snips.data;
% tind = Snips.data(:,8) > window1(2,2) & Snips.data(:,8) < window1(1,2)...
%     & Snips.data(:,15) > window2(2,2) & Snips.data(:,15) < window2(1,2);
% tind(bad) = false; 

load('Y:\Richy\S20161014_Root.mat');
% plx_ts = root.spike(1,12).ts;
% snip_ts = Snips.ts(Snips.chan == 12);

% temp = abs(snip_ts-plx_ts');
% [~,ind] = min(temp);

% trig = snip_ts(ind); 

chn = 62;
trigChn = chn;

Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu','Channel',chn);
Snips = Snips.snips.eNeu;

trig = root.spike(1,chn).ts; trig = trig.*fs;
trig = round(trig)';
trig = trig(end-10000:end);

trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];

% get stlfp
% LFPs.data = bpfilt(LFPs.data',filt,LFPs.fs,3)';

stLFPs = zeros(size(LFPs.data,1),length(range));
for j = 1:size(LFPs.data,1)
    d = LFPs.data(j,:);%bpfilt(LFPs.data(j,:)',[15,25],LFPs.fs,3)';
    d = d(floor(trialinds));
    d = d - mean(d);
    d = median(d,2);
    stLFPs(j,:) = d;%zscore(d);
end

% % plot
% figure('position', [0 0 1100 850], 'paperposition', [0 0 11 8.5], 'paperorientation', 'landscape');
% yl = nan(size(LFPs.data,1),2);
% for j = 1:size(LFPs.data,1)
%     [c,r,~] = GetWadeChannelPosition(j);
%     subplot(10,10,(r-1)*10+c);
%     
%     hold on;
%     plot(range/fs,stLFPs(j,:),'k'); hold on;
%     yl(j,:) = ylim; xlim([-window,window]);
% end

%% binned
data = LFPs.data(chn,:);
data = data(trialinds);
data = data-mean(data); 

bins = linspace(-3e-4,3e-4,50);

binned = [];
for i = 1:size(data,1)
    binned(:,i) = histcounts(data(:,i),bins);
end

figure; imagesc(binned);

%% Clean plot
YLIM = [-6e-6,6e-6];%nanmedian(yl)*1.5;
% plot
figure('position', [0 0 1100 850], 'paperposition', [0 0 11 8.5], 'paperorientation', 'landscape');
yl = nan(size(LFPs.data,1),2);
for j = 1:size(LFPs.data,1)
    [c,r,~] = GetWadeChannelPosition(j);
    subplot(10,10,(r-1)*10+c);
    
    hold on;
    if(j~=trigChn)
    plot(range/fs,stLFPs(j,:),'k'); hold on;
    else
            plot(range/fs,stLFPs(j,:),'r'); hold on;
    end
    xlim([-0.05,0.05]);
    yl(j,:) = ylim; 
end

for j = 1:size(LFPs.data,1)
    
    [c,r,~] = GetWadeChannelPosition(j);
    subplot(10,10,(r-1)*10+c);
    if(j~=trigChn)
        ylim(YLIM);
    end
%     xlim([-0.01,0.04])
    hold on;
    line([0 0], yl(j,:), 'linestyle', ':', 'color', [.5 .5 .5]);
    axis off
    sub_pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
%     if(j == trigChn)
%         title(num2str(j),'fontsize',7,'Color','r')
%     else
%         title(num2str(j),'fontsize',7,'Color','k')
%     end
    
end

annotation('textbox', [0.1, 0.95, 0.5, 0.01], 'linestyle', 'none', 'string', 'Recreation, xlim [-0.5,0,5], ylim [-6e-6,6e-6], trig all')

%% spectragram

data = LFPs.data(chn,:);
data = data(trialinds);
data = data-mean(data); 

% Define parameters for Chronux
params.tapers = [3,5]; params.Fs = fs; params.fpass = [0,100]; params.trialave = 1;
movingwin = [60,10];

% Moving window spectrum
[S,t,f] = mtspecgramc(LFPs.data(chn,:),movingwin,params);

[S,f] = mtspectrumc(LFPs.data(chn,1:8e5),params);
figure; plot(f,S);


%% Plot snips and autocorrelation (from last epoch, shouldn't matter)
ind = Snips.chan == trigChn & Snips.sortcode == 0;
snips = Snips.data(ind,:); sample = floor(linspace(1,size(snips,1), 100));
subplot(10,10,91);
plot(snips','Color',[0.7,0.7,0.7]); hold on; plot(mean(snips),'k','LineWidth',2);
ylim([min(min(snips(sample,:))),max(max(snips(sample,:)))]); axis off;
sub_pos = get(gca,'position'); % get subplot axis position
set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
subplot(10,10,100);
CrossCorr(cluster_class(cluster_class(:,1) == 1,2), 'ts2',cluster_class(cluster_class(:,1) == 1,2),'binsize', 0.001,'lag',[-window,window],'suppress_plot',0); axis off;
sub_pos = get(gca,'position'); % get subplot axis position
set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height

%% Title
axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
txt = sprintf('%s; Trigger Channel %d, Code %d; %d Triggers, +-%dms window, Ylim %e, %e',...
    blockname,trigChn,0,sum(ind),round(window*1000),YLIM(1),YLIM(2));
text(0.5, 0.1,txt,...
    'HorizontalAlignment' ,'center','VerticalAlignment', 'top','fontsize',9)










