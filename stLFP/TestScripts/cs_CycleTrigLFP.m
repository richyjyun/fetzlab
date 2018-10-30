clear; close all;

%% Load in proper times
tankpath = 'Y:\~NeuroWest\Spanky\RandomStim-180314-124242\';
blockname = 'Spanky-180914-145727';
TT = TDT2mat([tankpath,blockname],'TYPE',2);
Dscm = TT.epocs.Dscm;
[val,ind] = findpeaks(Dscm.data); 
val(end+1) = Dscm.data(end); ind(end+1) = length(Dscm.data);
ind = ind(val>1000); val = val(val>1000); 
tests =1 ;%1:2:length(val); 
times = [ind(tests)-val(tests),ind(tests)];%-val(tests)+val(2)];  % times = [ind(1) - val(1),ind(1)-val(1)+15000;ind(3)-val(3),ind(3)-val(3)+10000];
times(times==0) = 1;
times = Dscm.onset(times);

stimChn = 64; T1 = times(1); T2 = times(2);
Snips1 = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu','Channel',stimChn,'verbose',false);
Snips1 = Snips1.snips.eNe1;

blockname = 'Spanky-180802-101923';
TT = TDT2mat([tankpath,blockname],'TYPE',2);
Dscm = TT.epocs.Dscm;
[val,ind] = findpeaks(Dscm.data); 
val(end+1) = Dscm.data(end); ind(end+1) = length(Dscm.data);
ind = ind(val>1000); val = val(val>1000); 
tests =1 ;%1:2:length(val); 
times = [ind(tests)-val(tests),ind(tests)];%-val(tests)+val(2)];  % times = [ind(1) - val(1),ind(1)-val(1)+15000;ind(3)-val(3),ind(3)-val(3)+10000];
times(times==0) = 1;
times = Dscm.onset(times);

stimChn = 64; T1 = times(1); T2 = times(2);
Snips2 = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu','Channel',stimChn,'verbose',false);
Snips2 = Snips2.snips.eNe1;

%% 2 window Pre
figure; plot(Snips1.data(ceil(rand(1,1000)*length(Snips1.data)),:)');

win1 = getrect;
xlimit = [win1(1),win1(1)+win1(3)];
ylimit = [win1(2),win1(2)+win1(4)];
xbox = xlimit([1 1 2 2 1]);
ybox = ylimit([1 2 2 1 1]);
% 
% win2 = getrect;
% xlimit = [win2(1),win2(1)+win2(3)];
% ylimit = [win2(2),win2(2)+win2(4)];
% xbox2 = xlimit([1 1 2 2 1]);
% ybox2 = ylimit([1 2 2 1 1]);
% 
[xbox2,ybox2] = ginput;
xbox2(end+1) = xbox2(1);
ybox2(end+1) = ybox2(1);

idx = false(1,length(Snips1.data));
for i = 1:length(Snips1.data)
    [x1,y1] = polyxpoly(1:50,Snips1.data(i,:),xbox,ybox);
    [x2,y2] = polyxpoly(1:50,Snips1.data(i,:),xbox2,ybox2);
    if(~isempty(x1) && ~isempty(x2))
        idx(i) = true;
    end
end

cell1 = Snips1.data(idx,:);
cell2 = Snips1.data(~idx,:);

% split cell 2
figure; plot(cell2(ceil(rand(1,1000)*length(cell2)),:)');

win1 = getrect;
xlimit = [win1(1),win1(1)+win1(3)];
ylimit = [win1(2),win1(2)+win1(4)];
xbox = xlimit([1 1 2 2 1]);
ybox = ylimit([1 2 2 1 1]);

[xbox2,ybox2] = ginput;
xbox2(end+1) = xbox2(1);
ybox2(end+1) = ybox2(1);

idx2 = false(1,length(cell2));
for i = 1:length(cell2)
    [x1,y1] = polyxpoly(1:50,cell2(i,:),xbox,ybox);
    [x2,y2] = polyxpoly(1:50,cell2(i,:),xbox2,ybox2);
    if(~isempty(x1) && ~isempty(x2))
        idx2(i) = true;
    end
end

temp2 = cell2;
cell2 = temp2(idx2,:);
cell3 = temp2(~idx2,:);



figure; plot(cell1(ceil(rand(1,1000)*length(cell1)),:)')
figure; plot(cell2(ceil(rand(1,1000)*length(cell2)),:)')

ts1 = Snips1.ts(idx); figure;
CrossCorr(ts1, 'ts2',ts1,'binsize', 0.001,'lag',[-0.1,0.1],'suppress_plot',0); 

ind = find(~idx); ind = ind(idx2);
ts2 = Snips1.ts(ind); figure;
CrossCorr(ts2, 'ts2',ts2,'binsize', 0.001,'lag',[-0.1,0.1],'suppress_plot',0); 

% get lfps and plot stLFP
LFPs1 = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','verbose',false);
LFPs1 = LFPs1.streams.LFPs; fs = LFPs1.fs;

trig = round((ts1-T1)*fs)';
window = [-0.05,0.05];
range = round(window(1)*fs:1:window(2)*fs);

trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(LFPs1.data)) = [];

stLFPs1 = zeros(size(LFPs1.data,1),length(range));
for j = 1:size(LFPs1.data,1)
    d = bpfilt(LFPs1.data(j,:),[15,25],LFPs.fs,3)';
    d = d(floor(trialinds));
    d = d - mean(d);
    d = mean(d,2);
    stLFPs1(j,:) = d;%zscore(d);
end

fig = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]);
yl(i,:,:) = nan(size(LFPs1.data,1),2);
for j = 1:size(LFPs1.data,1)
    [c,r,~] = GetWadeChannelPosition(j);
    subplot(10,10,(r-1)*10+c);
    
    hold on;
    plot(range/fs,stLFPs1(j,:),'k'); hold on;
    yl(i,j,:) = ylim; xlim([window(1),window(2)]);
end

hold on; plot(mean(cell1),'k','linewidth',2);
figure; CrossCorr(ts1, 'ts2',ts1,'binsize', 0.001,'lag',[-0.1,0.1],'suppress_plot',0); 

%% 2 window post
figure; plot(Snips2.data(ceil(rand(1,1000)*length(Snips2.data)),:)'); %ylim([-1e-4,1e-4])
win1 = getrect;
xlimit = [win1(1),win1(1)+win1(3)];
ylimit = [win1(2),win1(2)+win1(4)];
xbox = xlimit([1 1 2 2 1]);
ybox = ylimit([1 2 2 1 1]);

% win2 = getrect;
% xlimit = [win2(1),win2(1)+win2(3)];
% ylimit = [win2(2),win2(2)+win2(4)];
% xbox2 = xlimit([1 1 2 2 1]);
% ybox2 = ylimit([1 2 2 1 1]);

[xbox2,ybox2] = ginput;
xbox2(end+1) = xbox2(1);
ybox2(end+1) = ybox2(1);


idx = false(1,length(Snips2.data));
for i = 1:length(Snips2.data)
    [x1,y1] = polyxpoly(1:50,Snips2.data(i,:),xbox,ybox);
    [x2,y2] = polyxpoly(1:50,Snips2.data(i,:),xbox2,ybox2);
    if(~isempty(x1) && ~isempty(x2))
        idx(i) = true;
    end
end

cell2 = Snips2.data(idx,:);
cell2 = Snips2.data(~idx,:);

figure; plot(cell2(ceil(rand(1,1000)*length(cell1)),:)')
figure; plot(cell2(ceil(rand(1,1000)*length(cell2)),:)')

ts2 = Snips2.ts(idx); figure;
CrossCorr(ts1, 'ts2',ts1,'binsize', 0.001,'lag',[-0.1,0.1],'suppress_plot',0); 

ts2 = Snips2.ts(~idx); figure; 
CrossCorr(ts2, 'ts2',ts2,'binsize', 0.001,'lag',[-0.1,0.1],'suppress_plot',0); 


% get lfps and plot stLFP
LFPs2 = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','verbose',false);
LFPs2 = LFPs2.streams.LFPs; fs = LFPs2.fs;

trig = round((ts2-T1)*fs)';
window = [-0.05,0.05];
range = round(window(1)*fs:1:window(2)*fs);

trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(LFPs2.data)) = [];

stLFPs2 = zeros(size(LFPs2.data,1),length(range));
for j = 1:size(LFPs2.data,1)
    d = bpfilt(LFPs2.data(j,:),[15,25],LFPs.fs,3)';
    d = d(floor(trialinds));
    d = d - mean(d);
    d = mean(d,2);
    stLFPs2(j,:) = d;%zscore(d);
end

hold on;
yl(i,:,:) = nan(size(LFPs.data,1),2);
for j = 1:size(LFPs.data,1)
    [c,r,~] = GetWadeChannelPosition(j);
    subplot(10,10,(r-1)*10+c);
    
    hold on;
    plot(range/fs,stLFPs2(j,:),'r'); hold on;
    yl(i,j,:) = ylim; xlim([window(1),window(2)]);
end

figure; plot(mean(cell1));
figure; CrossCorr(ts1, 'ts2',ts1,'binsize', 0.001,'lag',[-0.1,0.1],'suppress_plot',0); 



%% chn32
Snips3 = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu','Channel',32,'verbose',false);
Snips3 = Snips3.snips.eNe1;
ts3 = Snips3.ts(Snips3.sortcode==1);

figure; CrossCorr(ts2, 'ts2',ts3,'binsize', 0.002,'lag',[-0.1,0.1],'suppress_plot',0); 




%% Clean figure

YLIM = [-2e-6,2e-6];
for j = 1:size(LFPs.data,1)
    
    [c,r,~] = GetWadeChannelPosition(j);
    subplot(10,10,(r-1)*10+c);
    ylim(YLIM);
    hold on;
    line([0 0], YLIM, 'linestyle', ':', 'color', [.5 .5 .5]);
    axis off
    sub_pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
    if(j == 91)
        title(num2str(j),'fontsize',7,'Color','r')
    else
        title(num2str(j),'fontsize',7);
    end
end






