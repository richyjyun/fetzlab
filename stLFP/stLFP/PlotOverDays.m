clear; figure;

%% Load in proper times
tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';
% tankpath = 'Y:\~NeuroWest\Spanky\CycleTriggered-170710-143939\';
% tankpath = 'Y:\~NeuroWest\Spanky\Connectivity-180207-131758\';
blockname = 'Spanky-180204-155626';
TT = TDT2mat([tankpath,blockname],'TYPE',2);
Dscm = TT.epocs.Dscm;
[val,ind] = findpeaks(Dscm.data); 
ind = ind(val>1000); val = val(val>1000); 
tests = 1:2:length(val); 
times = [ind(tests)-val(tests),ind(tests)-val(tests)+val(1)];  % times = [ind(1) - val(1),ind(1)-val(1)+15000;ind(3)-val(3),ind(3)-val(3)+10000];
times = Dscm.onset(times);
time = times(1,:);


%% Pre
window = 0.05; %50 ms window, change as needed

T1 = time(1) - window;
T2 = time(2) + window;
LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs;
Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1'); Snips = Snips.snips.eNe1;

    LFPs.data = bpfilt(LFPs.data',[10,50],LFPs.fs,3)';


fs = LFPs.fs;
range = round(-window*fs:1:window*fs);
yl = [];
trigChn = 32;
trig = (Snips.ts(Snips.chan == trigChn & Snips.sortcode == 1)' - T1)*fs;

trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];

d = LFPs.data(32,:);
d = d(floor(trialinds));
d = d - mean(d);
d = mean(d,2);

subplot(2,5,1); plot(d,'k','linewidth',1.5)
subplot(2,5,5); plot(d,'linewidth',1.5); hold on;

d = LFPs.data(64,:);
d = d(floor(trialinds));
d = d - mean(d);
d = mean(d,2);

subplot(2,5,6); plot(d,'k','linewidth',1.5)
subplot(2,5,10); plot(d,'linewidth',1.5); hold on;


%% After stim
time = times(2,:);
T1 = time(1) - window;
T2 = time(2) + window;
LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs;
Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1'); Snips = Snips.snips.eNe1;
    LFPs.data = bpfilt(LFPs.data',[10,50],LFPs.fs,3)';

fs = LFPs.fs;
range = round(-window*fs:1:window*fs);
% yl = zeros(length(trigChns),2);

trigChn = 32;
trig = (Snips.ts(Snips.chan == trigChn & Snips.sortcode == 1)' - T1)*fs;

trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];

d = LFPs.data(32,:);
d = d(floor(trialinds));
d = d - mean(d);
d = mean(d,2);

subplot(2,5,2); plot(d,'k','linewidth',1.5)
subplot(2,5,5); plot(d,'linewidth',1.5); hold on;

d = LFPs.data(64,:);
d = d(floor(trialinds));
d = d - mean(d);
d = mean(d,2);

subplot(2,5,7); plot(d,'k','linewidth',1.5)
subplot(2,5,10); plot(d,'linewidth',1.5); hold on;


%% Load in proper times
tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';
% tankpath = 'Y:\~NeuroWest\Spanky\CycleTriggered-170710-143939\';
% tankpath = 'Y:\~NeuroWest\Spanky\Connectivity-180207-131758\';
blockname = 'Spanky-180206-144732';
TT = TDT2mat([tankpath,blockname],'TYPE',2);
Dscm = TT.epocs.Dscm;
[val,ind] = findpeaks(Dscm.data); 
ind = ind(val>1000); val = val(val>1000); 
tests = 1:2:length(val); 
times = [ind(tests)-val(tests),ind(tests)-val(tests)+val(1)];  % times = [ind(1) - val(1),ind(1)-val(1)+15000;ind(3)-val(3),ind(3)-val(3)+10000];
times = Dscm.onset(times);
time = times(1,:);


%% Pre
window = 0.05; %50 ms window, change as needed

T1 = time(1) - window;
T2 = time(2) + window;
LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs;
Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1'); Snips = Snips.snips.eNe1;
    LFPs.data = bpfilt(LFPs.data',[10,50],LFPs.fs,3)';

fs = LFPs.fs;
range = round(-window*fs:1:window*fs);
trigChn = 32;
trig = (Snips.ts(Snips.chan == trigChn & Snips.sortcode == 1)' - T1)*fs;

trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];

d = LFPs.data(32,:);
d = d(floor(trialinds));
d = d - mean(d);
d = mean(d,2);

subplot(2,5,3); plot(d,'k','linewidth',1.5)
subplot(2,5,5); plot(d,'linewidth',1.5); hold on;

d = LFPs.data(64,:);
d = d(floor(trialinds));
d = d - mean(d);
d = mean(d,2);

subplot(2,5,8); plot(d,'k','linewidth',1.5)
subplot(2,5,10); plot(d,'linewidth',1.5); hold on;

%% Load in proper times
tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';
% tankpath = 'Y:\~NeuroWest\Spanky\CycleTriggered-170710-143939\';
% tankpath = 'Y:\~NeuroWest\Spanky\Connectivity-180207-131758\';
blockname = 'Spanky-180208-133323';
TT = TDT2mat([tankpath,blockname],'TYPE',2);
Dscm = TT.epocs.Dscm;
[val,ind] = findpeaks(Dscm.data); 
ind = ind(val>1000); val = val(val>1000); 
tests = 1:2:length(val); 
times = [ind(tests)-val(tests),ind(tests)-val(tests)+val(1)];  % times = [ind(1) - val(1),ind(1)-val(1)+15000;ind(3)-val(3),ind(3)-val(3)+10000];
times = Dscm.onset(times);
time = times(1,:);


%% Pre
window = 0.05; %50 ms window, change as needed

T1 = time(1) - window;
T2 = time(2) + window;
LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs;
Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1'); Snips = Snips.snips.eNe1;
    LFPs.data = bpfilt(LFPs.data',[10,50],LFPs.fs,3)';

fs = LFPs.fs;
range = round(-window*fs:1:window*fs);
trigChn = 32;
trig = (Snips.ts(Snips.chan == trigChn & Snips.sortcode == 1)' - T1)*fs;

trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];

d = LFPs.data(32,:);
d = d(floor(trialinds));
d = d - mean(d);
d = mean(d,2);

subplot(2,5,4); plot(d,'k','linewidth',1.5)
subplot(2,5,5); plot(d,'linewidth',1.5); hold on;

d = LFPs.data(64,:);
d = d(floor(trialinds));
d = d - mean(d);
d = mean(d,2);

subplot(2,5,9); plot(d,'k','linewidth',1.5)
subplot(2,5,10); plot(d,'linewidth',1.5); hold on;


%% Cleaning the figure
yl = [];
for i = 6:10
    subplot(2,5,i)
    yl(end+1,:) = ylim;
end

yl = [min(yl(:,1)),max(yl(:,2))];
z = ceil(length(d)/2);

for i = 1:10
    subplot(2,5,i)
    if(i<6)
        ylim(yl*2)
    else
        ylim(yl)
    end
    xlim([1,length(d)]);
    hold on;
    plot([z,z],yl,'--','Color',[.7,.7,.7],'linewidth',1.5);
    axis off;
end

