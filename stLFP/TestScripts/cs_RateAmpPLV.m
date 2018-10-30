clear; close all;

%% Load in proper times
tankpath = 'Y:\~NeuroWest\Spanky\RandomStim-180314-124242\';
% tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';
% tankpath = 'Y:\~NeuroWest\Spanky\IFNN\';
% tankpath = 'Y:\~NeuroWest\Spanky\CycleTriggered-170710-143939\';
% tankpath = 'Y:\~NeuroWest\Spanky\Connectivity-180207-131758\';
% tankpath = 'Y:\~NeuroWest\Spanky\SpankyUpstairs\';
blockname = 'Spanky-180806-163308';
TT = TDT2mat([tankpath,blockname],'TYPE',2);
Dscm = TT.epocs.Dscm;
[val,ind] = findpeaks(Dscm.data); 
val(end+1) = Dscm.data(end); ind(end+1) = length(Dscm.data);
ind = ind(val>1000); val = val(val>1000); 
tests = 1:2:length(val); 
times = [ind(tests)-val(tests),ind(tests)];%-val(tests)+val(2)];  % times = [ind(1) - val(1),ind(1)-val(1)+15000;ind(3)-val(3),ind(3)-val(3)+10000];
times(times==0) = 1;
times = Dscm.onset(times);

T1 = times(1,1); T2 = times(1,2);
T1 = 0; T2 = 15*60;
LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','verbose',false);
LFPs = LFPs.streams.LFPs;
Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu','verbose',false);
Snips = Snips.snips.eNe1;
filt = bpfilt(LFPs.data',[15,25],LFPs.fs,3)';

%% compare LFPs (beta oscillations?) 
% plot firing rate, beta amplitude, and plv

smth = 20;

chn = 93; code = 1;
ind = Snips.chan==chn & Snips.sortcode==code;
spk = Snips.ts(ind)-T1;
binwidth = 0.1; bins = 0:binwidth:((T2-T1));
rate = histcounts(spk,bins); rate = smooth(rate,smth);
figure; plot(bins(1:end-1)+binwidth/2,zscore(rate));

h = hilbert(filt(chn,:)); amp = abs(h); amp = smooth(amp,smth*round(binwidth*LFPs.fs)); 
% amp = amp(100:(end-100));
hold on; plot((1:length(amp))/LFPs.fs,zscore(amp));

binwidth = 0.1; bins = 0:binwidth:((T2-T1));
h_spk = h(round(spk*LFPs.fs)); h_spk = h_spk./abs(h_spk); label = discretize(spk,bins); 
plv = zeros(1,length(bins)-1); temp = accumarray(label,h_spk); div = accumarray(label,label); div = div./(1:length(div))';
plv(1:length(temp)) = abs(temp); plv(div~=0) = plv(div~=0)./div(div~=0)'; plv = smooth(plv,smth);
hold on; plot(bins(1:end-1)+binwidth/2,zscore(plv));

legend('Firing Rate','Amplitude','PLV');

% plot scatter plots (to see if there's any correlation) 
binwidth = 0.1; bins = 0:binwidth:((T2-T1)); %bins of firing rate
xfire = bins(1:end-1)+binwidth/2; xfire = round(xfire*LFPs.fs);
amp = amp(xfire);

% bad = find(amp > 2.5e-5);
% rate(bad) = []; amp(bad) = []; plv(bad) = [];
c = corrcoef([rate,amp,plv]);

f = figure; subplot(3,2,1); scatter(rate,amp,'k.');
p = polyfit(rate,amp,1);
yfit = polyval(p,rate);
hold on; plot(rate,yfit);
title({['Chn ',num2str(chn)],'Amplitude vs Firing Rate',['R: ',num2str(c(1,2))]});

subplot(3,2,3); scatter(rate,plv,'k.');
p = polyfit(rate,plv,1);
yfit = polyval(p,rate);
hold on; plot(rate,yfit);
title({'PLV vs Firing Rate',['R: ',num2str(c(1,3))]});

subplot(3,2,5); scatter(amp,plv,'k.');
p = polyfit(amp,plv,1);
yfit = polyval(p,amp);
hold on; plot(amp,yfit);
title({'PLV vs Amplitude',['R: ',num2str(c(2,3))]});

% other channel
chn = 28; code = 1;
ind = Snips.chan==chn & Snips.sortcode==code;
spk = Snips.ts(ind)-T1;
binwidth = 0.1; bins = 0:binwidth:((T2-T1));
rate2 = histcounts(spk,bins); rate2 = smooth(rate2,smth);
% figure; plot(bins(1:end-1)+binwidth/2,zscore(rate));

h = hilbert(filt(chn,:)); amp2 = abs(h); amp2 = smooth(amp2,smth*round(binwidth*LFPs.fs)); 
% amp = amp(100:(end-100));
% hold on; plot((1:length(amp))/LFPs.fs,zscore(amp));

binwidth = 0.1; bins = 0:binwidth:((T2-T1));
h_spk = h(round(spk*LFPs.fs)); h_spk = h_spk./abs(h_spk); label = discretize(spk,bins); 
plv2 = zeros(1,length(bins)-1); temp = accumarray(label,h_spk); div = accumarray(label,label); div = div./(1:length(div))';
plv2(1:length(temp)) = abs(temp); plv2(div~=0) = plv2(div~=0)./div(div~=0)'; plv2 = smooth(plv2,smth);
% hold on; plot(bins(1:end-1)+binwidth/2,zscore(plv));

% legend('Firing Rate','Amplitude','PLV');

% plot scatter plots (to see if there's any correlation) 
binwidth = 0.1; bins = 0:binwidth:((T2-T1)); %bins of firing rate
xfire = bins(1:end-1)+binwidth/2; xfire = round(xfire*LFPs.fs);
amp2 = amp2(xfire);

% bad = find(amp > 3e-5);
% rate(bad) = []; amp(bad) = []; plv(bad) = [];
c = corrcoef([rate2,amp2,plv2]);

figure(f); subplot(3,2,2); scatter(rate2,amp2,'k.');
p = polyfit(rate2,amp2,1);
yfit = polyval(p,rate2);
hold on; plot(rate2,yfit);
title({['Chn ',num2str(chn)],'Amplitude vs Firing Rate',['R: ',num2str(c(1,2))]});

subplot(3,2,4); scatter(rate2,plv2,'k.');
p = polyfit(rate2,plv2,1);
yfit = polyval(p,rate2);
hold on; plot(rate2,yfit);
title({'PLV vs Firing Rate',['R: ',num2str(c(1,3))]});

subplot(3,2,6); scatter(amp2,plv2,'k.');
p = polyfit(amp2,plv2,1);
yfit = polyval(p,amp2);
hold on; plot(amp2,yfit);
title({'PLV vs Amplitude',['R: ',num2str(c(2,3))]});

%% 
smth = 20;

chn = 93; code = 1;
ind = Snips.chan==chn & Snips.sortcode==code;
spk = Snips.ts(ind)-T1;
binwidth = 0.1; bins = 0:binwidth:((T2-T1));
rate = histcounts(spk,bins); rate = smooth(rate,smth);

h = hilbert(filt(chn,:)); amp = abs(h); amp = smooth(amp,smth*round(binwidth*LFPs.fs)); 
% amp = zscore(amp);

binwidth = 0.1; bins = 0:binwidth:((T2-T1));
h_spk = h(round(spk*LFPs.fs)); h_spk = h_spk./abs(h_spk); label = discretize(spk,bins); 
plv = zeros(1,length(bins)-1); temp = accumarray(label,h_spk); div = accumarray(label,label); div = div./(1:length(div))';
plv(1:length(temp)) = abs(temp); plv(div~=0) = plv(div~=0)./div(div~=0)'; plv = smooth(plv,smth);


chn = 28; code = 1;
ind = Snips.chan==chn & Snips.sortcode==code;
spk = Snips.ts(ind)-T1;
binwidth = 0.1; bins = 0:binwidth:((T2-T1));
rate2 = histcounts(spk,bins); rate2 = smooth(rate2,smth);

h = hilbert(filt(chn,:)); amp2 = abs(h); amp2 = smooth(amp2,smth*round(binwidth*LFPs.fs)); 
% amp2 = zscore(amp2);

binwidth = 0.1; bins = 0:binwidth:((T2-T1));
h_spk = h(round(spk*LFPs.fs)); h_spk = h_spk./abs(h_spk); label = discretize(spk,bins); 
plv2 = zeros(1,length(bins)-1); temp = accumarray(label,h_spk); div = accumarray(label,label); div = div./(1:length(div))';
plv2(1:length(temp)) = abs(temp); plv2(div~=0) = plv2(div~=0)./div(div~=0)'; plv2 = smooth(plv2,smth);


% figure; scatter(plv,plv2,'k.')

figure; a1 = subplot(3,1,1);  plot((1:length(amp))/LFPs.fs,amp); hold on; plot((1:length(amp2))/LFPs.fs,amp2); title('Amp')
a2 = subplot(3,1,2); plot(bins(1:end-1)+binwidth/2,rate); hold on;  plot(bins(1:end-1)+binwidth/2,rate2); title('Rate');
a3 = subplot(3,1,3); plot(bins(1:end-1)+binwidth/2,zscore(plv)); hold on; plot(bins(1:end-1)+binwidth/2,zscore(plv2)); title('PLV');
linkaxes([a1,a2,a3],'x');

binwidth = 0.5; bins = 0:binwidth:((T2-T1)); %bins of firing rate
xfire = bins(1:end-1)+binwidth/2; xfire = round(xfire*LFPs.fs);

amp_diff = amp2(xfire)-amp(xfire);
rate_diff = zscore(rate2)-zscore(rate);
plv_diff = zscore(plv2)-zscore(plv);

figure; scatter(plv_diff(inds),amp_diff(inds),'k.');

mdl = fitlm(plv_diff(inds),amp_diff(inds));
yfit = mdl.Coefficients{2,1} * plv_diff(inds) + mdl.Coefficients{1,1};

p = polyfit(plv_diff,amp_diff,1);
yfit = polyval(p,plv_diff);
hold on; plot(plv_diff(inds),yfit);

c = corrcoef(plv_diff,amp_diff);


%% contiguous diff

s = sign(plv_diff); edges = diff(s); 
lims = find(edges~=0); lims = [0;lims];
lims = [lims(1:end-1)+1,lims(2:end)];

d = []; 
for i = 1:length(lims)
    d(i,1) = mean(amp_diff(lims(i,1):lims(i,2)));
    d(i,2) = mean(plv_diff(lims(i,1):lims(i,2)));
end

figure; scatter(d(:,1),d(:,2),'k.')

mdl = fitlm(d(:,1),d(:,2));
yfit = mdl.Coefficients{2,1} * d(:,1) + mdl.Coefficients{1,1};
hold on; plot(d(:,1),yfit);

%% amplitude compare
chn = 28; code = 1;
ind = Snips.chan==chn & Snips.sortcode==code;
spk = Snips.ts(ind)-T1;
h = hilbert(filt(chn,:)); amp = abs(h); amp = smooth(amp,smth*round(binwidth*LFPs.fs)); 
amp = amp(round(spk*LFPs.fs));
ind = round(spk*10); ind(ind==0) = 1;
plv3 = plv(ind); plv4 = plv2(ind);
rate3 = rate(ind); rate4 = rate2(ind);

chn = 93;
h = hilbert(filt(chn,:)); amp2 = abs(h); amp2 = smooth(amp2,smth*round(binwidth*LFPs.fs)); 
amp2 = amp2(round(spk*LFPs.fs));


chn = 93; code = 1;
ind = Snips.chan==chn & Snips.sortcode==code;
spk = Snips.ts(ind)-T1;
h = hilbert(filt(chn,:)); amp3 = abs(h); amp3 = smooth(amp3,smth*round(binwidth*LFPs.fs)); 
amp3 = amp3(round(spk*LFPs.fs));
chn = 28;
h = hilbert(filt(chn,:)); amp4 = abs(h); amp4 = smooth(amp4,smth*round(binwidth*LFPs.fs)); 
amp4 = amp4(round(spk*LFPs.fs));
ind = round(spk*10); ind(ind==0) = 1;
plv5 = plv(ind); plv6 = plv2(ind);
rate5 = rate(ind); rate6 = rate2(ind);


% bins = linspace(0,1e-4,50);
figure; histogram(amp,'normalization','probability');
hold on; histogram(amp2,'normalization','probability');

figure; histogram(amp4,'normalization','probability');
hold on; histogram(amp3,'normalization','probability');

bins = linspace(0,1,50);
figure; histogram(plv3,bins,'normalization','probability');
hold on; histogram(plv4,bins,'normalization','probability');

figure; histogram(plv5,bins,'normalization','probability');
hold on; histogram(plv6,bins,'normalization','probability');


bins = linspace(0,6,50);
figure; histogram(rate3,bins,'normalization','probability');
hold on; histogram(rate4,bins,'normalization','probability');

figure; histogram(rate5,bins,'normalization','probability');
hold on; histogram(rate6,bins,'normalization','probability');


%% plv

chn = 93; code = 1;
ind = Snips.chan==chn & Snips.sortcode==code;
spk = Snips.ts(ind)-T1; ind = round(spk*10); ind(ind==0) = 1;

h = hilbert(filt(chn,:));
binwidth = 0.1; bins = 0:binwidth:((T2-T1));
h_spk = h(round(spk*LFPs.fs)); h_spk = h_spk./abs(h_spk); label = discretize(spk,bins); 
plv = zeros(1,length(bins)-1); temp = accumarray(label,h_spk); div = accumarray(label,label); div = div./(1:length(div))';
plv(1:length(temp)) = abs(temp); plv(div~=0) = plv(div~=0)./div(div~=0)'; plv = smooth(plv,smth);

chn = 28;
h = hilbert(filt(chn,:));
binwidth = 0.1; bins = 0:binwidth:((T2-T1));
h_spk = h(round(spk*LFPs.fs)); h_spk = h_spk./abs(h_spk); label = discretize(spk,bins); 
plv2 = zeros(1,length(bins)-1); temp = accumarray(label,h_spk); div = accumarray(label,label); div = div./(1:length(div))';
plv2(1:length(temp)) = abs(temp); plv2(div~=0) = plv2(div~=0)./div(div~=0)'; plv2 = smooth(plv2,smth);

bins = linspace(0,1,50);
figure; histogram(plv(ind),bins,'normalization','probability');
hold on; histogram(plv2(ind),bins,'normalization','probability');


chn = 28; code = 1;
ind = Snips.chan==chn & Snips.sortcode==code;
spk = Snips.ts(ind)-T1; ind = round(spk*10); ind(ind==0) = 1;

h = hilbert(filt(chn,:));
binwidth = 0.1; bins = 0:binwidth:((T2-T1));
h_spk = h(round(spk*LFPs.fs)); h_spk = h_spk./abs(h_spk); label = discretize(spk,bins); 
plv = zeros(1,length(bins)-1); temp = accumarray(label,h_spk); div = accumarray(label,label); div = div./(1:length(div))';
plv(1:length(temp)) = abs(temp); plv(div~=0) = plv(div~=0)./div(div~=0)'; plv = smooth(plv,smth);

chn = 93;
h = hilbert(filt(chn,:));
binwidth = 0.1; bins = 0:binwidth:((T2-T1));
h_spk = h(round(spk*LFPs.fs)); h_spk = h_spk./abs(h_spk); label = discretize(spk,bins); 
plv2 = zeros(1,length(bins)-1); temp = accumarray(label,h_spk); div = accumarray(label,label); div = div./(1:length(div))';
plv2(1:length(temp)) = abs(temp); plv2(div~=0) = plv2(div~=0)./div(div~=0)'; plv2 = smooth(plv2,smth);

bins = linspace(0,1,50);
figure; histogram(plv(ind),bins,'normalization','probability');
hold on; histogram(plv2(ind),bins,'normalization','probability');


% phase distribution
chn = 93; code = 1;
ind = Snips.chan==chn & Snips.sortcode==code;
spk = Snips.ts(ind)-T1; ind = round(spk*10); ind(ind==0) = 1;
h = hilbert(filt(chn,:)); p1 = angle(h); p1 = p1(ind);
plv1 = h(ind); plv1 = plv1./abs(plv1); plv1 = abs(sum(plv1)/length(plv1));
chn = 28;
h = hilbert(filt(chn,:)); p2 = angle(h); p2 = p2(ind);
plv2 = h(ind); plv2 = plv2./abs(plv2); plv2 = abs(sum(plv2)/length(plv2));

bins = linspace(-pi,pi,64);
figure; histogram(p1,bins,'normalization','probability');
hold on; histogram(p2,bins,'normalization','probability');

chn = 28; code = 1;
ind = Snips.chan==chn & Snips.sortcode==code;
spk = Snips.ts(ind)-T1; ind = round(spk*10); ind(ind==0) = 1;
h = hilbert(filt(chn,:)); p1 = angle(h); p1 = p1(ind);
chn = 93;
h = hilbert(filt(chn,:)); p2 = angle(h); p2 = p2(ind);

bins = linspace(-pi,pi,64);
figure; histogram(p1,bins,'normalization','probability');
hold on; histogram(p2,bins,'normalization','probability');

%% Jon
channel = 93;
figure; a1 = subplot(2,1,1); plot(filt(channel,:))
chn = 93; code = 1;
ind = Snips.chan==chn & Snips.sortcode==code;
spk = Snips.ts(ind)-T1;
hold on; s1 = scatter(round(spk*LFPs.fs),filt(channel,round(spk*LFPs.fs)),25,'filled','k');
chn = 28; code = 1;
ind = Snips.chan==chn & Snips.sortcode==code;
spk = Snips.ts(ind)-T1;
hold on; s2 = scatter(round(spk*LFPs.fs),filt(channel,round(spk*LFPs.fs)),25,'filled','r');
title('93 LFP'); legend([s1,s2],{'93','28'})

channel = 28;
a2 = subplot(2,1,2); plot(filt(channel,:))
chn = 93; code = 1;
ind = Snips.chan==chn & Snips.sortcode==code;
spk = Snips.ts(ind)-T1;
hold on; scatter(round(spk*LFPs.fs),filt(channel,round(spk*LFPs.fs)),25,'filled','k')
chn = 28; code = 1;
ind = Snips.chan==chn & Snips.sortcode==code;
spk = Snips.ts(ind)-T1;
hold on; scatter(round(spk*LFPs.fs),filt(channel,round(spk*LFPs.fs)),25,'filled','r')
title('28 LFP');

linkaxes([a1,a2],'xy');


