%% Load in relevant data
tankpath = 'R:\Fetz Lab\neurowest\ARB_spankybackup\OP_DT1_052915\';
blockname = 'S20161014';

% spike times (root)
load('Y:\Richy\S20161014_Root.mat');

times(1) = 240; times(2) = 110*60;

T1 = times(1); T2 = times(2);
BxID = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','bxID'); BxID = BxID.streams.bxID;

%% parse Box ID (see google sheets)
jitter = 100;
vals = [7182,5134,8206,3086,1055,4110,9230,6158,10250];
start = BxID.data > vals(5)-jitter & BxID.data < vals(5) + jitter;
trialStart = find(diff(start)== -1);

trialID = zeros(length(trialStart),1);
bad = [];
for i = 1:length(trialStart)
    temp = BxID.data(trialStart(i)+1:end) > 100;
    ind = find(diff(temp) == 1, 1);
    if(isempty(ind))
        bad(end+1) = i;
        continue;
    end
    %     ind = find(BxID.data(trialStart(i)+1:end) > 100 ,1);
    id = BxID.data(trialStart(i)+ind+1);
    ind = find(vals-jitter < id & vals+jitter > id);
    if(isempty(ind) || length(ind)>1)
        continue;
    end
    trialID(i) = ind;
end

trialStart(bad) = []; trialID(bad) = [];

trialStart = trialStart/BxID.fs+times(1); % get time of trial start

window = [-1,2]; % window around trial start
binwidth = 0.05;
binedges = window(1):binwidth:window(2);
bincenter = binedges(1:end-1)+binwidth/2;

%% find direction of each cell
fname = 'F:\S\Packets\S20161014\Tuning.ps';
for c = 1:length(root.cells)
    spk = root.spike(root.cells(c)).ts;
    spkrate = zeros(9,length(binedges)-1);
    for t = 1:length(trialStart)
        bins = binedges + trialStart(t);
        [counts,~] = histcounts(spk,bins);
        spkrate(trialID(t),:) = spkrate(trialID(t),:)+counts;
    end
    figure('visible','off');
    for t = 1:size(spkrate,1)
        if t==5
            continue;
        end
        subplot(3,3,t)
        plot(bincenter,spkrate(t,:),'k');
        xlim(window);
    end
    subplot(3,3,5);
    text(0,0.5,[{['Channel ',num2str(root.spike(root.cells(c)).channel)]},...
        {[num2str(length(root.spike(root.cells(c)).i)),' Spikes']}]);
    axis off;
    print('-painters',gcf, '-dpsc2', fname, '-append');
    close(gcf)
end

%% synchrony?
spk1 = root.spike(1,90).ts;
spk2 = root.spike(3,58).ts;

spkdiff = zeros(1,length(spk1));
for i = 1:length(spk1)
    [~,idx] = min(abs(spk2-spk1(i)));
    spkdiff(i) = spk2(idx)-spk1(i);
end

tspkdiff = cell(1,9);
for t = 1:length(trialStart)
    idx = find(spk1 > trialStart(t) + window(1) & spk1 < trialStart(t) + window(2));
    tspkdiff{trialID(t)} = [tspkdiff{trialID(t)},spkdiff(idx)];
end

% histogram
figure; bins = linspace(-0.05,0.05,100);
for t = 1:size(spkrate,1)
    if t==5
        subplot(3,3,t)
        histogram(spkdiff,bins,'normalization','probability');
    else
        subplot(3,3,t)
        histogram(tspkdiff{t},bins,'normalization','probability');
    end
end

% normalized histogram
figure; bins = linspace(-0.05,0.05,100);
binc = bins(1:end-1)+(bins(2)-bins(1))/2;
[base,~] = histcounts(spkdiff,bins,'normalization','probability');
for t = 1:size(spkrate,1)
    if t==5
        continue;
    end
    subplot(3,3,t)
    [N,~]=histcounts(tspkdiff{t},bins,'normalization','probability');
    bar(binc,N-base,1);
end









