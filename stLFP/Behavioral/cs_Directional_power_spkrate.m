tankpath = 'R:\Fetz Lab\neurowest\ARB_spankybackup\OP_DT1_052915\';
blockname = 'S20160617';

times(1) = 240; times(2) = 110*60;

% Load in relevant data
T1 = times(1); T2 = times(2);
% Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu'); Snips = Snips.snips.eNeu;
BxID = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','bxID'); BxID = BxID.streams.bxID;
% LFPs = TDT2mat([tankpath,blockname],'T1',0,'T2',1,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs;

% parse Box ID (see google sheets) 
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

% chn = 58; code = 0;
% spks = Snips.ts(Snips.chan == chn & Snips.sortcode == code);

%% calculate lfp power and spkrate according to trial type
lfp = cell(9,1);
for i = 1:length(lfp)
    lfp{i} = cell(96,1);
end

spkrate = cell(9,1);
for i = 1:length(spkrate)
    spkrate{i} = cell(96,1);
end

for t = 1:length(trialStart)
    
    fprintf('Trial %d: %d%% Done\n',t,round(t*100/length(trialStart)));

    if(trialID(t) == 0 || trialID(t) ==5)
        continue;
    end
    
    T1 = trialStart(t)+window(1); T2 = trialStart(t)+window(2);
    
    LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','VERBOSE',false); 
    LFPs = LFPs.streams.LFPs;
    Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu','VERBOSE',false); 
    Snips = Snips.snips.eNeu;
    
    LFPs.data = bpfilt(LFPs.data',[15,25],LFPs.fs,3)';
    
    for ch = 1:96
        h = hilbert(LFPs.data(ch,:));
        lfp{trialID(t)}{ch} = [lfp{trialID(t)}{ch};abs(h)];
    end
    
    for ch = 1:96
        chn = ch; code = 0;
        spks = Snips.ts(Snips.chan == chn & Snips.sortcode == code);
        spks = spks - trialStart(t);
        
        [counts,~] = histcounts(spks,binedges);
        
        spkrate{trialID(t)}{ch} = [spkrate{trialID(t)}{ch};counts];
    end
end

fs = LFPs.fs;

save(['F:\S\Directional\',blockname,'_powerspkrate'],'lfp','spkrate','window','-v7.3');

%% plot

% plot all power
figure;
yl = [];
for ch = 1:96
    [c,r,e] = GetWadeChannelPosition(ch);
    subplot(10,10,(r-1)*10+c);
    
    for i = 1:length(lfp)
        if(i==5)
            continue;
        end
        hold on;
        plot(mean(lfp{i}{ch}))
    end
        
    yl(ch,:) = ylim;
    if(ch==chn)
        title(num2str(ch),'fontsize',7,'Color','r');
    else
        title(num2str(ch),'fontsize',7);
    end
    axis off;
end
legend({'1','2','3','4','5','6','7','8'})

% plot single channel power
ch = 58;
figure;
for t = 1:length(spkrate)
    if t==5
        continue;
    end
    subplot(3,3,t)
    plot(linspace(window(1),window(2),size(lfp{t}{ch},2)),mean(lfp{t}{ch}));
    xlim([window(1),window(2)]); ylim([1,2.5]*1e-5)
    box off;
end

% plot spkrate
ch = 58;
figure;
ttl = 1:8;
ttlind = 1;
for t = 1:length(spkrate)
    if t==5
        continue;
    end
    subplot(3,3,t)
    plot(binedges(1:end-1)+(binedges(2)-binedges(1))/2,mean(spkrate{t}{ch}));
    xlim([window(1),window(2)]); ylim([0,4])
    title(num2str(ttl(ttlind)))
    ttlind = ttlind+1;
    box off;
end


