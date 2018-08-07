tankpath = 'R:\Fetz Lab\neurowest\ARB_spankybackup\OP_DT1_052915\';
blockname = 'S20160617';

times(1) = 10; times(2) = 100*60;
chn = 53; code = 0;

% Load in relevant data
T1 = times(1); T2 = times(2);
% Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu'); Snips = Snips.snips.eNeu;
% BxID = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','bxID'); BxID = BxID.streams.bxID;
LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','Channel',52); LFPs = LFPs.streams.LFPs;
Mani = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','Mani'); Mani = Mani.streams.Mani;
Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu'); Snips = Snips.snips.eNeu;

%% trying with instantaneous velocity
Mani.data = [-Mani.data(5,:);Mani.data(3,:)]; % 5 is FE, 3 is RU. same as current set up in needing to flip FE.
% Mani.data = interp1((1:length(Mani.data))/Mani.fs,Mani.data',(1:length(LFPs.data))/LFPs.fs)';
smoothwin = 1*Mani.fs;
Mani.data(1,:) = smooth(Mani.data(1,:)',smoothwin)';
Mani.data(2,:) = smooth(Mani.data(2,:)',smoothwin)';

% figure;
% for i = 10001:length(Mani.data(1,:))
%     subplot(2,1,1)
%     plot(Mani.data(1,i-5000:i),Mani.data(2,i-5000:i));
%     subplot(2,1,2)
%     plot(IDs(1:ceil(i*BxID.fs/Mani.fs)));
%     pause(0.005);
% end

vel = diff(Mani.data')';
ang = atan2(vel(2,:),vel(1,:));
spd = sqrt(vel(2,:).^2+vel(1,:).^2);

edges = -pi:pi/8:pi;

direction = discretize(ang,edges);
trialDir = [14,15;...
            12,13;...
            10,11;...
            1,16;...
            nan,nan;...
            8,9;...
            2,3;...
            4,5;...
            6,7];
        
trialID = [];
for i = 1:length(trialDir)
   trialID(direction == trialDir(i,1) | direction == trialDir(i,2)) = i; 
end

% trialID(spd<0.5*std(spd)) = 0;

trig = Snips.ts(Snips.chan == chn & Snips.sortcode == code);

trigID = trialID(round((trig-T1)*Mani.fs)');

trig = round((trig-T1)*LFPs.fs)';

win = 0.05; % 50ms windows for stLFPs
range = round(-win*LFPs.fs:1:win*LFPs.fs);

trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
badInd = floor(trialinds(1,:))<=0 | floor(trialinds(end,:))>length(LFPs.data);
trialinds(:,badInd) = []; trigID(badInd) = [];

%raw stlfp
d = LFPs.data;
d = d(floor(trialinds));
d = d - mean(d);

figure;     spkInd = find(range==0);
avg = []; plv = []; pow = [];
for i = 1:9
    subplot(3,3,i);
    inds = trigID==i;
    temp = d(:,inds);
    plot(range,mean(temp,2)); title(num2str(sum(inds)));
end

%filtered stlfp
Filt = bpfilt(double(LFPs.data),[15,25],LFPs.fs,4);
d = Filt;
d = d(floor(trialinds));
for i = 1:9
    subplot(3,3,i);
    inds = trigID==i;
    temp = d(:,inds);
    
    % phase locking value
    if i==5 || isempty(temp)
        continue;
    end
    
    plot(range,mean(temp,2));
    title(num2str(size(temp,2)));
    
    h = hilbert(temp);
    pow(:,i) = mean(abs(h),2);
    h = angle(h);
    phases = h(spkInd,:);
    S = sum(exp(1j*phases));
    avg(i) = angle(S);
    plv(i) = abs(S)/length(phases);
    
end

good = 1:9; good(5) = [];
theta = pi:-pi/4:-pi;
inds = [4,1,2,3,6,9,8,7,4];

figure; polarplot(theta,avg(inds)); title('Avg Phase');
% rlim([2.2,2.8])
% rlim([min(avg(inds)),max(avg(inds))]);
figure; polarplot(theta,plv(inds)); title('PLV');

figure;
for i = 1:9
    subplot(3,3,i);
    plot(pow(:,i));
end

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

chn = 58; code = 0;
% spks = Snips.ts(Snips.chan == chn & Snips.sortcode == code);
% neu = Snips.data(Snips.chan == 58 & Snips.sortcode == 0,:);
win = 0.05; % 50ms windows for stLFPs
fs = LFPs.fs;
range = round(-win*fs:1:win*fs);

%% Trial triggered spike hist
figure;
for i = 1:9
    trig = trialStart(trialID == i);
    
end

%% Trial triggered LFP


%% get stLFP
stLFP = cell(9,1);
for i = 1:length(stLFP)
    stLFP{i} = cell(96,1);
end

used = [];
for t = 1:length(trialStart)
    
    fprintf('Trial %d: %d%% Done\n',t,round(t*100/length(trialStart)));

    if(trialID(t) == 0 || trialID(t) ==5)
        continue;
    end
    
    T1 = trialStart(t)+window(1)-win; T2 = trialStart(t)+window(2)+win;

    LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','VERBOSE',false); 
    LFPs = LFPs.streams.LFPs;
    Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu','VERBOSE',false);
    Snips = Snips.snips.eNeu;
    
    LFPs.data = bpfilt(LFPs.data',[10,50],LFPs.fs,3)';
    
    trig = Snips.ts(Snips.chan == chn & Snips.sortcode == code);
    trig = round((trig-T1)*fs)';
    
    trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
    trialinds(:,floor(trialinds(1,:))<=0) = [];
    trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];
    
    for ch = 1:96
        d = LFPs.data(ch,:);
        d = d(floor(trialinds));
        d = d - mean(d);
        if(size(d,2) == length(range))
            d = d';
        end
        stLFP{trialID(t)}{ch} = [stLFP{trialID(t)}{ch},d];
    end
    
end

%% plot
figure;
yl = [];
for ch = 1:96
    [c,r,e] = GetWadeChannelPosition(ch);
    subplot(10,10,(r-1)*10+c);
    
    for i = 1:length(stLFP)
        if(i==5)
            continue;
        end
        hold on;
        if(i == 9)
            plot(range/fs,mean(stLFP{i}{ch},2),'k')
        else
            plot(range/fs,mean(stLFP{i}{ch},2))
        end
    end
    
    xlim([range(1)/fs,range(end)/fs]);
    
    yl(ch,:) = ylim;
    if(ch==chn)
        title(num2str(ch),'fontsize',7,'Color','r');
    else
        title(num2str(ch),'fontsize',7);
    end
    axis off;
end

legend({'1','2','3','4','5','6','7','8'})


%% other test code
ind = Snips.chan == chn & Snips.sortcode == code;
snips = Snips.data(ind,:); sample = floor(linspace(1,size(snips,1), 100));
subplot(10,10,91);
plot(snips(sample,:)','Color',[0.7,0.7,0.7]); hold on; plot(mean(snips),'k','LineWidth',2);
ylim([min(min(snips(sample,:))),max(max(snips(sample,:)))]); 
xlim([1,size(snips,2)]);axis off;

subplot(10,10,100);
CrossCorr(Snips.ts(ind), 'ts2',Snips.ts(ind),'binsize', 0.001,'lag',[-window,window],'suppress_plot',0); axis off;

figure;
stimChn = 83;
for i = 1:length(stLFP)
    if(i==5)
        continue;
    end
    subplot(3,3,i)
    plot(range/fs,mean(stLFP{i}{stimChn},2))
end




test = [];
for i = 1:length(stLFP)
    if(i==5) 
        continue;
    end
    test = [test,stLFP{i}{stimChn}];
    
end



%% cell props
bad = []; window = 0.05;
for c = 1:96
    ind = Snips.chan == c & Snips.sortcode == 0;
    
    spkLim = 10000;
    
    if(sum(ind) < spkLim)
        bad(end+1) = c;
        continue;
    end
    
    ttl = sprintf('%s, %d s,Chn%d, %d Spks',blockname,round(T2-T1),c,sum(ind));
    disp(ttl); %#ok<DSPS>
    
    fig = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(fig,'visible','off');
    
    % spike waveform
    subplot(3,2,1);
    snips = Snips.data(ind,:); sample = floor(linspace(1,size(snips,1), 1000));
    plot(snips(sample,:)','Color',[0.7,0.7,0.7]); hold on; plot(mean(snips),'k','LineWidth',2);
    ylim([min(min(snips(sample,:))),max(max(snips(sample,:)))]);
    xlim([1,size(snips,2)]); axis off;
    title(ttl);
    
    % autocorrelation
    subplot(3,2,2);
    CrossCorr(Snips.ts(ind), 'ts2',Snips.ts(ind),'binsize', 0.002,'lag',[-window,window],'suppress_plot',0); axis off;
    title('Autocorrelation');
    
    % directional tuning
    subplot(3,2,3:6)
    spk = (Snips.ts(ind)' - T1)*fs;
    
    [direction,~] = DirectionalTuning(Mani,spk,fs,window,dt,1);
    
    edges = -180:6:180; edges = edges*pi/180;
    [N,~] = histcounts(direction,edges);
    polarplot(edges,[N,N(1)]);
    title('Directional Tuning');
    
    print('-painters',fig, '-dpsc2', packet, '-append');
    close(fig);
    
end

% % plot cross correlations
% fig = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(fig,'visible','off');
% disp([blockname,', Cross Correlations'])
% chns = 1:96; codes = zeros(1,96);
% Chns = chns; Chns(bad) = [];
% Codes = codes; Codes(bad) = [];
% sp = length(Chns)+1;
% for i = 1:length(Chns)
%     ind1 = Snips.chan == Chns(i) & Snips.sortcode == Codes(i);
%     
%     snips = Snips.data(ind1,:); sample = floor(linspace(1,size(snips,1), 100));
%     subaxis(sp, sp, 1, i+1, 'spacing', 0, 'padding', 0.001)
%     plot(snips(sample,:)','Color',[0.7,0.7,0.7]); hold on; plot(mean(snips),'k','LineWidth',2);
%     axis off; title([num2str(Chns(i)),',',num2str(Codes(i)),',',num2str(size(snips,1))],'fontsize',5)
%     subaxis(sp, sp, i+1, 1, 'spacing', 0, 'padding', 0.001)
%     plot(snips(sample,:)','Color',[0.7,0.7,0.7]); hold on; plot(mean(snips),'k','LineWidth',2);
%     axis off; title({[num2str(Chns(i)),',',num2str(Codes(i))],num2str(size(snips,1))},'fontsize',5)
%     
%     for j = i:length(Chns)
%         ind2 = Snips.chan == Chns(j) & Snips.sortcode == Codes(j);
%         
%         subaxis(sp, sp, i+1, j+1, 'spacing', 0, 'padding', 0.001)
%         window = 0.2;
%         bin = 0.002;
%         [cor,lags] = CrossCorr(Snips.ts(ind1), 'ts2',Snips.ts(ind2),'binsize', bin,'lag',[-window,window],'suppress_plot',0);
%         axis off;
%         if i~=j
%             ylim([min(cor),max(cor)])
%         end
%     end
% end
% 
% subaxis(sp,sp,1,1,'spacing', 0, 'padding', 0.001)
% str = sprintf('Win %dms\nBin %dms',round(window*1000),round(bin*1000));
% text(0,1,str,'HorizontalAlignment','left','VerticalAlignment','top','fontsize',7);
% axis off;
% 
% print('-painters',fig, '-dpsc2', packet, '-append');
% close(fig);

callps2pdf(packet);


