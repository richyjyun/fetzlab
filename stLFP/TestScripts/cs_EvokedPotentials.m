block = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\Spanky-181005-153115';

T1 = 480; T2 = 610;

% SUAs = TDT2mat(block,'T1',T1,'T2',T2,'Type',4,'Store','SUAs','Channel',1); 
% SUAs = SUAs.streams.SUAs;
% 
% [pks,loc] = findpeaks(double(SUAs.data));
% 
% bad = pks < 0.0015;
% pks = pks(~bad); loc = loc(~bad);
% 
% bad = find(diff(loc/SUAs.fs)<0.05);
% bad = bad+1;
% 
% pks(bad) = []; loc(bad) = [];

% figure; plot(SUAs.data); hold on; scatter(loc,zeros(1,length(loc)));

Stim = TDTbin2mat(block,'Type',5);
Stim = Stim.scalars.Test.ts;

LFPs = TDTbin2mat(block,'T1',T1,'T2',T2,'Type',4,'Store','LFPs');
LFPs = LFPs.streams.LFPs; 

trig = (Stim-T1)*LFPs.fs;
trig = round(trig);

window = [-0.01,0.1];
range = round(window(1)*LFPs.fs):1:round(window(2)*LFPs.fs);

trialinds = repmat(trig(1:end), length(range), 1) + repmat(range(:), 1, size(trig(1:end),2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];

stimChn = 57;

temp = medfilt1(LFPs.data(81,:),4);
noise = findNoise(temp,LFPs.fs);

% bad = 6.611e4:6.734e4;

% figure;
for i = 1:96
    
    
    [c,r,e] = GetWadeChannelPosition(i);
    subplot(10,10,(r-1)*10+c)
    
    hold on;
    
    d = LFPs.data(i,:); d = d-mean(d);
    d(noise) = nan;
    d = d(trialinds); hold on;
    plot(range/LFPs.fs,nanmean(d,2),'r');

    xlim(window);
    ylim([-5e-5,5e-5])

    if(i == stimChn)
        title(num2str(i),'color','k','fontsize',7);
    else
        title(num2str(i),'fontsize',7)
    end
    axis off;
end


for i = 1:96

    [c,r,e] = GetWadeChannelPosition(i);
    subplot(10,10,(r-1)*10+c)
     xlim([-0.01,0.1]);
    ylim([-5e-5,5e-5])
    
end

i = 88;
d = LFPs.data(i,:); d = d-mean(d);
d = d(trialinds);
figure; plot(mean(d(:,1:1000),2))
hold on; plot(mean(d(:,end-1000:end),2))


chns = 1:96;
bad = [1,3,16,18,20,22,15,25,54,21,29,26,31,28,51,30,95];
chns(bad) = [];

total = [];
for i = chns
    
    d = LFPs.data(i,:); 
    d(noise) = nan;
    d = d-nanmean(d);
    d = d(trialinds);
    
    avg = movmean(d,200,2,'omitnan');
    
    avg = avg-mean(avg(1:30,:));
    
%     figure; imagesc(avg');
%     caxis([-12e-5,2e-5])
    
    [trough,t] = min(avg(45:75,:));
    peak = max(avg(65:170,:));
    amp = peak-trough;
    total(end+1,:) = amp;
    
end

figure; plot(mean(zscore(total'),2));

y = mean(zscore(total'),2)'; x = 1:length(y);
dy = nanstd(zscore(total')');

figure; fill([x;flipud(x)],[y-dy;flipud(y+dy)],[0.7,.7,.7]);
hold on; line(x,y)

% compare to stlfp (or ctlfp)
Snips = TDT2mat(block,'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu','Channel',91,'verbose',false);
Snips = Snips.snips.eNe1; 

trig =round((Snips.ts(Snips.chan==91 & Snips.sortcode==2) - T1)*LFPs.fs);
trig = trig';

window = [-0.05,0.05];
range = round(window(1)*LFPs.fs):1:round(window(2)*LFPs.fs);

trialinds = repmat(trig(1:end), length(range), 1) + repmat(range(:), 1, size(trig(1:end),2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];

stLFP = medfilt1(LFPs.data,4);
stLFP = stLFP(trialinds);


avg = movmean(stLFP,500,2,'omitnan');
figure; imagesc(avg');
caxis([-1e-5,1.5e-5])



