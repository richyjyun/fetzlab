function [stLFP,win,dt] = PlotstLFPTime(tankpath,blockname,trigChn,code,stimChn,times,win,dt)

% win of 30 and dt of 1 seems to work fairly well

window = 0.05; 

T1 = times(1) - window;
T2 = times(2) + window;% get all LFPs
LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','CHANNEL',stimChn,'VERBOSE',0); LFPs = LFPs.streams.LFPs;
Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu','VERBOSE',0); Snips = Snips.snips.eNeu;

LFPs.data = bpfilt(LFPs.data',[15,50],LFPs.fs,3)';

fs = LFPs.fs;
range = round(-window*fs:1:window*fs);

trig = (Snips.ts(Snips.chan == trigChn & Snips.sortcode == code)' - T1)*fs;
dt = round(dt*fs);
win = round(win*fs);

stLFP = [];
w1 = 1; w2 = w1+win-1;
while(w2 < max(trig))
    t = trig(trig>w1 & trig<w2);
    trialinds = repmat(t, length(range), 1) + repmat(range(:), 1, size(t,2));
    trialinds = round(trialinds);
    trialinds(:,trialinds(1,:)<=0) = [];
    trialinds(:,trialinds(end,:)>length(LFPs.data)) = [];
    
    d = LFPs.data(trialinds);
    temp = mean(d,2);
    stLFP(end+1,:) = temp - mean(temp); %mean subtracting in case not filtering
    w1 = w1+dt; w2 = w2+dt;
end

% figure;
% imagesc(stLFP);
% xticks([1:size(stLFP,2)/8:size(stLFP,2),size(stLFP,2)])
% xticklabels(-window*1000:window*1000/4:window*1000)
% hold on; yl = ylim; plot([size(stLFP,2)/2+1,size(stLFP,2)/2+1],yl,'k--');
% colorbar;

% figure; plot(mean(stLFP));

end