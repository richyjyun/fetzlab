function PlotAmpPhase(tankpath,blockname,trigChn,stimChn,times)

fig = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); %set(gcf,'visible','off');
Theta = cell(1,2);
PLV = cell(1,2);

for t = 1:2
    
    window = 0.2;
    
    T1 = times(t,1) - window;
    T2 = times(t,2) + window;% get all LFPs
    LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','CHANNEL',stimChn); LFPs = LFPs.streams.LFPs;
    %     Dscm = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',2,'VERBOSE',0); Dscm = Dscm.epocs.Dscm;
    Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1'); Snips = Snips.snips.eNe1;
    
    fs = LFPs.fs;
    range = round(-window*fs:1:window*fs);
    
    %     trig = (Dscm.onset'-T1)*fs;
    trig = (Snips.ts(Snips.chan == trigChn & Snips.sortcode == 1)' - T1)*fs;
    
    trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
    trialinds(:,floor(trialinds(1,:))<=0) = [];
    trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];
    
    [freq,amp,phase] = getAmpPhase(LFPs.data,trialinds,fs,2,90,2);
    [Theta{t},PLV{t}] = getPercentiles(amp,phase,0:10:100);
    
end

T = Theta{2}-Theta{1};
bad = abs(T)>pi;
T(bad) = 2*pi-abs(T(bad));

P = PLV{2};%-PLV{1};

map = hsv;

M = max(max(abs(T)));
minv = -M; maxv = M;

ncol = size(map,1);
scaled = round(1+(ncol-1)*(T'-minv)/(maxv-minv));
RGB = ind2rgb(scaled,map);
HSV = rgb2hsv(RGB);
minv = min(P(:)); maxv = max(P(:));
newVal = (P'-minv)/(maxv-minv);
HSV(:,:,2) = newVal;

imagesc(hsv2rgb(HSV));
colormap hsv; set(gca,'YDir','normal'); c = colorbar;
xticks(1.5:1:10.5); xticklabels(10:10:100); xlabel('Amplitude Envelope (Percentile)');
yticks(freq-0.5); yticklabels(freq*2); ylabel('Frequency (Hz)');
set(c,'Ticks',0:1/8:1,'TickLabels',{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'});

end



