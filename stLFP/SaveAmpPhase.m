function SaveAmpPhase(tankpath,blockname,trigChns,Codes,plotChns,times,rep)

%% Load
% 2. Epochs - Thd1 (threshold) and all corresponding values
% 3. Snips - Beta (snippets, 24kHz)
% 4. Streams - Mani (manipulandum, 3kHz), LFPs (3kHz), SUAs (24kHz), Filt
% (3kHz)
% 5. Scalars (all at 1Hz) - Trig (trigger channel), SUAc (SUA channels), Stim (Stim
% params and times)
% tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';


window = 0.2; %100 ms window, change as needed

T1 = times(1) - window;
T2 = times(2) + window;% get all LFPs
LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','CHANNEL',plotChns); LFPs = LFPs.streams.LFPs;
% get all snippets from spike sorting
Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1'); Snips = Snips.snips.eNe1;
%     Discrim = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',2,'VERBOSE',0); Discrim = Discrim.epocs.Dscm;


% define variables
fs = LFPs.fs;
range = round(-window*fs:1:window*fs);
% yl = zeros(length(trigChns),2);

%     for i = 1:length(plotChns)

trig = (Snips.ts(Snips.chan == trigChns & Snips.sortcode == Codes)' - T1)*fs;

% define all indices to get data from
trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];

if(isempty(trialinds))
    return;
end

sweeps = LFPs.data(floor(trialinds));

% Loops through frequency range in 2hz periods and get amplitude
% and phase
freq = 2:2:90; amp = []; phase = []; spikes = []; ratio = [];
prcnt = 0:10:100;
for f = 1:length(freq)-1
    filt = bpfilt(data',[freq(f),freq(f+1)],fs,3)';
    d = filt;
    d = d(floor(trialinds));
    h = hilbert(d);
    temp = abs(h);
    amp(f,:) = temp(ceil(size(temp,1)/2),:);
    temp = abs(hilbert(filt));
    for p = 1:length(prcnt)-1
        p1 = prctile(amp(f,:),prcnt(p)); p2 = prctile(amp(f,:),prcnt(p+1));
        % calculate firing probabilty (get left and right of regions using
        % diff, check if there is at least one spike within that)
        inds = temp > p1 & temp <= p2;
        Lbound = find(diff(inds) == 1); Rbound = find(diff(inds) == -1);
        if(inds(1) == 1), Lbound = [1,Lbound]; end
        if(inds(end) == 1), Rbound = [Rbound,length(inds)]; end
        spikes = 0;
        for s = 1:length(Lbound) 
            if(any(trig>Lbound(s) & trig < Rbound(s)))
                spikes = spikes + 1;
            end
        end
        ratio(f,p) = spikes/length(Lbound);
    end
    temp = angle(h);
    phase(f,:) = temp(ceil(size(temp,1)/2),:);
end

path = ['F:\S\Spectra\Post\',num2str(trigChns),'_',num2str(plotChns),'\'];
save([path,blockname,'_',num2str(rep),'.mat'],'sweeps','amp','phase','ratio');

end

