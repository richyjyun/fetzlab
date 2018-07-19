function [Pre, Post] = getCycleTrigLFP(tankpath,blockname,stimChns,times)

for t = 1:size(times,1)
    
    T1 = times(t,1); %- window;
    T2 = times(t,2); %+ window;% get all LFPs
    LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs;
    % get all snippets from spike sorting
    Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1'); Snips = Snips.snips.eNe1;
    
    % define variables
    fs = LFPs.fs;
    
    % get cycle trigger
    bp = [15,25];
    stim = bpfilt(LFPs.data(stimChns,:)',bp,fs,3)'; %filter one channel
    testtime = round(2*160*fs); % 2 minutes of test
    thresh = std(stim(1:testtime)); % set threshold to be a standard deviation
    
    % find threshold crossings
    trig = stim<=-thresh; trig = diff(trig); trig(trig<0) = 0;
    trig = find(trig);
    
    % window discrimination
    period = 1/mean(bp);
    delay1 = round((period/5)*fs); % 10ms delay
    delay2 = round((period/5+period/2)*fs); % 25ms delay after delay1
    trig((trig+delay2)>length(stim)) = [];
    trig = trig((stim(trig+delay1) <= -thresh) & (stim(trig+delay2) >= thresh));
    
    if(isempty(trig))
        Pre = []; Post = [];
        return;
    end
    
    % save oscillations
    range = round(-0.05*fs:1:0.05*fs);
    trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
    trialinds(:,floor(trialinds(1,:))<=0) = [];
    trialinds(:,floor(trialinds(end,:))>length(stim)) = [];
    
    Oscill = cell(1,96);
    for chn = 1:96
        temp = bpfilt(LFPs.data(chn,:)',[15,25],fs,3)'; %filter one channel
        Oscill{chn} = temp(floor(trialinds));
    end
    
    if(t == 1)
        Pre = Oscill;
    else
        Post = Oscill;
    end
end

if(~exist('Post'))
    Post = [];
end

end




