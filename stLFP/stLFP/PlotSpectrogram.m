tankpath = 'Y:\~NeuroWest\Spanky\CycleTriggered-170710-143939\';
blockname = 'Spanky-171214-143604';

trigChn = 54; stimChn = 25; Code = 1;

window = 0.05;
T1 = times(1,1)-window; T2 = times(1,2)+window;
LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs;
Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1'); Snips = Snips.snips.eNe1;

fs = LFPs.fs;
range = round(-window*fs:1:window*fs);

allTrigs = (Snips.ts(Snips.chan == trigChn & Snips.sortcode == Code)' - T1)*fs;

bad = [trigChn];
good = 1:96; good(bad) = [];
% remove bad channels, and average it using those around
temp = LFPs.data;
temp = bpfilt(temp',[3,300],fs,3)';

[E,D] = eig(cov(temp(good,:)'));
W = E*diag(diag(D).^(-1/2))*E';

        
win = 1000; step = 200;
inds = 1:win;
Spect = [];
while inds(end) < length(allTrigs)
    trig = allTrigs(inds);
    % define all indices to get data from
    trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
    trialinds(:,floor(trialinds(1,:))<=0) = [];
    trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];

%     times = floor(trialinds(1,1):trialinds(end,end));
%     [E,D] = eig(cov(temp(good,times)'));
%     W = E*diag(diag(D).^(-1/2))*E';
 
    stLFPs = zeros(size(LFPs.data,1),length(range));
    % loop through all LFP channels and get stLFPs
    for j = 1:size(LFPs.data,1)
        d = temp(j,:);
        d = d(floor(trialinds));
        d = mean(d,2);
        stLFPs(j,:) = d;
    end
    wstLFPs = W*stLFPs(good,:);
                
    Spect(end+1,:) = wstLFPs(find(good == stimChn),:);
    
    
    inds = inds+step;
end

figure; imagesc(Spect); title('Post');





