fname = 'Y:\~NeuroWest\Spanky\Neurochip\S_20180613_07\S_20180613_07.mat';
events = nc3events(fname);
stim = events.stim{1};

% [sta, x, n, names] = nc3sta(stim, [17:32], -20, 200, 20000, [0,0], fname);
% 
% save('Y:\~NeuroWest\Spanky\Neurochip\S_20180613_07\sta.mat','sta','x','n','names','-v7.3');
% 
% figure; plot(x,sta(:,10)); hold on; plot(x,sta(:,12));
% 
% figure; plot(x,sta)
% 
% [Data, end_sec] = nc3chan(0, 600, 20000, [], 20000, 'Y:\~NeuroWest\Spanky\Neurochip\S_20180613_07\S_20180613_07_Chan26.i16');


window = [-20,200]; % in ms
window = window/1000;
channels = 12;
fs = 20000;
filt = [0,0];

sweeps = [];
for s = 1:length(stim)
    if(stim(s)+window(1) <=0 || stim(s)+window(2) > events.session_time)
        continue;
    end
    
    [data, ~, ~] = nc3data(28, stim(s)+window(1), window(2)-window(1), 20000, filt, fname);
    
    sweeps(:,:,s) = data';
    
end

save('Y:\~NeuroWest\Spanky\Neurochip\S_20180613_07\sta_sweeps.mat','sweeps','window','channels','fs','-v7.3');

% mean subtract
sweeps_ms = zeros(size(sweeps));
for i = 1:size(sweeps,3)
    sweeps_ms(:,:,i) = sweeps(:,:,i) - mean(sweeps(:,:,i));
end

% artifact blanked. arbitrary windows for now
sweeps_b = sweeps_ms;
sweeps_b(:,400:420,:) = 0;

% do a moving average
binwidth = 600; %seconds
bins = stim(1):binwidth:stim(end);

bin = discretize(stim,bins);
bin(isnan(bin)) = 0;

movavg = [];
for i = 1:max(bin)
    movavg(i,:) = mean(sweeps_b(1,:,bin==i),3);
end

% plot
imagesc(movavg)

Accel = nc3data(38, 0, session_time-mod(session_time,1000), 100, [], fname);
Accel = interp(Accel,10);

% get moving average
a = movmean(Accel,100*1000);
thresh = 300; % arbitrary, always going to be the same with accelerometer

% find times below threshold
lowa = a <= thresh;
