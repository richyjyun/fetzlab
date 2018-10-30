block = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\Spanky-180126-154107';

TT = TDT2mat(block,'TYPE',2);
Dscm = TT.epocs.Dscm;
[val,ind] = findpeaks(Dscm.data); 
val(end+1) = Dscm.data(end); ind(end+1) = length(Dscm.data);
ind = ind(val>1000); val = val(val>1000); 
tests = 2; 
times = [ind(tests)-val(tests),ind(tests)];%-val(tests)+val(2)];  % times = [ind(1) - val(1),ind(1)-val(1)+15000;ind(3)-val(3),ind(3)-val(3)+10000];
times(times==0) = 1;
times = Dscm.onset(times);

% stim = TDT2mat(block,'T1',times(1),'T2',times(2),'Type',5);
% stim = stim.scalars.Stim.ts;

cell = TDT2mat(block,'Type',3,'Store','eNe1','Channel',83);
cell = cell.snips.eNe1;

lfp = TDT2mat(block,'Type',4,'T1',times(1),'T2',times(2),'Store','LFPs','Channel',32);
lfp = lfp.streams.LFPs; fs = lfp.fs;

% [pks,loc] = findpeaks(double(lfp.data));
% stim = loc(pks > 3e-3)/fs + times(1);

% stim = cell.ts+0.005;

delay = zeros(1,length(stim));
spike = zeros(1,length(stim));
for i = 1:length(stim)
    closest = cell.ts(find(cell.ts>stim(i),1));
    lag = closest-stim(i);
    if(isempty(lag))
        continue;
    end
    delay(i) = lag;
end

figure; histogram(delay,1000);

spike(delay<=0.005) = 1;

temp = [0,spike,0];
d = diff(temp);
up = find(d==1);
down = find(d==-1);

runs = down-up;
runts = stim(up); 

figure; histogram(runs)

binwidth = 1;
bins = stim(1):binwidth:stim(end);

bin = discretize(runts,bins);
bin(isnan(bin)) = 0;

movavg = [];
for i = 1:max(bin)
    movavg(i) = mean(runs(bin==i));
end
figure; plot(movavg);



