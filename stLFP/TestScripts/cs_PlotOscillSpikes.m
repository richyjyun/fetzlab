tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';
blockname = 'Spanky-180122-130528';

TT = TDT2mat([tankpath,blockname],'TYPE',2);
Dscm = TT.epocs.Dscm;
[val,ind] = findpeaks(Dscm.data); 
val(end+1) = Dscm.data(end); ind(end+1) = length(Dscm.data);
ind = ind(val>1000); val = val(val>1000); 
tests = 1:2:length(val); 
times = [ind(tests)-val(tests),ind(tests)];%-val(tests)+val(2)];  % times = [ind(1) - val(1),ind(1)-val(1)+15000;ind(3)-val(3),ind(3)-val(3)+10000];
times(times==0) = 1;
times = Dscm.onset(times);

trigChn = 32;

T1 = times(1,1); T2 = times(1,2);
Snip = TDT2mat([tankpath,blockname(1,:)],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1','Channel',32);
Snip = Snip.snips.eNe1;
LFPs = TDT2mat([tankpath,blockname(1,:)],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','Channel',32);
LFPs = LFPs.streams.LFPs; fs = LFPs.fs;

LFP2 = TDT2mat([tankpath,blockname(1,:)],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','Channel',64);
LFP2 = LFP2.streams.LFPs; fs = LFP2.fs;

spks = Snip.ts(Snip.sortcode==1);
spks = round((spks-T1)*fs);

filt = bpfilt(LFPs.data',[15,25],LFPs.fs,3)';
h = hilbert(filt);
phase = angle(h);

spkPhase = phase(spks);
trough = spkPhase < -pi/2 | spkPhase > pi/2;

trig = spks(~trough)'; 
window = 0.1;
range = round(-window*fs:1:window*fs);

trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];

filt = bpfilt(LFPs.data',[15,25],LFPs.fs,3)';
d = filt(trialinds);
d = d-mean(d);
stLFP = mean(d,2);

figure; plot(range/fs,stLFP);

filt = bpfilt(LFP2.data',[15,25],LFPs.fs,3)';
d = filt(trialinds);
d = d-mean(d);
stLFP2 = mean(d,2);
figure; plot(range/fs,stLFP2);
