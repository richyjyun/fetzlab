% load

T1 = 10;
TT = TDT2mat('Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\Spanky-180802-093217',...
    'T1',T1,'T2',T1+180,'TYPE',4,'STORE','LFPs');
preLFPs = TT.streams.LFPs;

T1 = 23*60+10;
TT = TDT2mat('Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\Spanky-180802-093217',...
    'T1',T1,'T2',T1+180,'TYPE',4,'STORE','LFPs');
stimLFPs = TT.streams.LFPs;

% plot

R = corrcoef(preLFPs.data');

bad = [1,3,15,16,18,20,22,24,29,31,43,51,87,89,95];
good = 1:96; good(bad) = [];

time = 2; time = round(time*preLFPs.fs);
space = 5e-4;

figure; ax1 = subplot(2,2,1);
for i = good
    hold on; plot((1:time)/preLFPs.fs,preLFPs.data(i,1:time)-i*space);
end
yticks((-96:-1)*space); yticklabels(96:-1:1); ylabel('Channel'); 
xlim([0,time/preLFPs.fs]); xlabel('Time (s)');
title('Good LFP');

ax2 = subplot(2,2,2); 
for i = good
    hold on; plot((1:time)/preLFPs.fs,stimLFPs.data(i,1:time)-i*space);
end
yticks((-96:-1)*space); yticklabels(96:-1:1); ylabel('Channel'); 
xlim([0,time/preLFPs.fs]); xlabel('Time (s)');
title('Stimulator LFP');

linkaxes([ax1,ax2],'xy')


params.tapers = [5,9]; params.Fs = preLFPs.fs; params.fpass = [9,100]; params.trialave = 1;
movingwin = [0.25,0.01];

[S1,t1,f1] = mtspecgramc(preLFPs.data(good,:)',movingwin,params);
S1 = mean(S1);

[S2,t2,f2] = mtspecgramc(stimLFPs.data(good,:)',movingwin,params);
S2 = mean(S2);

subplot(2,2,3); plot(f1,S1); xlim([10,100]); xlabel('Frequency (Hz)'); ylabel('Power'); title('Aversge Power Spectra')
subplot(2,2,4); plot(f2,S2); xlim([10,100]); xlabel('Frequency (Hz)'); ylabel('Power'); title('Average Power Spectra')




