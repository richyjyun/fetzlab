% load (NeuroWest)
% 
% T1 = 10;
% TT = TDT2mat('Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\Spanky-180802-093217',...
%     'T1',T1,'T2',T1+180,'TYPE',4,'STORE','LFPs');
% preLFPs = TT.streams.LFPs;
% 
% T1 = 23*60+10;
% TT = TDT2mat('Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\Spanky-180802-093217',...
%     'T1',T1,'T2',T1+180,'TYPE',4,'STORE','LFPs');
% stimLFPs = TT.streams.LFPs;

% load (Irene)

T1 = 10;
TT = TDT2mat('Y:\~Neuro1\transfer\Gomez-180111\rec8',...
    'T1',T1,'T2',T1+180,'TYPE',4,'STORE','RawD');
preLFPs = TT.streams.RawD;


T1 = 10;
TT = TDT2mat('Y:\~Neuro1\transfer\Gomez-180111\rec7',...
    'T1',T1,'T2',T1+120,'TYPE',4,'STORE','RawD');
stimLFPs = TT.streams.RawD;


% plot
% bad = [1,3,15,16,18,20,22,24,29,31,43,51,87,89,95];

R = corrcoef(stimLFPs.data');
bad = find(mean(R) < 0.6);

R = corrcoef(preLFPs.data');
bad = union(bad,find(mean(R) < 0.6));

good = 1:96; good(bad) = [];

% good = '17/26/30/32/39/49/53/55/57/58/61/63/69/71/74/77/79/81/83/84/88/90/91/93/94';
% good = strsplit(good,'/'); good = cellfun(@str2num,good);

time = 2; time = round(time*preLFPs.fs);
% space = 5e-4;
space = 500;

figure; ax1 = subplot(2,2,1);
for i = good
    hold on; plot((1:time)/preLFPs.fs,preLFPs.data(i,(time+1):2*time)-i*space);
end
yticks((-96:-1)*space); yticklabels(96:-1:1); ylabel('Channel'); 
xlim([0,time/preLFPs.fs]); xlabel('Time (s)');
title('Good LFP');

ax2 = subplot(2,2,2); 
for i = good
    hold on; plot((1:time)/preLFPs.fs,stimLFPs.data(i,(time+1):2*time)-i*space);
end
yticks((-96:-1)*space); yticklabels(96:-1:1); ylabel('Channel'); 
xlim([0,time/preLFPs.fs]); xlabel('Time (s)');
title('Stimulator LFP');

linkaxes([ax1,ax2],'xy')


params.tapers = [3,5]; params.Fs = preLFPs.fs; params.fpass = [10,100]; params.trialave = 1;
movingwin = [5,1];

[S1,t1,f1] = mtspecgramc(preLFPs.data(good,:)',movingwin,params);
S1 = mean(S1);

[S2,t2,f2] = mtspecgramc(stimLFPs.data(good,:)',movingwin,params);
S2 = mean(S2);

subplot(2,2,3); plot(f1,S1,'k','linewidth',2); 
xlim([10,100]); xlabel('Frequency (Hz)'); ylabel('Power'); title('Aversge Power Spectra')
subplot(2,2,4); plot(f2,S2,'k','linewidth',2);  
xlim([10,100]); xlabel('Frequency (Hz)'); ylabel('Power'); title('Average Power Spectra')




