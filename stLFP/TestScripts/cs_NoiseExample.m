block = 'Y:\~NeuroWest\Spanky\RandomStim-180314-124242\Spanky-180914-145727';

TT = TDT2mat(block,'T1',15*60,'T2',20*60);

Data = TT.streams.LFPs.data;

fs = TT.streams.LFPs.fs;

T1 = 155; T2 = 185;

Sample = Data(:,round(T1*fs):round(T2*fs));

figure; plot(Sample(93,:))

Data = Data(:,round(T1*fs):round(T2*fs));

save('F:\S\NoiseSample\NoStim.mat','Data','fs','-v7.3');






