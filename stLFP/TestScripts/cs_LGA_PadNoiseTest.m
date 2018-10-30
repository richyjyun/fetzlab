clear; close all;

%% Load in proper times
tankpath = 'Y:\~NeuroWest\Spanky\RandomStim-180314-124242\';
blockname = 'Spanky-180806-163308';
T1 = 23*60; T2 = T1+140; stimChn = 91;

LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','verbose',false);
LFPs = LFPs.streams.LFPs; fs = LFPs.fs;

figure; plot(LFPs.data(81,:));

params.tapers = [3,5]; params.Fs = fs; params.fpass = [1,100]; params.trialave = 1;
movingwin = [5,1];

figure; yl = [];
stimSpectrum = [];
for i = 1:96
    [pos,~] = GetLGAChannelPosition(i);
    subplot(19,20,pos);
    [S,t,f] = mtspecgramc(LFPs.data(i,:),movingwin,params);
    data = mean(S);
    stimSpectrum(i,:) = data;
    if(i==stimChn)
        plot(f,data/sum(data),'r');
        xlim([50,70])
    else
        plot(f,data/sum(data),'k');
        xlim([50,70])
    end
    axis off;
    yl(i,:) = ylim;
end

YLIM = [0,0.002];
for i = 1:96
    [pos,~] = GetLGAChannelPosition(i);
    subplot(19,20,pos);
        xlim([50,70])
        ylim(YLIM);
end

%% without stim
tankpath = 'Y:\~NeuroWest\Spanky\RandomStim-180314-124242\';
blockname = 'Spanky-180806-163308';
T1 = 10*60; T2 = T1+140; stimChn = 91;

LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','verbose',false);
LFPs = LFPs.streams.LFPs; fs = LFPs.fs;

figure; plot(LFPs.data(81,:));

params.tapers = [3,5]; params.Fs = fs; params.fpass = [1,100]; params.trialave = 1;
movingwin = [5,1];

figure; yl = [];
preSpectrum = [];
for i = 1:96
    [pos,~] = GetLGAChannelPosition(i);
    subplot(19,20,pos);
    [S,t,f] = mtspecgramc(LFPs.data(i,:),movingwin,params);
    data = mean(S);
    preSpectrum(i,:) = data;
    if(i==stimChn)
        plot(f,data/sum(data),'r');
        xlim([50,70])
    else
        plot(f,data/sum(data),'k');
        xlim([50,70])
    end
    axis off;
    yl(i,:) = ylim;
end

YLIM = [0,4e-10];
for i = 1:96
    [pos,~] = GetLGAChannelPosition(i);
    subplot(19,20,pos);
    xlim([1,100])
    %     ylim(YLIM);
end

figure;
ratio = stimSpectrum./preSpectrum;
yl = []; YLIM = [0,0.01];
for i = 1:96
    [pos,~] = GetLGAChannelPosition(i);
    subplot(19,20,pos);
    if(i==stimChn)
        plot(f,ratio(i,:)./sum(ratio(i,:)),'r');
        xlim([50,70])
    else
        plot(f,ratio(i,:)./sum(ratio(i,:)),'k');
        xlim([50,70])
    end
    axis off;
    yl(i,:) = ylim;
    ylim(YLIM);
end

%% irene data
TT = TDT2mat('Y:\~Neuro1\rec5');
LFPs = TT.streams.RawD;
filt = bpfilt(LFPs.data',[55,65],LFPs.fs,3)';

figure; yl = [];
preSpectrum = [];
for i = 1:96
    [pos,~] = GetLGAChannelPosition(i);
    subplot(19,20,pos);
    [S,t,f] = mtspecgramc(LFPs.data(i,:),movingwin,params);
    data = mean(S);
    preSpectrum(i,:) = data;
    if(i==stimChn)
        plot(f,data/sum(data),'r');
        xlim([50,70])
    else
        plot(f,data/sum(data),'k');
        xlim([50,70])
    end
    axis off;
    yl(i,:) = ylim;
end

figure;
for i = 1:96
    [pos,~] = GetLGAChannelPosition(i);
    subplot(19,20,pos);
    h = hilbert(filt(i,:));
    if(i == 66 || i ==68)
        bar(1,median(abs(h)),'r');
    else
        bar(1,median(abs(h)));
    end
    yl(i,:) = ylim;
end
for i = 1:96
    [pos,~] = GetLGAChannelPosition(i);
    subplot(19,20,pos);
    axis off;
    ylim([0,20])
end

%plot correlation
C = corrcoef(filt');
figure;
for i = 1:96
    [pos,~] = GetLGAChannelPosition(i);
    subplot(19,20,pos);
    c = C(i,70);
    if(i == 70)
        bar(1,c,'r');
    else
        bar(1,c);
    end
    yl(i,:) = ylim;
end
for i = 1:96
    [pos,~] = GetLGAChannelPosition(i);
    subplot(19,20,pos);
    axis off;
    ylim([min(yl(:,1)),max(yl(:,2))])
end


% plot coherence
params.tapers = [3,5]; params.Fs = fs; params.fpass = [1,100]; params.trialave = 1;
movingwin = [5,1];
figure;
for i = 1:96
    [pos,~] = GetLGAChannelPosition(i);
    subplot(19,20,pos);
    [C,phi,S12,S1,S2,f]=coherencyc(LFPs.data(i,:),LFPs.data(66,:),params);
    if(i == 70)
        plot(f,C,'r');
    else
        plot(f,C);
    end
    yl(i,:) = ylim;
end
for i = 1:96
    [pos,~] = GetLGAChannelPosition(i);
    subplot(19,20,pos);
    axis off;
    ylim([min(yl(:,1)),max(yl(:,2))])
    xlim([5,70]);
end








