function fig = stLFP(varargin)
% Summary
%   Plots the spike-triggered average of LFPs, triggering from the provided
%   channel/code pairs at the epochs specified.
%
% Inputs
%   tankpath    path of TDT tank
%   blockname   path of TDT block
%   trigChns    all channels to trigger off of
%   codes       corresponding sort codes (must have one for each channel)
%   stimChn     stimulated channel
%   times       time epochs to plot ([T1,T2] per row)
%   filt        frequency range for filtering (empty for no filter)
%   packet      print to packet or not (default true)
%   whiten      whiten data or not (default false)
%   window      +-time from spiking to look at, in seconds
%
% Outputs
%   fig         array of figures, returned if packet is empty
%
% RJY 07/16/2018

%% Set up name/value pairs

p = inputParser;

addParameter(p,'tankpath',[]);
addParameter(p,'blockname',[]);
addParameter(p,'trigChns',[]);
addParameter(p,'codes',[]);
addParameter(p,'times',[]);
addParameter(p,'filt',[]);
addParameter(p,'packet',true);
addParameter(p,'whiten',false);
addParameter(p,'window',0.05);
addParameter(p,'stimChn',[]);

parse(p,varargin{:});

tankpath = p.Results.tankpath;
blockname = p.Results.blockname;
trigChns = p.Results.trigChns;
codes = p.Results.codes;
times = p.Results.times;
filt = p.Results.filt;
packet = p.Results.packet;
whiten = p.Results.whiten;
window = p.Results.window;
stimChn = p.Results.stimChn;

%% Initialize variables

% figures
fig = zeros(length(trigChns),1);
for i = 1:length(trigChns)
    fig(i) = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12],'visible','off'); 
end

% plot colors for different times
colors = 'krbgy';

% packet path
date = blockname(8:13);
if(isempty(filt))
    f = 'raw';
else
    f = sprintf('filt%d-%d',filt(1),filt(2));
end
fname = fullfile('F:\S\Packets\stLFP', sprintf('%s_%s.ps',date,f));

%% Loop through each epoch, calculate, and plot
for t = 1:size(times,1)
    color = colors(t);
    
    % load LFPs
    T1 = times(t,1) - window;
    T2 = times(t,2) + window;
    LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs;
    
    % get all snippets from spike sorting
    Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu'); Snips = Snips.snips.eNe1;

    % filter
    if(~isempty(filt))
        LFPs.data = bpfilt(LFPs.data',filt,LFPs.fs,3)';
    end
    
%     LFPs.data = zscore(LFPs.data')';
    
    % define variables
    fs = LFPs.fs;
    range = round(-window*fs:1:window*fs);
    
    %% Loop through each trigger channel
    for i = 1:length(trigChns)
        
        set(0, 'CurrentFigure', fig(i));
        
        %% Get trigger times and define indices
        trig = (Snips.ts(Snips.chan == trigChns(i) & Snips.sortcode == codes(i))' - T1)*fs;
        
        trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
        trialinds(:,floor(trialinds(1,:))<=0) = [];
        trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];
        
        if(isempty(trialinds))
            continue;
        end
        
        %% Loop through all LFP channels and get stLFPs
        stLFPs = zeros(size(LFPs.data,1),length(range));
        for j = 1:size(LFPs.data,1)
            d = LFPs.data(j,:);
            d = d(floor(trialinds));
%             d = zscore(d);
            d = d - mean(d);
            d = mean(d,2);
            stLFPs(j,:) = zscore(d);
        end
        
        %% Whitening
        if(whiten)
            [E,D] = eig(cov(LFPs.data'));
            W = E*diag(diag(D).^(-1/2))*E';
            wstLFPs = W*stLFPs(good,:);
            
%             % maximum method
%             dif = max(wstLFPs')-min(wstLFPs');
%             ind = find(dif==max(dif));
            
            % Correlation method for matching y scale to compare 
            Corr = [];
            for j = 1:length(good)
                corr = corrcoef(stLFPs(good(j),:),wstLFPs(j,:));
                Corr(j) = corr(1,2);
            end
            ind = find(Corr == max(Corr));
            ind = find(good==19);
            ratio = (max(stLFPs(good(ind),:)) - min(stLFPs(good(ind),:))) /...
                (max(wstLFPs(ind,:)) - min(wstLFPs(ind,:)));
            wstLFPs = ratio*wstLFPs;
    
        else
            wstLFPs = stLFPs;
        end
 
        
        %% Plot
        yl(i,:,:) = nan(size(LFPs.data,1),2);
        for j = 1:size(LFPs.data,1)
            [c,r,~] = GetWadeChannelPosition(j);
            subplot(10,10,(r-1)*10+c);
            
            hold on;
            plot(range/fs,wstLFPs(j,:),color); hold on;
            yl(i,j,:) = ylim; xlim([-window,window]);
            
        end
                
    end
end

%% Cleaning up figures
for i = 1:length(trigChns)
    
    set(0, 'CurrentFigure', fig(i));

    %% Ylim, titles, axes
    YLIM = nanmedian(squeeze(yl(i,:,:)))*1.5;
    
    for j = 1:size(LFPs.data,1)
        
        [c,r,~] = GetWadeChannelPosition(j);
        subplot(10,10,(r-1)*10+c);
        ylim(YLIM);
        hold on;
        line([0 0], YLIM, 'linestyle', ':', 'color', [.5 .5 .5]);
        axis off
        sub_pos = get(gca,'position'); % get subplot axis position
        set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
        if(j == trigChns(i))
            title(num2str(j),'fontsize',7,'Color','r')
        elseif(j==stimChn)
            title(num2str(j),'fontsize',7,'Color','b')
        else
            title(num2str(j),'fontsize',7,'Color','k')
        end
        
    end
    
    %% Plot snips and autocorrelation (from last epoch, shouldn't matter)
    ind = Snips.chan == trigChns(i) & Snips.sortcode == codes(i);
    snips = Snips.data(ind,:); sample = floor(linspace(1,size(snips,1), 100));
    subplot(10,10,91);
    plot(snips(sample,:)','Color',[0.7,0.7,0.7]); hold on; plot(mean(snips),'k','LineWidth',2);
    ylim([min(min(snips(sample,:))),max(max(snips(sample,:)))]); axis off;
    sub_pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
    subplot(10,10,100);
    CrossCorr(Snips.ts(ind), 'ts2',Snips.ts(ind),'binsize', 0.001,'lag',[-window,window],'suppress_plot',0); axis off;
    sub_pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
    
    %% Title
    axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    txt = sprintf('%s; Trigger Channel %d, Code %d; %d Triggers, +-%dms window, Ylim %e, %e',...
        blockname,trigChns(i),codes(i),sum(ind),round(window*1000),YLIM(1),YLIM(2));
    text(0.5, 0.1,txt,...
        'HorizontalAlignment' ,'center','VerticalAlignment', 'top','fontsize',9)
    
    %% Print to packet or make visible
    if(packet)
        print('-painters',gcf, '-dpsc2', fname, '-append');
        close(gcf);
    else
        set(gcf,'visible','on');
    end
    
end

callps2pdf(fname);

