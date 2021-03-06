function plotSpikeStats(fpath)
% Summary
%   Saves spike characteristics including waveform, ISI, autocorrelation,
%   and spikerate
%
% Inputs
%   fpath - path of the experiment
%
% RJY 07/02/2018

packet = fullfile(fpath,'SpikeStats.ps');

load(fullfile(fpath,'SpkChannels'));
load(fullfile(fpath,'Sleep'));

fs = 20000; % may need to change this);

for chn = 2%1:length(SpkChannels)
    
    % load spike data
    file = fullfile(fpath,['spk_',num2str(SpkChannels(chn)),'.mat']);
    if(~exist(file))
        warning([file,' does not exist']);
        continue;
    end
    load(file);
    
    % plot spike waveform
    figure('visible','off');
    subplot(3,1,1);
    sample = Spike(:,1:1000:size(Spike,2));
    x = (0:1:(size(sample,1)-1))/fs*1000;
    plot(sample,'color',[0.5,0.5,0.5]);
    hold on; plot(mean(sample,2),'k','linewidth',2);
    title(['Channel ',num2str(SpkChannels(chn)),' Spike Waveform']);
    xlim([0,size(Spike,1)]);
    xlabel('Time (ms)'); ylabel('Amplitude (\muV)');
    
    %     % plot ISI histogram
    %     subplot(4,1,2);
    %     ISI = diff(spkTime)*1000;
    %     edges = min(ISI):1:max(ISI);
    %     histogram(ISI,edges,'Facecolor','k','facealpha',1,'edgealpha',0);
    %     xlim([0,100]); xlabel('ISI (ms)');
    %     ylabel('Counts'); title('ISI histogram');
    
    % plot autocorrelation
    subplot(3,1,2);
    CrossCorr(spkTime, 'ts2',spkTime,'binsize', 0.001,'lag',[-0.1,0.1],'suppress_plot',0);
    title('Auotcorrelation'); yticks([]); xlabel('Time (ms)');
    
    % plot spikerate across experiment
    subplot(3,1,3);
    binwidth = 60; % in seconds
    binedges = 0:binwidth:spkTime(end);
    histogram(spkTime,binedges,'facecolor','k','edgecolor','k');
    
    % overlay sleep times
    start = find(diff(sleep) == 1)/sleepfs;
    finish = find(diff(sleep) == -1)/sleepfs;
    yl = ylim;
    for i = 1:length(start)
        hold on;
        area([start(i),finish(i)],[yl(2),yl(2)],...
            'facecolor','r','facealpha',0.3,'edgealpha',0);
    end
    
    title('Spike rate with sleep times'); xlim([0,binedges(end)]);
    ylabel(['Spikes per ',num2str(binwidth),'s'])
    hours = 0:3600:binedges(end);
    xticks(hours); xticklabels(0:(length(hours)-1))
    xlabel('Experiment Time (h)');
    
    % print packet
    print('-opengl','-fillpage',gcf, '-dpsc2', packet, '-append');
    close(gcf);
    
end

callps2pdf(packet);

end