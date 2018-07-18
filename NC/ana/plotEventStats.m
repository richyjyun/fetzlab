function plotEventStats(fpath,fname)
% Summary
%   Saves event characteristics including ISI, autocorrelation, and rate
%
% Inputs
%   fpath - path of the experiment
%
% RJY 07/02/2018

% initilialize and load data
packet = fullfile(fpath,'EventStats.ps');

load(fullfile(fpath,'Sleep'));

Events = nc3events(fname);
EventChns = find(Events.ndiscrim>0);

% loop through all event channels
for chn = 1:length(EventChns)
    
    figure('visible','off');
    
    % set event times
    spkTime = Events.discrim{EventChns(chn)};
    
    % plot autocorrelation
    subplot(3,1,1);
    CrossCorr(spkTime, 'ts2',spkTime,'binsize', 0.001,'lag',[-0.1,0.1],'suppress_plot',0);
    title({['Event ',num2str(EventChns(chn))],'Autocorrelation'}); 
    yticks([]); xlabel('Time (ms)');
    
    % plot spikerate across experiment
    subplot(3,1,[2,3]);
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
    
    % print packet (can't do vector for sleep overlay)
    print('-opengl','-fillpage',gcf, '-dpsc2', packet, '-append');
    close(gcf);
    
end

callps2pdf(packet);

end