function stimTrigHist(fpath,fname,binwidth,window,fig)
% Summary
%   Plots and saves stimulus triggered event histograms
%
% Inputs
%   fpath - path of the experiment
%   fname - path of the Matlab settings file
%   binwidth - width of bins for histogram (s)
%   window - +-time to look at (s)
%   fig - figure
%
% RJY 07/06/2018

% default values
if(~exist(binwidth))
    binwidth = 0.001;
end
if(~exist(window))
    window = 0.1;
end

% setup of figure
if(~exist('fig'))
    fig = figure;
else
    set(0, 'CurrentFigure', fig)
end

% get all events
Events = nc3events(fname);

% get stimulation event. ASSUMES ONE STIMULUS
stim = find(Events.ndiscrim == Events.nstim(1));

% get all events that were on
inds = find(Events.ndiscrim ~= 0);
inds(inds==stim) = []; % dont keep stim in this list

% set triggeres to be when stimulation occured
trig = Events.stim{stim};

% set up the bin edges 
binedges = -window:binwidth:window;
edges = repmat(trig,1,length(binedges)) + repmat(binedges,length(trig),1);
edges = mat2cell(edges,ones(1,size(edges,1)),size(edges,2));

% loop through the events and get all histogram counts
Hist = {};
for i = 1:length(inds)
    [N,~] = cellfun(@(x) histcounts(Events.discrim{inds(i)},x),edges,'Uni',0);
    subplot(length(inds),1,i);
    temp = cell2mat(N); temp = sum(temp); temp = temp/sum(temp);
    bar((binedges(1:end-1)+diff(binedges)/2)*1000,temp,1,'facecolor','k','facealpha',1,'edgealpha',0);
    title(['Event ',num2str(inds(i))]); ylabel('Probability of Spiking');
    Hist(:,i) = N;
end
xlabel('Time from Stim (ms)');

% save data for later
save(fullfile(fpath,'StimTrigHist'),'trig','binwidth','binedges','window','Hist','-v7.3');

% print figure to packet
packet = fullfile(fpath,'StimTrigHistogram.ps');
print('-painters','-fillpage',fig, '-dpsc2', packet, '-append');
close(fig);
callps2pdf(packet);

end

