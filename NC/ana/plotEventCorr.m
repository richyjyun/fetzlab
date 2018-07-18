function plotEventCorr(fpath,fname,fig)
% Summary
%   Plots auto- and cross-correlations of all events in NC data.
%
% Inputs
%   fpath - path of the experiment
%   fname - path of the Matlab settings file
%
% RJY 07/05/2018

% setup of figure
if(~exist('fig'))
    fig = figure;
else
    set(0, 'CurrentFigure', fig)
end

% get all events
Events = nc3events(fname);

% get all events that were on
inds = find(Events.ndiscrim ~= 0); 

% loop through and plot 
n = length(inds);
for i = 1:length(inds)
    for j = 1:length(inds)
        subplot(n,n,(i-1)*n+j)
        CrossCorr(Events.discrim{inds(i)}, 'ts2',Events.discrim{inds(j)},'binsize', 0.002,'lag',[-0.2,0.2],'suppress_plot',0);
        if(i==1)
            title({['Event ',num2str(j)],[num2str(length(Events.discrim{inds(j)})),' Occurances']});
        end
    end
end

orient(fig,'landscape')

packet = fullfile(fpath,'EventCorrelation.ps');
print('-painters','-fillpage',fig, '-dpsc2', packet, '-append');
close(fig);
callps2pdf(packet);

end