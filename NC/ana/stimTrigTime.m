function stimTrigTime(fpath,win,dt)
% Summary
%   Plots and saves stimulus triggered event histograms as a moving sum
%   over time. Run stimTrigHist first.
%
% Inputs
%   fpath - path of the experiment
%   win - +-time to look at (s)
%   dt - time steps for moving sum (s)
%
% RJY 07/09/2018

% default values
if(~exist('win'))
    win = 100;
end
if(~exist('dt'))
    dt = 1;
end

% load data
load(fullfile(fpath,'StimTrigHist'));

start = 0;

movBin = {};

% loop through and load into movBins
while start+win <= max(trig)
    
    inds = find(trig >= start & trig < start+win);
    
    if(isempty(inds))
        
        start = start+dt;
        
        continue;
        
    end
    
    for i=1:size(Hist,2)
        
        if( length(inds) == 1)
            
            if( i==1 )
                
                movBin{end+1,i} = cell2mat(Hist(inds,i));
                
            else
                
                movBin{end,i} = cell2mat(Hist(inds,i));
                
            end
            
        else
            
            if( i==1 )
                
                movBin{end+1,i} = sum(cell2mat(Hist(inds,i)));
                
            else
                
                movBin{end,i} = sum(cell2mat(Hist(inds,i)));
                
            end
            
        end
        
    end
    
    start = start+dt;
    
end

% print to packet
packet = fullfile(fpath,'StimTrigOverTime.ps');

for i = 1:size(movBin,2)
    
    fig = figure('visible','off');
    
    temp = cell2mat(movBin(1:18500,i));
    
    imagesc(temp(:,round(size(temp,2)/2):size(temp,2)));
    
    title({['Stim Triggered Event ',num2str(i)],...
        ['Window ',num2str(win),'s, dt ',num2str(dt),'s']});
    
    ylabel(['Experiment Time (',num2str(dt),'s)']);
    
    xlabel('Time Since Stim (ms)');
    
    c = colorbar;
    
    ylabel(c,'Event Counts');
    
    print('-painters','-fillpage',fig, '-dpsc2', packet, '-append');
    
    close(fig);
    
end

callps2pdf(packet);

end

