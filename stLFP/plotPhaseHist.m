function plotPhaseHist(tankpath,blockname,trigChns,Codes,times,packet)

%% Load
% 2. Epochs - Thd1 (threshold) and all corresponding values
% 3. Snips - Beta (snippets, 24kHz)
% 4. Streams - Mani (manipulandum, 3kHz), LFPs (3kHz), SUAs (24kHz), Filt
% (3kHz)
% 5. Scalars (all at 1Hz) - Trig (trigger channel), SUAc (SUA channels), Stim (Stim
% params and times)
% tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';
avg = figure;%('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); 
set(gcf,'visible','off');
yl = nan(96,2);

for t = 1:size(times,1)
    if(t == 1)
        color = 'b';
    elseif(t == 2)
        color = 'o';
    elseif(t == 3)
        color = 'g';
    elseif(t == 4)
        color = 'r';
    else
        color = 'y';
    end
    
    %     window = 0.05; %50 ms window, change as needed
    
    T1 = times(t,1); %- window;
    T2 = times(t,2); %+ window;% get all LFPs
    LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs;
    % get all snippets from spike sorting
    Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1'); Snips = Snips.snips.eNe1;
    
    % define variables
    fs = LFPs.fs;
    
    % filter
    Beta = bpfilt(LFPs.data',[15,25],fs,3)';
    
    set(0, 'CurrentFigure', avg);
    
    % get timestamps
    trig = (Snips.ts(Snips.chan == trigChns & Snips.sortcode == Codes)-T1)*fs;
    trig = round(trig);
    
    if(isempty(trig))
        return;
    end
    
    bins = -pi:pi/9:pi;
    % plot
    for j = 1:size(Beta,1)
        [c,r,e] = GetWadeChannelPosition(j);
        subplot(10,10,(r-1)*10+c);
        
        h = hilbert(Beta(j,:));
        ang = angle(h);
        h = histogram(ang(trig),bins,'EdgeAlpha',0); hold on;
        yl(j,:) = ([min([yl(j,1),h.BinCounts]),max([yl(j,2),h.BinCounts])]);
    end
end

%         YLIM = [min(yl(:,1)),max(yl(:,2))];
YLIM = [nanmedian(yl)]*1.5;

if t == size(times,1)
    % set y lims and other parts
    for j = 1:size(LFPs.data,1)
        
        [c,r,e] = GetWadeChannelPosition(j);
        subplot(10,10,(r-1)*10+c);
        ylim(yl(j,:));
        xlim([-pi,pi])
        hold on;
%         line([0 0], YLIM, 'linestyle', ':', 'color', [.5 .5 .5]);
        axis off
%         sub_pos = get(gca,'position'); % get subplot axis position
%         set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
        if(j == trigChns)
            title(num2str(j),'fontsize',7,'Color','r')
        else
            title(num2str(j),'fontsize',7)
        end
    end
    
    % just titling. can include ylim (if normalizing) here to be able to
    % tell since axis is turned off.
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    %             text(0.5,1,'Spike-triggered LFPs', 'HorizontalAlignment' ,'center','VerticalAlignment', 'top','fontsize',9);
    text(0.5, 0.1,[blockname,'; Trigger Channel ',num2str(trigChns)],...
        'HorizontalAlignment' ,'center','VerticalAlignment', 'top','fontsize',9)
    
    % Print to packet
    %     if(exists('packet','var'))
    %     print('-painters',gcf, '-dpsc2', packet, '-append');
    set(avg,'PaperOrientation','landscape');
    print(gcf, '-dpsc2', packet, '-append','-fillpage');
    %     close(gcf);
    %     end
end

end




