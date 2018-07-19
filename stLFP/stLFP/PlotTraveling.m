function PlotTraveling(tankpath,blockname,trigChns,Codes,times,packet)

% set(0, 'CurrentFigure', fig);

for t = 1:size(times,1)
    if(t == 1)
        color = 'k';
    elseif(t == 2)
        color = 'r';
    elseif(t == 3)
        color = 'b';
    elseif(t == 4)
        color = 'g';
    else
        color = 'y';
    end
    
    window = 0.05; %50 ms window, change as needed
    
    T1 = times(t,1) - window;
    T2 = times(t,2) + window;% get all LFPs
    LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs;
    % get all snippets from spike sorting
    Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1'); Snips = Snips.snips.eNe1;
    %     Discrim = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',2,'VERBOSE',0); Discrim = Discrim.epocs.Dscm;
    
    LFPs.data = bpfilt(LFPs.data',[15,50],LFPs.fs,3)';
    
    
    % define variables
    fs = LFPs.fs;
    range = round(-window*fs:1:window*fs);
    yl = zeros(length(trigChns),2);
    
    %     set(0, 'CurrentFigure', avg);
    for i = 1:length(trigChns)
        %         avg = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(gcf,'visible','off');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Getting all snippets from a single channel. Can use the PCA function
        % to apply sortcodes too and use Snips.sortcode to determine which ones
        % you want.
        
        trig = (Snips.ts(Snips.chan == trigChns(i) & Snips.sortcode == Codes(i))' - T1)*fs;
        %         trig = (Discrim.onset-T1)*fs; trig = trig';
        
        % define all indices to get data from
        trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
        trialinds(:,floor(trialinds(1,:))<=0) = [];
        trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];
        
        if(isempty(trialinds))
            continue;
        end
        
        %% Whitening
        %% Need to add another input and figure out bad channels for each day. this is for 20171208
     
        
        stLFPs = zeros(size(LFPs.data,1),length(range));
        % loop through all LFP channels and get stLFPs
        for j = 1:size(LFPs.data,1)
            d = LFPs.data(j,:);
            d = d(floor(trialinds));
            d = d - mean(d);
            d = mean(d,2);
            stLFPs(j,:) = d;
        end
        
        figure;
        yl = nan(96,2);
        % plot together
        for j = 1:size(LFPs.data,1)
            [c,r,e] = GetWadeChannelPosition(j);
            subplot(10,10,(r-1)*10+c);
            hold on;
            plot(range/fs,stLFPs(j,:),color); hold on;
            yl(j,:) = ylim; xlim([-window,window]);
            
        end
        
        YLIM = [nanmedian(yl)]*2;
        close(gcf)
        
        trough = [];
        bad = [18,20,52,31,51];
        for j = 1:size(LFPs.data,1)
            if( any(bad==j))
                continue;
            end
            [c,r,e] = GetWadeChannelPosition(j);
            ind = find(stLFPs(j,:) == min(stLFPs(j,:)));
            time = range(ind)/fs;
            trough(r,c) = time;
        end
        
        trough(trough==0) = nan;
        trough = trough*1000;
        figure;
        b = imagesc(trough);
        colormap(parula);
        set(b,'AlphaData',~isnan(trough))
        c = colorbar;
        ylabel(c, 'Trough Latency (ms)')
        
        if t == size(times,1)
            % set y lims and other parts
            for j = 1:size(LFPs.data,1)
                
                [c,r,e] = GetWadeChannelPosition(j);
                subplot(10,10,(r-1)*10+c);
                ylim(YLIM);
                hold on;
                line([0 0], YLIM, 'linestyle', ':', 'color', [.5 .5 .5]);
                axis off
                sub_pos = get(gca,'position'); % get subplot axis position
                set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
                if(j == trigChns(i))
                    title(num2str(j),'fontsize',7,'Color','r')
                else
                    title(num2str(j),'fontsize',7)
                end
            end
            
            % Plot snips and autocorrelation
            ind = Snips.chan == trigChns(i) & Snips.sortcode == Codes(i);
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
            
            % just titling. can include ylim (if normalizing) here to be able to
            % tell since axis is turned off.
            ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
            %             text(0.5,1,'Spike-triggered LFPs', 'HorizontalAlignment' ,'center','VerticalAlignment', 'top','fontsize',9);
            text(0.5, 0.1,[blockname,'; Trigger Channel ',num2str(trigChns(i)),'; ',num2str(length(trig)),' Triggers; +-',num2str(window*1000),'ms window, Ylim ',...
                num2str(YLIM(1)),', ',num2str(YLIM(2)) ],...
                'HorizontalAlignment' ,'center','VerticalAlignment', 'top','fontsize',9)
            
            % Print to packet
            %     if(exists('packet','var'))
            %             print('-painters',gcf, '-dpsc2', packet, '-append');
            %             close(gcf);
            %     end
        end
        
    end
end



