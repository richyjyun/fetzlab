function PlotWstLFP(tankpath,blockname,trigChns,Codes,times,filt,packet,whiten)

%% Load
% 2. Epochs - Thd1 (threshold) and all corresponding values
% 3. Snips - Beta (snippets, 24kHz)
% 4. Streams - Mani (manipulandum, 3kHz), LFPs (3kHz), SUAs (24kHz), Filt
% (3kHz)
% 5. Scalars (all at 1Hz) - Trig (trigger channel), SUAc (SUA channels), Stim (Stim
% params and times)
% tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';

for i = 1:length(trigChns)
    fig(i) = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); %set(gcf,'visible','off');
end
    
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
    Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNeu'); Snips = Snips.snips.eNe1;
%     Discrim = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',2,'VERBOSE',0); Discrim = Discrim.epocs.Dscm;

    if(~isempty(filt))
        LFPs.data = bpfilt(LFPs.data',filt,LFPs.fs,3)';
    end
    
    % define variables
    fs = LFPs.fs;
    range = round(-window*fs:1:window*fs);
%     yl = zeros(length(trigChns),2);
    
%     set(0, 'CurrentFigure', avg);
    for i = 1:length(trigChns)
        
        set(0, 'CurrentFigure', fig(i));

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

        temp = LFPs.data;
%         temp = bpfilt(temp',[3,300],fs,3)';
        
        stLFPs = zeros(size(LFPs.data,1),length(range));
        % loop through all LFP channels and get stLFPs
        for j = 1:size(LFPs.data,1)
            d = temp(j,:);
            d = d(floor(trialinds));
            d = d - mean(d);
            d = mean(d,2);
            stLFPs(j,:) = d;
        end
        
        % whiten
        if(whiten)
            [E,D] = eig(cov(temp(good,:)'));
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
        
        yl(i,:,:) = nan(96,2);
        % plot together
        for j = 1:size(LFPs.data,1)
            [c,r,e] = GetWadeChannelPosition(j);
            subplot(10,10,(r-1)*10+c);
            
            if(whiten && ismember(j,good))
                plot(range/fs,wstLFPs(find(good==j),:),color); hold on;
                yl(i,j,:) = ylim; xlim([-window,window]);
            elseif ~whiten
                hold on;
                plot(range/fs,stLFPs(j,:),color); hold on;
                yl(i,j,:) = ylim; xlim([-window,window]);
            end
             
%             % plot st and wst
%             plot(range/fs,stLFPs(j,:),'k'); hold on;
%             if(whiten && ismember(j,good))
%                 plot(range/fs,wstLFPs(find(good==j),:),'g');
%                 yl(j,:) = ylim;
%             end
            
        end
        

        
    end
end


%% Cleaning up figures
for i = 1:length(trigChns)
    
    set(0, 'CurrentFigure', fig(i));

    %         YLIM = [min(yl(:,1)),max(yl(:,2))];
    YLIM = [nanmedian(squeeze(yl(i,:,:)))]*1.5;
    
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
        txt = sprintf('%s; Trigger Channel %d, Code %d; %d Triggers, +-%dms window, Ylim %e, %e',...
            blockname,trigChns(i),Codes(i),sum(ind),round(window*1000),YLIM(1),YLIM(2));
        text(0.5, 0.1,txt,...
            'HorizontalAlignment' ,'center','VerticalAlignment', 'top','fontsize',9)
        
        % Print to packet or make visible
        if(exists('packet','var'))
            print('-painters',gcf, '-dpsc2', packet, '-append');
            close(gcf);
        else
            set(gcf,'visible','on');
        end
    end
end

