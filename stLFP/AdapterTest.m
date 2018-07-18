function [allChns,avgDist] = AdapterTest(blockname,trigChns,stimChn,times)

%clear; %close all;
%% Load
% 2. Epochs - Thd1 (threshold) and all corresponding values, Dscm -
% Discrimination times
% 3. Snips - Beta (snippets, 24kHz)
% 4. Streams - Mani (manipulandum, 3kHz), LFPs (3kHz), SUAs (24kHz), Filt
% (3kHz)
% 5. Scalars (all at 1Hz) - Trig (trigger channel), SUAc (SUA channels), Stim (Stim
% params and times)
tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';
% blockname = 'Spanky-170831-145709';

%% Define regions in seconds
DSCM = 1;  % trigger off of discrimination that's been saved (1), or find spikes using eNe1 (0)
Whiten = 0;   % whether to perform whitening or not
% trigChns = [56];%[20,44,70];%56;%[20,44,56,70]; % channel to trigger off of
dist = []; [target(1),target(2),~] = GetWadeChannelPosition(trigChns(1));
% stimChn = 5;
% times = [1,18;30,47];%;97,115];%;48,62;75,90];%;37,56;78,97];%[21,41;57,75;91,110];
badChns = [];

%% Load in each and plot spike triggered LFP
window = 0.05; back = 0.05; % in seconds
yl = zeros(96,2); % saving all the ylims
ap = zeros(length(trigChns),size(times,1)); % number of spikes in each phase
amp = [];
phase = [];
avgDist = figure;%('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); %set(avgDist,'visible','off');

% if(DSCM)
    allChns = figure;%('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); %set(allChns,'visible','off');
% end

bad = [];

for i = 1:size(times,1)
    if(i == 1)
        color = 'k';
    elseif(i == 2)
        color = 'b--';
    elseif(i == 3)
        color = 'r';
    else
        color = 'g';
    end
    T1 = times(i,1) - back;
    T2 = times(i,2) + window;
    LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs','VERBOSE',0); 
    if(isempty(LFPs.streams))
        title(['COULD NOT ACCESS FILES FOR ',blockname]);
        return;
    end
    LFPs = LFPs.streams.LFPs;
    fs = LFPs.fs; range = round(-back*fs:1:window*fs); %50 ms window for capturing CCEPs
    
    if(~DSCM)
        D = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'VERBOSE',0);
    else
        D = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',2,'VERBOSE',0);
    end
    
    for chn = 1:length(trigChns)
        
        trigChn = trigChns(chn);
        
        if(~DSCM)
            Discrim = D.snips.eNe1.ts(D.snips.eNe1.chan==trigChn & D.snips.eNe1.sortcode == 1);
%             allChns = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); %set(allChns,'visible','off');
        else
            Discrim = D.epocs.Dscm.onset;
        end
        
        trig = ((Discrim - T1)*fs)'; ap(chn,i) = length(trig);
%         trig = trig(1:5000);
%         
%         if length(trig)>15000
%             trig = trig(1:15000);
%         end
        
        trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
        trialinds(:,floor(trialinds(1,:))<=0) = [];
        trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];
        
%         bpf = [15,300]; % bandpass filter (Hz)
%         [b1,a1] = butter(2,bpf(2)/(fs/2)); % 5th order lowpass
%         [b2,a2] = butter(2,bpf(1)/(fs/2),'high'); % 1st order highpass
        %       figure; freqz(b1,a1);
        
        
        stLFPs = zeros(size(LFPs.data,1),length(range));
        
        for j = 1:size(LFPs.data,1)
            d = LFPs.data(j,:);
            
%             d = filtfilt(b1,a1,double(d)); % zero phase filtering
%             d = filtfilt(b2,a2,double(d));
%             
            d = d(floor(trialinds));
            

            d = d - mean(d);
            
            % remove those with too large of std, or zero std
            bad = std(d) > 1e-3 | std(d)==0;
            d(:,bad) = [];
            
            d = mean(d,2);
            
            stLFPs(j,:) = d;
            
            %             %probably a bad channel, remove from whitening matrix
            %             if(max(abs(d))<=0.5e-6)
            %                 W(j,:) = 0; W(:,j) = 0;
            %                 badChns(end+1) = j;
            %             end
            
            % removing trigger channel
%             if(j == trigChn)
%                 W(j,1:j-1) = 0; W(j,j+1:end) = 0;
%                 W(1:j-1,j) = 0; W(j+1:end,j) = 0;
%             end
            
        end
        
        if(Whiten)
            W = sqrtm(inv(cov(LFPs.data')));
            wstLFPs = W*stLFPs;
        else
            wstLFPs = stLFPs;
        end
        
        set(0, 'CurrentFigure', allChns);
        for j =1:size(LFPs.data,1)
            [c,r,e] = GetWadeChannelPosition(j);
            dist(j) = abs(c-target(1))+abs(r-target(2));
            subplot(10,10,(r-1)*10+c); hold on;
            plot(range/fs,wstLFPs(j,:),color);
            xlim([-back,window]);
            yl(j,:) = ylim;
            axis off;
            title(num2str(j),'fontsize',7)
        end
        
        if(i == size(times,1))
            ylNorm = nanmedian(yl,1).*1.5;
            
            for j = 1:size(LFPs.data,1)
                [c,r,e] = GetWadeChannelPosition(j);
                subplot(10,10,(r-1)*10+c);  hold on;
                ylim(ylNorm);
                line([0 0], ylNorm, 'linestyle', ':', 'color', [.7 .7 .7]);
                axis off;
%                 sub_pos = get(gca,'position'); % get subplot axis position
%                 set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
                %             plot([0,0],[-ylNorm,ylNorm],'k--');
            end
            
            % Plot averages of manhattan distance from trig. channel
            set(0, 'CurrentFigure', avgDist);
            checkDist = max(4,dist(stimChn));
            for k = 1:checkDist
                ind = find(dist == k);
                subplot(checkDist,1,k); hold on;
                plot(range/fs,mean(wstLFPs(ind,:)),color);
                xlim([-back,window]);
                ylim(ylNorm)
                title([num2str(k),' chns away from Trig'])
            end
        end
        
%         if(~DSCM)
%             ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off',...
%                 'Units','normalized', 'clipping' , 'off');
%             text(0.5, 1,[blockname,', Trigger ',num2str(trigChn),', Stim ',num2str(stimChn),', ylim ',num2str(ylNorm),...
%                 ', Window ',num2str(-back*1000),' +',num2str(window*1000),'ms'],...
%                 'HorizontalAlignment' ,'center','VerticalAlignment', 'top')
%         end
        
    end
end

set(0, 'CurrentFigure', allChns);
if(DSCM)
    subplot(10,10,1); axis off;
    text(0,1.5,[blockname,', Trigger ',num2str(trigChns),', Stim ',num2str(stimChn),', ylim ',num2str(ylNorm),...
        ', Window ',num2str(-back*1000),' +',num2str(window*1000),'ms'],...
        'fontsize',7)
%     ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off',...
%         'Units','normalized', 'clipping' , 'off');
%     text(0.5, 1,[blockname,', Trigger ',num2str(trigChn),', Stim ',num2str(stimChn),', ylim ',num2str(ylNorm),...
%         ', Window ',num2str(-back*1000),' +',num2str(window*1000),'ms'],...
%         'HorizontalAlignment' ,'center','VerticalAlignment', 'top')
end

end
