clear; %close all;

%% Load
% 2. Epochs - Thd1 (threshold) and all corresponding values
% 3. Snips - Beta (snippets, 24kHz)
% 4. Streams - Mani (manipulandum, 3kHz), LFPs (3kHz), SUAs (24kHz), Filt
% (3kHz)
% 5. Scalars (all at 1Hz) - Trig (trigger channel), SUAc (SUA channels), Stim (Stim
% params and times) 
tankpath = 'Y:\~NeuroWest\Spanky\CycleTriggered-170710-143939\';
blockname = 'Spanky-171207-145138';
% T1 = 60; T2 = 0;  % in seconds. 0 to denote start or end of entire recording
% get all LFPs
% LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs; 
% fs = LFPs.fs;

% Define good single unit channels for the day, can be any number of
% channels.
SUAc = [69]; 
% if wanting to use sortcodes. Should have one code per channel
Codes = [1];
DSCM = 1;

if(DSCM)
    D = TDT2mat([tankpath,blockname],'TYPE',2,'VERBOSE',0);
    ind1 = find(D.epocs.Dscm.data == 10)-20; ind2 = find(D.epocs.Dscm.data == 3700);
    T1 = D.epocs.Dscm.onset(ind1);
    T2 = D.epocs.Dscm.onset(ind2);
    Discrim = D.epocs.Dscm.onset([ind1:ind2]);
    LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs;
    fs = LFPs.fs;
    trig = ((Discrim - T1)*fs)'; 
else
    % get all snippets from spike sorting
    Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1'); Snips = Snips.snips.eNe1;
    trig = (Snips.ts(Snips.chan == SUAc(i) & Snips.sortcode == Codes(i))' - T1)*fs;
end



% define variables
window = 0.03; %100 ms window, change as needed
range = round(-window*fs:1:window*fs); 
yl = zeros(length(SUAc),2);
for i = 1:length(SUAc)
%     figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); %set(gcf,'visible','off');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Getting all snippets from a single channel. Can use the PCA function
    % to apply sortcodes too and use Snips.sortcode to determine which ones
    % you want. 
    
%     trig = (Snips.ts(Snips.chan == SUAc(i) & Snips.sortcode == Codes(i))' - T1)*fs; 

%         D = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',2,'VERBOSE',0);
%         Discrim = D.epocs.Dscm.onset;
%         trig = ((Discrim - T1)*fs)';
    
    % define all indices to get data from
    trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
    trialinds(:,floor(trialinds(1,:))<=0) = [];
    trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];
    
    % loop through all LFP channels
    for j = 1:size(LFPs.data,1)
        [c,r,e] = GetWadeChannelPosition(j);
        subplot(10,10,(r-1)*10+c); hold on;
        
%         bandpass = [15,300];
%         Wn_theta = [bandpass(1)/(fs/2) bandpass(2)/(fs/2)]; % normalized by the nyquist frequency
%         [btheta,atheta] = butter(3,Wn_theta);
%         filtered = filtfilt(btheta,atheta,double(LFPs.data(j,:)));
%         d = filtered;
        
        d = bpfilt([15,300],fs,LFPs.data(j,:),3);
        d = d(floor(trialinds));
        
        %%%%%%%%%%%%%%%%%%%
        % Might want to add thresholding code here? Or just use the PCA
        % function?
        %
        
        % get average and plot
        d = mean(d,2);
        if(j == SUAc(i))
            plot(range/fs,d,'r');
        else
            plot(range/fs,d);
        end
        if(j ~= SUAc(i) && j ~=25)
            yl(j,:) = ylim;
        end
        
    end
      
    YLIM = [nanmedian(yl(:,1)),nanmedian(yl(:,2))].*3/2;
    for j = 1:size(LFPs.data,1)
        [c,r,e] = GetWadeChannelPosition(j);
        subplot(10,10,(r-1)*10+c);
%         ylim(YLIM);
        xlim([-window,window]);
%         axis off
%         sub_pos = get(gca,'position'); % get subplot axis position
%         set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
        title(num2str(j),'fontsize',7)
    end
    
    % just titling. can include ylim (if normalizing) here to be able to
    % tell since axis is turned off.
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,['Trigger Channel ',num2str(SUAc(i))],'HorizontalAlignment' ,'center','VerticalAlignment', 'top')
    
end

%% Spike stats
Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1'); Snips = Snips.snips.eNe1;

Chn = SUAc;
ind = Snips.chan == Chn & Snips.sortcode == 1;

% Plot snippets
snips = Snips.data(ind,:); subplot(10,10,1); title('Snips','fontsize',7)
plot(snips','Color',[0.7,0.7,0.7]); hold on; plot(mean(snips),'k','LineWidth',2); axis off;

% Plot autocorrelation
[r,lags] = xcorr(mean(snips),mean(snips));
subplot(10,10,10); plot(lags,r,'k'); title('AutoCorr','fontsize',7); axis off;

% Plot ISI
intervals = diff(Snips.ts(ind)); intervals = intervals(intervals < 0.5); 
subplot(10,10,91); H = histogram(intervals); title('ISI','fontsize',7); axis off;

% % Raster plot
% figure;
% for i = 1:96
%     ind = Snips.chan == i & Snips.sortcode == 1;
%     times = Snips.ts(ind); hold on;
%     scatter(times,ones(length(times),1)*i,3,'k','filled');
% end


