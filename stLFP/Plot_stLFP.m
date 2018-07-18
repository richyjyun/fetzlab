function AdapterTest(blockname,trigChns,stimChn,times,packet)

%% Load
% 2. Epochs - Thd1 (threshold) and all corresponding values
% 3. Snips - Beta (snippets, 24kHz)
% 4. Streams - Mani (manipulandum, 3kHz), LFPs (3kHz), SUAs (24kHz), Filt
% (3kHz)
% 5. Scalars (all at 1Hz) - Trig (trigger channel), SUAc (SUA channels), Stim (Stim
% params and times) 
tankpath = 'Y:\~NeuroWest\Spanky\CycleTriggered-170710-143939\';
T1 = 500; T2 = 1500;  % in seconds. 0 to denote start or end of entire recording
% get all LFPs
LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs; 
% get all snippets from spike sorting
Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1'); Snips = Snips.snips.eNe1;

% Define good single unit channels for the day, can be any number of
% channels.
SUAc = [69]; 
% if wanting to use sortcodes. Should have one code per channel
Codes = [1];

% define variables
fs = LFPs.fs;
window = 0.1; %100 ms window, change as needed
range = round(-window*fs:1:window*fs); 
% yl = zeros(length(trigChns),2);
for i = 1:length(trigChns)
    figure; % new figure for each trigger channel
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Getting all snippets from a single channel. Can use the PCA function
    % to apply sortcodes too and use Snips.sortcode to determine which ones
    % you want. 
    
    trig = (Snips.ts(Snips.chan == trigChns(i))' - T1)*fs; 

%         D = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',2,'VERBOSE',0);
%         Discrim = D.epocs.Dscm.onset;
%         trig = ((Discrim - T1)*fs)';
    
    % define all indices to get data from
    trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
    trialinds(:,floor(trialinds(1,:))<=0) = [];
    trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];
    
    % loop through all LFP channels
    for j = 1:size(LFPs.data,1)
        [r,c,e] = GetWadeChannelPosition(j);
        subplot(10,10,(r-1)*10+c); hold on;
        d = LFPs.data(j,:);
        d = d(floor(trialinds));
        
        %%%%%%%%%%%%%%%%%%%
        % Might want to add thresholding code here? Or just use the PCA
        % function?
        %
        
        % get average and plot
        d = mean(d,2);
        if(j == trigChns(i))
            plot(range/fs,d,'r');
        else
            plot(range/fs,d);
        end
        
    end
      
%     YLIM = [min(yl(:,1)),max(yl(:,2))];
    for j = 1:size(LFPs.data,1)
        [r,c,e] = GetWadeChannelPosition(j);
        subplot(10,10,(r-1)*10+c);
%         ylim(YLIM);
        axis off
        sub_pos = get(gca,'position'); % get subplot axis position
        set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
        title(num2str(j),'fontsize',7)
    end
    
    % just titling. can include ylim (if normalizing) here to be able to
    % tell since axis is turned off.
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,['Trigger Channel ',num2str(SUAc(i))],'HorizontalAlignment' ,'center','VerticalAlignment', 'top')
    
    % Print to packet
    
end


