%% Load
% 2. Epochs - Thd1 (threshold) and all corresponding values
% 3. Snips - Beta (snippets, 24kHz)
% 4. Streams - Mani (manipulandum, 3kHz), LFPs (3kHz), SUAs (24kHz), Filt
% (3kHz)
% 5. Scalars (all at 1Hz) - Trig (trigger channel), SUAc (SUA channels), Stim (Stim
% params and times) 
tankpath = 'Y:\~NeuroWest\Spanky\CycleTriggered-170710-143939\';
blockname = 'Spanky-170726-151731';
T1 = 300; T2 = 900;  % in seconds. 0 to denote start or end of entire recording
LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs; 
Triggers = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',2,'STORE','LFPs'); 
Trig(1) = Triggers.epochs.Cl1_; Trig(2) = Triggers.epochs.Cl2_;

fs = LFPs.fs;
window = 0.1; range = round(-0.01*fs:1:window*fs); %50 ms window for capturing CCEPs
for i = 1:length(Trig)
    figure;
    trig = (Trig(i).onset - T1)*fs; % might be onset
    trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
    trialinds(:,floor(trialinds(1,:))<=0) = []; 
    trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];
    
    for j = 1:size(LFPs.data,1)
        [r,c,e] = GetWadeChannelPosition(j);
        subplot(10,10,(r-1)*10+c);
        d = LFPs.data(j,:);
        d = mean(d(floor(trialinds)),2);
        
        plot(range/fs,d);
        
        axis off
        sub_pos = get(gca,'position'); % get subplot axis position
        set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
        title(num2str(j),'fontsize',7)
    end
   
end
