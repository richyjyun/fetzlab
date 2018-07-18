clear; close all;

%% Load
% 2. Epochs - Thd1 (threshold) and all corresponding values, Dscm -
% Discrimination times
% 3. Snips - Beta (snippets, 24kHz)
% 4. Streams - Mani (manipulandum, 3kHz), LFPs (3kHz), SUAs (24kHz), Filt
% (3kHz)
% 5. Scalars (all at 1Hz) - Trig (trigger channel), SUAc (SUA channels), Stim (Stim
% params and times)
tankpath = 'Y:\~NeuroWest\Spanky\CycleTriggered-170710-143939\';
blockname = 'Spanky-170725-145504';


%% Find when stimulation started
% For some reason Stim in scalars did not save. We can find when
% stimulation began by looking at snips saved by the spike sorter, then
% count the number of discrims from there on out. START is when the first
% stimulation occured. 
T1 = 1620; T2 = 1860;  % in seconds. 0 to denote start or end of entire recording
% plot(TT.snips.eNe1.data(TT.snips.eNe1.chan == 20,:)');
TT = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2);
trigSnips = TT.snips.eNe1.data(TT.snips.eNe1.chan == 20,:)';
stim = zeros(size(trigSnips,2),1);
for i = 1:length(stim)
    if(max(trigSnips(:,i))>1e-4 && min(trigSnips(:,i))<-2e-4)
        stim(i) = 1;
    end
end
stimStart = find(stim,1);
trigChns = find(TT.snips.eNe1.chan == 20); snipNum = trigChns(stimStart);
START = TT.snips.eNe1.ts(snipNum); 

%% Load in all of epoch data to have the discrims, and figure out time ranges
Epochs = TDT2mat([tankpath,blockname],'TYPE',2);
Discrim = Epochs.epocs.Dscm;
[~,ind] = min(abs(Discrim.onset-START));
ind = ind-1; % last discrimination index before stimulation began

testPhase = [ind-4999,ind;ind+10001,ind+15000;ind+25001,ind+30000;ind+40001,ind+45000];
% testPhase = [ind-14999,ind;ind+40001,ind+55000];
testTime = Discrim.onset(testPhase);

%% Load in appropriate data and plot
window = 0.1; 
yl = zeros(size(testTime,1),96,2);
for i = 1:size(testTime,1)
    figure;
    T1 = testTime(i,1) - 0.01;
    T2 = testTime(i,2) + window;
    LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs;
    fs = LFPs.fs; range = round(-0.01*fs:1:window*fs); %50 ms window for capturing CCEPs
    trig = (Discrim.onset(testPhase(i,1):testPhase(i,2))' - T1)*fs;
    trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
    trialinds(:,floor(trialinds(1,:))<=0) = [];
    trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];
    for j = 1:size(LFPs.data,1)
        [r,c,e] = GetWadeChannelPosition(j);
        subplot(10,10,(r-1)*10+c);
        d = LFPs.data(j,:);
        d = mean(d(floor(trialinds)),2);
%         yl = [-15,5];
        if(j == 20)
            plot(range/fs,d,'r');
        else
            plot(range/fs,d);
        end
        %         ylim(yl*1e-6);
        yl(i,j,:) = ylim;
        axis off
        sub_pos = get(gca,'position'); % get subplot axis position
        set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
        title(num2str(j),'fontsize',7)
    end
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,['Test Phase ',num2str(i)],'HorizontalAlignment' ,'center','VerticalAlignment', 'top')
end


%%
LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs;
Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1'); Snips = Snips.snips.eNe1;
Epochs = TDT2mat([tankpath,blockname],'TYPE',2);
Discrim = Epochs.epocs.Dscm;%Stim = Epochs.epocs.
Snips = TDT2mat([tankpath,blockname],'TYPE',3,'NODATA',true);
Scalars = TDT2mat([tankpath,blockname],'TYPE',5);




figure;
trig = (Discrim.onset' - T1)*fs; %trig = trig(1:3000);
trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
trialinds(:,floor(trialinds(1,:))<=0) = [];
trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];
for j = 1:size(LFPs.data,1)
    [r,c,e] = GetWadeChannelPosition(j);
    subplot(10,10,(r-1)*10+c);
    d = LFPs.data(j,:);
    d = mean(d(floor(trialinds)),2);
    yl = [-15,5];
    if(j == 20)
        plot(range/fs,d,'r');
    else
        plot(range/fs,d);
    end
    %         ylim(yl*1e-6);
    axis off
    sub_pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
    title(num2str(j),'fontsize',7)
end
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,['Channel ',num2str(SUAc(i)),', ylim: ',num2str(yl),'uV'],'HorizontalAlignment' ,'center','VerticalAlignment', 'top')