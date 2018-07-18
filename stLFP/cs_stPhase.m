clear
clc
close all

%% load in epocs
tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';
blockname = 'Spanky-180309-132433';
% blockname = 'Spanky-180226-150458';
T1 = 1; T2 = 0;
Epocs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',2); Onset = Epocs.epocs.Dscm.onset; Epocs = Epocs.epocs.Dscm.data;
fs = TDT2mat([tankpath,blockname],'T1',1,'T2',1.01,'TYPE',4); fs = fs.streams.LFPs.fs;
SUAc = 32;
stim = [64];

%% Search for intervals
dEpocs = cat(1,diff(Epocs),0);
Interends = find(dEpocs < -500);
for i = 1:length(Interends)
    for j = (Interends(i)-1):-1:1
        if dEpocs(j) < 1
            Interstart(i,1) = j + 1;
            break
        end
    end
end

%% Load data
window = 0.05;
range = round(-window*fs:1:window*fs);
nostimtrial = [1,3]; % trials with no stim
% fbands = [5:3:101];
fbands = [15,25];
y = cell(length(nostimtrial),96);
yy = y;
% yl = nan(length(nostimtrial)*96,2);
yl = [];
sc = zeros(length(nostimtrial),1);
phase = cell(length(fbands)-1,96,length(nostimtrial));

figure

for k = 1:length(nostimtrial)
    onset = Onset(Interstart(nostimtrial(k)):Interends(nostimtrial(k)))+0.195;
    sc(k) = (Interends(nostimtrial(k))-Interstart(nostimtrial(k)))/(max(onset)-min(onset));
    LFPs = TDT2mat([tankpath,blockname],'T1',min(onset),'T2',max(onset),'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs.data; 
    Raw = LFPs;
%     LFPs = transpose(bpfilt(LFPs',[5,150],fs,3)); 
    trig = (onset'-min(onset))*fs; 
    trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
    trig(floor(trialinds(end,:))>length(LFPs)) = [];
    trig(floor(trialinds(1,:))<=0) = []; 
    trialinds(:,floor(trialinds(1,:))<=0) = []; 
    trialinds(:,floor(trialinds(end,:))>length(LFPs)) = [];
    %% plot mean
%     for j = 1:size(LFPs,1)
%         [c,r,e] = GetWadeChannelPosition(j); 
% %         d = LFPs(j,:);
%         d = Raw(j,:);
%         d = d(floor(trialinds));
%         y{k,j} = d;
%         subplot(10,10,(r-1)*10+c);
%         plot(mean(d,2)-mean(mean(d,2))); hold on
%         if j == SUAc
%             title(num2str(j),'fontsize',7,'color','r')
%         elseif ismember(j,stim)
%             title(num2str(j),'fontsize',7,'color','b')
%         else
%             title(num2str(j),'fontsize',7)
%         end
%         if(j ~= SUAc)
%             yl((k-1)*96+j,:) = ylim;
%         end
%     end
    %% find bands  
    for i = 1:length(fbands)-1
        disp(i)
        bands = transpose(bpfilt(Raw',[fbands(i),fbands(i+1)],fs,3));
        ang = angle(hilbert(bands')');
        for j = 1:size(LFPs,1)
            dd = bands(j,:);
            dd = dd(floor(trialinds));
            yy{k,j} = dd;
            angj = ang(j,:);
            phase{i,j,k} = angj(floor(trialinds));
            [c,r,e] = GetWadeChannelPosition(j); 
            subplot(10,10,(r-1)*10+c);
            
            
            r = sum(exp(1j*phase{i,j,k}),2);
            
            plot(angle(r))
%             plot(mean(dd,2) - mean(mean(dd,2)))
            hold on
            if j == SUAc
                title(num2str(j),'fontsize',7,'color','r')
            elseif ismember(j,stim)
                title(num2str(j),'fontsize',7,'color','b')
            else
                title(num2str(j),'fontsize',7)
            end
            yl = [yl;ylim];
        end
    end
end

YLIM = [nanmedian(yl(:,1)),nanmedian(yl(:,2))]*2;
for j = 1:size(LFPs,1)
    [c,r,e] = GetWadeChannelPosition(j);
    subplot(10,10,(r-1)*10+c);
    hold on
    line([154,154],[-10,10],'Color','black','LineStyle','--')
    ylim(YLIM);
    axis off
    sub_pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
    title(num2str(j),'fontsize',7)
end

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,['Trigger Channel ',num2str(SUAc),', ',num2str(length(trig)),' Triggers, Ylim ',num2str(YLIM)],'HorizontalAlignment' ,'center','VerticalAlignment', 'top')







% A = phase{1,34};
% B = A(154,:);
% [z,I] = sort(B);
% d = Raw(64,:);
% d = d(floor(trialinds));
% newA = d(:,I);
% newA = newA - mean(newA,1);
% plot(mean(newA(:,1:5000),2))
% hold on
% plot(mean(newA(:,1:7000),2))
% hold on
% plot(mean(newA(:,1:9000),2))