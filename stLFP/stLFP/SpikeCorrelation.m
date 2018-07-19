close all; clear all; pack

%% Loading in appropriate variables
vals = GetGoogleSpreadsheet('1WLfx_3Zq1MdA2T0S6-LUTU0QqMe68vDLM8caA5_lEMc');

% Arrange into a struct. First row is the fields
SL = struct([]);
for i = 2:size(vals,1)
%     if(~strcmp(vals(i,1),'Spike Trig'))
%         continue;
%     end

    cur = length(SL)+1;
    for j = 1:size(vals,2)
        SL(cur).(char(vals(1,j))) = char(vals{i,j});
    end
end

% Loading in data
% Find correct block
date = '20180713';
tankpath = 'Y:\~NeuroWest\Spanky\IFNN\';

% tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';

ind = find(strcmp(date,extractfield(SL,'Date')));

% Get channels and corresponding sort codes
Channels = split(SL(ind).Channels,'/');
Codes = split(SL(ind).SortCodes,'/');

chns = []; codes = [];
for i = 1:length(Channels)
    c = split(Codes(i),',');
    for j = 1:length(c)
        chns(end+1) = str2double(Channels(i));
        codes(end+1) = str2double(c(j));
    end
end

blocks = dir(tankpath);
blocks = blocks([blocks.isdir]);
blocks = extractfield(blocks,'name')';

date = date(3:end);
blockname = char(blocks(find(~cellfun(@isempty,strfind(blocks,date)))));
% blockname = blockname(2,:);


%% Plotting cross correlations
% load in block 
% ADD IN CODE TO FIND TIMES

TT = TDT2mat([tankpath,blockname],'TYPE',2);
Dscm = TT.epocs.Dscm;
[val,ind] = findpeaks(Dscm.data); 
val(end+1) = Dscm.data(end); ind(end+1) = length(Dscm.data);
ind = ind(val>1000); val = val(val>1000); 
tests = 3; 
times = [ind(tests)-val(tests),ind(tests)-val(tests)+val(1)];  % times = [ind(1) - val(1),ind(1)-val(1)+15000;ind(3)-val(3),ind(3)-val(3)+10000];
times = Dscm.onset(times);

mlim = 2000;

yl = [];

% for i = 1:size(times,1)
    T1 = times(1); T2 = times(2);
    Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3);
    Snips = Snips.snips.eNe1;
    
    bad = [];
    for i = 1:length(chns)
        ind1 = Snips.chan == chns(i) & Snips.sortcode == codes(i);
        if(sum(ind1)<mlim)
            bad(end+1) = i;
        end
    end
    chns(bad) = []; codes(bad) = [];
    
    figure;
    sp = length(chns)+1;
    for i = 1:length(chns)
        
        ind1 = Snips.chan == chns(i) & Snips.sortcode == codes(i);
               
        snips = Snips.data(ind1,:); sample = floor(linspace(1,size(snips,1), 100));
        subaxis(sp, sp, 1, i+1, 'spacing', 0, 'padding', 0.001)
        plot(snips(sample,:)','Color',[0.7,0.7,0.7]); hold on; plot(mean(snips),'k','LineWidth',2);
        axis off; title([num2str(chns(i)),',',num2str(codes(i)),',',num2str(size(snips,1))],'fontsize',7)
        subaxis(sp, sp, i+1, 1, 'spacing', 0, 'padding', 0.001)
        plot(snips(sample,:)','Color',[0.7,0.7,0.7]); hold on; plot(mean(snips),'k','LineWidth',2);
        axis off; title([num2str(chns(i)),',',num2str(codes(i)),',',num2str(size(snips,1))],'fontsize',7)
        
        for j = i:length(chns)
            ind2 = Snips.chan == chns(j) & Snips.sortcode == codes(j);
            subaxis(sp, sp, i+1, j+1, 'spacing', 0, 'padding', 0.001)
            window = 0.2;
            bin = 0.005;
            [cor,lags] = CrossCorr(Snips.ts(ind1), 'ts2',Snips.ts(ind2),'binsize', bin,'lag',[-window,window],'suppress_plot',0); 
            axis off;
%             yl(end+1,:) = [min(cor),max(cor)];
            if i~=j 
               ylim([min(cor),max(cor)]) 
            end
        end
    end
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,[date,', ' ,num2str(T2-T1),'s, Test ',num2str(i)])
% end

%% Plot only autocorrelations
chns = 54; codes = 1;
for i = 1:length(chns)
    figure;

    ind1 = Snips.chan == chns(i) & Snips.sortcode == codes(i);
    snips = Snips.data(ind1,:); 
    if(size(snips,1) == 0)
        continue;
    end
    sample = floor(linspace(1,size(snips,1), 100));
    subplot(2,1,1)
    plot(snips(sample,:)','Color',[0.7,0.7,0.7]); hold on; plot(mean(snips),'k','LineWidth',2); 
    axis off; title([num2str(chns(i)),',',num2str(codes(i)),',',num2str(size(snips,1))],'fontsize',7)
    subplot(2,1,2)
    CrossCorr(Snips.ts(ind1), 'ts2',Snips.ts(ind1),'binsize', 0.002,'lag',[-0.2,0.2],'suppress_plot',0); axis off;

end


%% Plot a specific cell
for i = 1:size(times,1)
    % load in block
    T1 = times(i,1); T2 = times(i,2); tankpath = 'Y:\~NeuroWest\Spanky\CycleTriggered-170710-143939\';
    Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3);
    Snips = Snips.snips.eNe1;
    
    figure;
    subplot(2,1,1);
    ind1 = Snips.chan == 32 & Snips.sortcode == 2;
    snips = Snips.data(ind1,:); sample = floor(linspace(1,size(snips,1), 100));
    plot(snips(sample,:)','Color',[0.7,0.7,0.7]); hold on; plot(mean(snips),'k','LineWidth',2);
    title(num2str(i))
    
    subplot(2,1,2);
    CrossCorr(Snips.ts(ind1), 'ts2',Snips.ts(ind1),'binsize', 0.002,'lag',[-0.2,0.2],'suppress_plot',0); axis off;
end








