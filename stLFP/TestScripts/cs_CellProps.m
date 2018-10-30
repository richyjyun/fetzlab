close all; clear all; pack

%% Loading in appropriate variables
vals = GetGoogleSpreadsheet('1WLfx_3Zq1MdA2T0S6-LUTU0QqMe68vDLM8caA5_lEMc');

% choose day
date = '20180914';

% Arrange into a struct. First row is the fields
SL = struct;
for i = 2:size(vals,1)
    if(strcmp(vals(i,2),date))
        for j = 1:size(vals,2)
            SL.(char(vals(1,j))) = char(vals{i,j});
        end
    end
end

%% Get channels and corresponding sort codes
Channels = split(SL.Channels,'/');
Codes = split(SL.SortCodes,'/');

chns = []; codes = [];
for i = 1:length(Channels)
    c = split(Codes(i),',');
    for j = 1:length(c)
        chns(end+1) = str2double(Channels(i));
        codes(end+1) = str2double(c(j));
    end
end

%% Find correct block
tankpath = 'Y:\~NeuroWest\Spanky\RandomStim-180314-124242\';
blocks = dir(tankpath);
blocks = blocks([blocks.isdir]);
blocks = extractfield(blocks,'name')';
date = date(3:end);
blockname = char(blocks(find(~cellfun(@isempty,strfind(blocks,date)))));

%% Loop through blocks, make new packet for each, and plot everything 
for b = 1:size(blockname,1)
    packet = ['F:\S\Packets\CellProps\',blockname(b,end-12:end),'.ps'];
    
    if(exist(packet,'file'))
        delete(packet);
    end
    
%     % find times
    blockpath = [tankpath,blockname(b,:)];
%     TT = TDT2mat(blockpath,'TYPE',2);
%     Dscm = TT.epocs.Dscm;
%     [val,ind] = findpeaks(Dscm.data);
%     ind = ind(val>1000); val = val(val>1000);
%     times = [ind(1)-val(1),ind(1)];  %just get first epoch
%     times = Dscm.onset(times);
    
    times = [1,20*60];

    % load data
    window = 0.1; dt = 0.02; %time window and steps to look at 
    T1 = times(1) - window; T2 = times(2) + window;
    Mani = TDT2mat(blockpath,'T1',T1,'T2',T2,'TYPE',4,'STORE','Mani'); Mani = Mani.streams.Mani;
    Snips = TDT2mat(blockpath,'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1'); Snips = Snips.snips.eNe1;
    fs = Mani.fs;

    % Plot spike waveform, autocorrelation, and directional tuning
    spkLim = 2000; bad = [];
    for c = 1:length(chns)
        ind = Snips.chan == chns(c) & Snips.sortcode == codes(c);
        
        if(sum(ind) < spkLim)
            bad(end+1) = c;
            continue;
        end
        
        ttl = sprintf('%s, %d s,Chn%d, Code%d, %d Spks',blockname(b,:),round(T2-T1),chns(c),codes(c),sum(ind));
        disp(ttl); %#ok<DSPS>
        
        fig = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(fig,'visible','off');
        
        % spike waveform
        subplot(3,2,1);
        snips = Snips.data(ind,:); sample = floor(linspace(1,size(snips,1), 100));
        plot(snips(sample,:)','Color',[0.7,0.7,0.7]); hold on; plot(mean(snips),'k','LineWidth',2);
        ylim([min(min(snips(sample,:))),max(max(snips(sample,:)))]); 
        xlim([1,size(snips,2)]); axis off;
        title(ttl);

        % autocorrelation
        subplot(3,2,2);
        CrossCorr(Snips.ts(ind), 'ts2',Snips.ts(ind),'binsize', 0.002,'lag',[-window,window],'suppress_plot',0); axis off;
        title('Autocorrelation');
        
        % directional tuning
        subplot(3,2,3:6)
        spk = (Snips.ts(ind)' - T1)*fs;
        DirectionalTuning(Mani,spk,fs,window,dt,1);
        title('Directional Tuning');
        
        print('-painters',fig, '-dpsc2', packet, '-append');
        close(fig);
        
    end
    
    % plot cross correlations
    fig = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(fig,'visible','off');
    disp([blockname,', Cross Correlations'])
    Chns = chns; Chns(bad) = []; 
    Codes = codes; Codes(bad) = [];
    sp = length(Chns)+1;
    for i = 1:length(Chns)
        ind1 = Snips.chan == Chns(i) & Snips.sortcode == Codes(i);
        
        snips = Snips.data(ind1,:); sample = floor(linspace(1,size(snips,1), 100));
        subaxis(sp, sp, 1, i+1, 'spacing', 0, 'padding', 0.001)
        plot(snips(sample,:)','Color',[0.7,0.7,0.7]); hold on; plot(mean(snips),'k','LineWidth',2);
        axis off; title([num2str(Chns(i)),',',num2str(Codes(i)),',',num2str(size(snips,1))],'fontsize',5)
        subaxis(sp, sp, i+1, 1, 'spacing', 0, 'padding', 0.001)
        plot(snips(sample,:)','Color',[0.7,0.7,0.7]); hold on; plot(mean(snips),'k','LineWidth',2);
        axis off; title({[num2str(Chns(i)),',',num2str(Codes(i))],num2str(size(snips,1))},'fontsize',5)
        
        for j = i:length(Chns)
            ind2 = Snips.chan == Chns(j) & Snips.sortcode == Codes(j);
            
            subaxis(sp, sp, i+1, j+1, 'spacing', 0, 'padding', 0.001)
            window = 0.2;
            bin = 0.002;
            [cor,lags] = CrossCorr(Snips.ts(ind1), 'ts2',Snips.ts(ind2),'binsize', bin,'lag',[-window,window],'suppress_plot',0); 
            axis off;
            if i~=j 
               ylim([min(cor),max(cor)]) 
            end
        end
    end
    
    subaxis(sp,sp,1,1,'spacing', 0, 'padding', 0.001)
    str = sprintf('Win %dms\nBin %dms',round(window*1000),round(bin*1000));
    text(0,1,str,'HorizontalAlignment','left','VerticalAlignment','top','fontsize',7);
    axis off;
    
    print('-painters',fig, '-dpsc2', packet, '-append');
    close(fig);
    
end

callps2pdf(packet);

disp('DONE')



