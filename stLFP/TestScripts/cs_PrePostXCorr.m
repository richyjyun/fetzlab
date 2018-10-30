tankpath = 'Y:\~NeuroWest\Spanky\RandomStim-180314-124242\';
blockname = 'Spanky-180817-140524';

times = [0,15;55,70]; times = times*60;
T1 = times(1,1); T2 = times(1,2);

chn1 = 32; chn2 = 55;
Snip32 = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1','Channel',chn1);
Snip32 = Snip32.snips.eNe1;
Snip30 = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1','Channel',chn2);
Snip30 = Snip30.snips.eNe1;

Snips.data = [Snip32.data;Snip30.data];
Snips.chan = [Snip32.chan;Snip30.chan];
Snips.sortcode = [Snip32.sortcode;Snip30.sortcode];
Snips.ts = [Snip32.ts;Snip30.ts];

% plot cross correlations
chns = [32,55]; codes = [1,1];

fig = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); 
disp([blockname,', Cross Correlations'])
Chns = chns;
Codes = codes; 
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


%% after stim
T1 = times(2,1); T2 = times(2,2);

Snip32 = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1','Channel',chn1);
Snip32 = Snip32.snips.eNe1;
Snip30 = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1','Channel',chn2);
Snip30 = Snip30.snips.eNe1;

Snips.data = [Snip32.data;Snip30.data];
Snips.chan = [Snip32.chan;Snip30.chan];
Snips.sortcode = [Snip32.sortcode;Snip30.sortcode];
Snips.ts = [Snip32.ts;Snip30.ts];

% plot cross correlations
fig = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); 
disp([blockname,', Cross Correlations'])
Chns = chns;
Codes = codes; 
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
