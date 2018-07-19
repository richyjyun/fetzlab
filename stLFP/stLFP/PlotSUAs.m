function fig = PlotSUAs(tankpath,blockname,times)
    Snips = TDT2mat([tankpath,blockname],'T1',times(1),'T2',times(2),'TYPE',3,'STORE','eNe1'); 
    Snips = Snips.snips.eNe1;
    
    fig = figure;
    for chn = 1:96
        [c,r,e] = GetWadeChannelPosition(chn);
        subplot(10,10,(r-1)*10+c);
        title(num2str(chn));
        axis off;
        ChnInd = Snips.chan == chn;
        Codes = unique(Snips.sortcode(ChnInd));
        if(isempty(Codes))
            continue;
        end
        for sc = 1:length(Codes)
            inds = ChnInd & Snips.sortcode == Codes(sc);
            if(Codes(sc) == 31 || sum(inds) < 1000)
                continue;
            end
            
            spikes = Snips.data(inds,:);
            
            plot(mean(spikes),'linewidth',2); hold on;
        end
        title(num2str(chn));
        axis off;
    end
    
end