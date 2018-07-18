function [Pre,Post] = getPhaseHist(tankpath,blockname,trigChns,Codes,times)

for t = 1:size(times,1)
    
    T1 = times(t,1); %- window;
    T2 = times(t,2); %+ window;% get all LFPs
    LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs;
    % get all snippets from spike sorting
    Snips = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',3,'STORE','eNe1'); Snips = Snips.snips.eNe1;
    
    % define variables
    fs = LFPs.fs;
    
    % filter
    Beta = bpfilt(LFPs.data',[30,50],fs,3)';
        
    % get timestamps
    trig = (Snips.ts(Snips.chan == trigChns & Snips.sortcode == Codes)-T1)*fs;
    trig = round(trig);
    
    if(isempty(trig))
        return;
    end
    
    bins = -pi:pi/9:pi;
    % plot
    Counts = zeros(96,length(bins)-1);
    for j = 1:size(Beta,1)
        h = hilbert(Beta(j,:));
        ang = angle(h);
        [Counts(j,:),~] = histcounts(ang(trig),bins);
    end
    if(t == 1)
        Pre = Counts;
    else
        Post = Counts;
    end
end

end




