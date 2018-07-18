function PlotPLV(tankpath,blockname,trigChn,stimChn,times,packet)

fig = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); %set(gcf,'visible','off');
plv = [];

for t = 1:size(times,1)
    if(t == 1)
        color = 'k';
    elseif(t == 2)
        color = 'r';
    elseif(t == 3)
        color = 'b';
    elseif(t == 4)
        color = 'g';
    else
        color = 'y';
    end
    
    window = 0.1;
    
    T1 = times(t,1) - window;
    T2 = times(t,2) + window;% get all LFPs
    LFPs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',4,'STORE','LFPs'); LFPs = LFPs.streams.LFPs;
    Dscm = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',2,'VERBOSE',0); Dscm = Dscm.epocs.Dscm;
    
    LFPs.data = bpfilt(LFPs.data',[15,25],LFPs.fs,3)';
    
    fs = LFPs.fs;
    range = round(-window*fs:1:window*fs);
    
    trig = (Dscm.onset'-T1)*fs;
    trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
    trialinds(:,floor(trialinds(1,:))<=0) = [];
    trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];
    
    %     if(isempty(trialinds))
    %         continue;
    %     end
    
    d = LFPs.data(trigChn,:);
    d = d(floor(trialinds));
    d = d - mean(d);
    trigSite = angle(hilbert(d));
    
    yl = nan(96,2);
    % loop through all LFP channels and get spike coherence
    for j = 1:size(LFPs.data,1)
        
%         if(j == trigChn)
%             continue;
%         end
        
        d = LFPs.data(j,:);
        d = d(floor(trialinds));
        d = d - mean(d);
        testSite = angle(hilbert(d));
        
        PLV = abs(trigSite-testSite);
        PLV(PLV>pi) = 2*pi-PLV(PLV>pi);
        PLV = exp(1i*PLV);
        
        PLV = mean(PLV,2);
        PLV = abs(PLV);
        PLV = PLV-mean(PLV);
        
        [c,r,~] = GetWadeChannelPosition(j);
        subplot(10,10,(r-1)*10+c); hold on;
        
        plot(range/fs,PLV,'Color',color); xlim([-0.05,0.05]);
        
        yl(j,:) = ylim;
    end
    
end

YLIM = nanmedian(yl);
y = diff(YLIM)/2;
YLIM = [YLIM(1)-y,YLIM(2)+y];

% set y lims and other parts
for j = 1:size(LFPs.data,1)
    [c,r,~] = GetWadeChannelPosition(j);
    subplot(10,10,(r-1)*10+c);
%     ylim(YLIM);
    hold on;
    line([0 0], yl(j,:), 'linestyle', ':', 'color', [.5 .5 .5]);
    axis off
    sub_pos = get(gca,'position'); % get subplot axis position
    set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
    if(j == trigChn)
        title(num2str(j),'fontsize',7,'Color','r')
    elseif(j == stimChn)
        title(num2str(j),'fontsize',7,'Color','b')
    else
        title(num2str(j),'fontsize',7)
    end
end

end