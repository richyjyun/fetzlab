tankpath = 'Y:\~Neuro1\transfer\Hugo-180904';
blockname = 'rec6';

% for finding pre/post times
Stim = TDT2mat(fullfile(tankpath,blockname),'Type',4,'Store','Syn2');
Stim = Stim.streams.Syn2;
times = [0,find(Stim.data > 0.5,1)-1;...
    find(Stim.data >0.5,1,'last')+1,length(Stim.data)];
times = times./Stim.fs;

%% Pre & Post
trigChn = 90; stimChn = [74,76];
for t = 1:2
    T1 = times(t,1); T2 = times(t,2);
    LFPs = TDT2mat(fullfile(tankpath,blockname),'T1',T1,'T2',T2,'Type',4,'Store','RawD');
    LFPs = LFPs.streams.RawD; fs= LFPs.fs;
    
    tLFP = bpfilt(LFPs.data(trigChn,:),[15,25],fs,3);
    thresh = std(tLFP);
    h = hilbert(tLFP); phase = angle(h); amp = abs(h);
    
    trig = find([0,diff(phase)]<-pi & amp > thresh);
    
    window = [-0.05,0.05];
    range = round(window(1)*fs:1:window(2)*fs);
    
    trialinds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
    trialinds(:,floor(trialinds(1,:))<=0) = [];
    trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];
    
    %% Loop through all LFP channels and get stLFPs
    stLFPs = zeros(size(LFPs.data,1),length(range));
    for j = 1:size(LFPs.data,1)
        d = LFPs.data(j,:);
        d = d(floor(trialinds));
        d = d - mean(d);
        d = mean(d,2);
        stLFPs(j,:) = d;%zscore(d);
    end
    
    for j = 1:size(LFPs.data,1)
        [c,r,~] = GetWadeChannelPosition(j);
        subplot(10,10,(r-1)*10+c);
        
        if(t == 1)
            plot(range/fs,stLFPs(j,:),'k');
        else
            plot(range/fs,stLFPs(j,:),'r');
        end
        
        if(j == trigChn)
            title(num2str(j),'fontsize',7,'Color','b');
        elseif(any(j==stimChn))
            title(num2str(j),'fontsize',7,'Color','r');
        else
            title(num2str(j),'fontsize',7);
        end
        xlim([window(1),window(2)]); axis off;
        hold on;
        yl = ylim;
        plot([0,0],yl,'k:');
    end
end

