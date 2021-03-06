clear; close all;

fpath = 'Y:\~NeuroWest\Spanky\Neurochip\S_20180613_07\';
chn = 28; fs = 20000;

% [data, names, session_time] = nc3data(chn, 0, 100, 20000, [300,3000], path);

index = strfind(fpath, '\');
fprefix = fpath(index(end-1)+1:end-1);
fname = [fpath,fprefix,'_Chan' num2str(chn,'%.02d') '.i16'];

%% initialize sorting parameters
init_time = 600; % time to use 
[Data, end_sec] = nc3chan(0, init_time, fs, [], fs, fname);

[p,spk,noise] = getSortingParams(Data,fs);

%% sort and find spike times
step = 3600; % look at an hour at a time
times = step:step:end_sec;
times = [times,end_sec];

spkTime = [];
for t = times
    tic
    disp(['Loading ',num2str(t-step),'s to ',num2str(t),'s']);
    [Data, ~] = nc3chan(t-step, step, fs, [], fs, fname);

    disp('Sorting...');
    [s,~,~] = sortChannel(Data,p);
    spkTime = [spkTime,s/fs+t-step];
    
    toc
    
    disp([num2str(round(find(times==t)/length(times)*100)),'% Done']);

end

save([fpath,'spk_',num2str(chn)],'chn','init_time','p','spkTime','-v7.3');

%% stLFP over time

load([fpath,'spk_',num2str(chn)])

tic
window = 5*60; %5 minute window
dt = 1*60; %1 minute shifts

packet = ['F:\S\Packets\NC\',fprefix,'_chn',num2str(chn),'_2.ps'];
if(exist(packet))
    delete(packet);
end

frame = 0.05;
range = round(-frame*fs):round(frame*fs);

[sleep,sleepfs] = getSleepTimes(fpath,'Y:\~NeuroWest\Spanky\Neurochip\S_20180605_01\S_20180605_01.mat');
marker = nan(1,length(sleep));
marker(sleep) = find(range==0);
sleepy = (1:length(marker))/sleepfs/dt; % put it into minutes
shift = round(window/2*sleepfs);
marker = marker(shift:end); sleepy = sleepy(shift:end);

end_sec = length(sleep)/sleepfs;

for chn = 17:32
    
    start = 0;
    stLFP = [];
    
    % chn = 25;
    fname = [fpath,fprefix,'_Chan' num2str(chn,'%.02d') '.i16'];

    
    while start+window < end_sec
        
        disp(['Chn',num2str(chn),', Time ',num2str(start),'s, ',num2str(round(start/end_sec*100)),'%']);
        
        st = start-frame;
        if(st<0)
            st = start;
        end
        wd = window+2*frame;
        if(wd>end_sec)
            wd = window;
        end
        
        [Data, ~] = nc3chan(st, wd, fs, [], fs, fname);
        Data = bpfilt(Data,[15,50],fs,3); %3rd order bp filter
        
        trig = spkTime-st;
        trig = trig(trig>=0 & trig<=wd);
        trig = round(trig*fs);
        
        inds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
        inds(:,inds(1,:)<=0) = [];
        inds(:,inds(end,:)>length(Data)) = [];
        
        d = Data(inds); d = d - mean(d);
        stLFP = [stLFP,mean(d,2)];
        
        start = start+dt;
    end
    
    % map1 = hot; map1 = map1(24:end,:);
    % map2 = fliplr(map1);
    % map = [map1(1:end-1,:);flipud(map2)];
    
    figure('visible','off'); imagesc(stLFP'); %colormap(map);
    xticks(1:((length(range)-1)/4):length(range))
    xticklabels(-frame*1000:(frame*2000/4):frame*1000);
    xlabel('Snippet Time (ms)');
    ylabel('Experiment Time (min)');
    box off; c = colorbar;
    ylabel(c,'Amplitude (uV)');
    title(['NC Channel ',num2str(chn)]);
    
    hold on; yl = ylim;
    plot(marker,sleepy,'Color',[0,0,0],'linewidth',2);
    
    print('-painters','-fillpage',gcf, '-dpsc2', packet, '-append');
    close(gcf);
    
end

callps2pdf(packet);

toc



