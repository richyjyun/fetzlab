function SL = AppendNeuralDelay(SL)
for i = 1:length(SL)
    D = SL(i).Date;
    S = SL(i).Session_Guger_Train;
    Session = [char(D),'_',char(S(2))];
    disp(['Session ',Session])
    trainfile = [Session,'.f32'];
    gugfile = Session;
    if(~exist(trainfile))
        continue;
    end
    
    dwn = 10;
    
    % they should have the same fs
    [~, trig1, ~, lefttrials, righttrials, lefttrialsuccess, righttrialsuccess] = u.LoadTrain(trainfile,dwn);
    [data, fs, chnm, ~] = u.LoadGug(gugfile, dwn);
   
    bpf = [200 479]; % chL/R transform: 1. bandpass filter (Hz) (around beta wave range)
    [bbpf,abpf] = butter(1,bpf/(fs/2)); % 1st order bandpass (see daq_sapi_*)
    Filter = filtfilt(bbpf,abpf,double(data));
    
    % Assuming stimulating in left hemisphere
    Trials = righttrials;
    Hemi = find(any(char(chnm)' == 'L')');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First, determine when stimulation actually occurs
    %
    window = 0.6;
    window = window*fs;
    
    % need to figure out when stimulation occured
    stim = zeros(length(trig1),1);
    back = floor(window/10);
    fore = window;
    ratio = 2;
    for j = 1:length(trig1)
        ind = ceil(trig1(j)/dwn);
        if ind+fore < length(Filter) && ind-back > 0
            ind = ceil(trig1(j)/dwn);
            snip = Filter(ind-back:ind+fore,33);
            % assume stimulation begins at 2*rms. CHECK WITH SMALLER
            % AMPLITUDE STIMULATION TO MAKE SURE THIS IS THE CASE
            if isempty(find(snip > rms(snip)*ratio,1))
                stim(j) = NaN;
            else
                spike = find(snip > rms(snip)*ratio);
                interval = diff(spike);
                inrow = 4;
                m = movsum(interval,inrow);
                s = spike(find(m(floor(inrow/2):end-1)<50,1)+2);
                if isempty(s)
                    s = spike(1);
                end
                stim(j) = s- back + ind;
            end
        end
    end
    
    difference = stim - ceil(trig1/dwn);
    nanmedian(difference)/fs
    nanmean(difference)/fs
    
    figure
    plot(snip)
    hold on
    plot([0,length(snip)],[rms(snip)*ratio,rms(snip)*ratio]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract TRPs from each channel
    %
    window = 0.3;
    window = window*fs;
    
    bpf = [1 10]; % chL/R transform: 1. bandpass filter (Hz) (around beta wave range)
    [bbpf,abpf] = butter(1,bpf/(fs/2)); % 1st order bandpass (see daq_sapi_*)
    Filter2 = filtfilt(bbpf,abpf,double(data));
    
    close all;
    ch = 23;
    figure(1); subplot(2,1,1); title('Conditioning: Filtered')
    figure(2); subplot(2,1,1); title('Conditioning: Raw')
    for j = 1:length(trig1)
        ind = find(Trials(:,1)<trig1(j)+50,1,'last');
        if(~SL(i).righttrialsuccess(ind))
            continue;
        end
        %         rt = SL(i).rts_r(ind); rt = rt*fs/1000;
        ind = ceil(Trials(ind,1)/dwn);
        if ind+window < length(Filter) && (ind-window/2) > 0
            %         ind = find(Trials(:,1)<SL.trig1(i)+50,1,'last');
            snip = Filter2(ind-window/2:ind+window,ch);
            figure(1);subplot(2,1,1); hold on; plot((-window/2):(window),snip);
            %             plot([window/2+rt,window/2+rt],[-max(snip),max(snip)]);
            snip = Filter(ind-window/2:ind+window,ch);
            figure(2);subplot(2,1,1); hold on; plot((-window/2):(window),snip);
            %             plot([window/2+rt,window/2+rt],[-max(snip),max(snip)]);
        end
    end
    
    figure(1); subplot(2,1,2); title('Pre-conditioning: Filtered')
    figure(2); subplot(2,1,2); title('Pre-conditioning: Raw')
    for j = 1:200
        ind = ceil(Trials(j,1)/dwn);
        if ((ind+window) < length(Filter)) && ((ind-window/2) > 0)
            snip = Filter2(ind-window/2:ind+window,ch);
            figure(1);subplot(2,1,2); hold on; plot((-window/2):(window),snip);
            snip = Filter(ind-window/2:ind+window,ch);
            figure(2);subplot(2,1,2); hold on; plot((-window/2):(window),snip);
        end
    end
    
    frame = [0.1,0.3];
    frame = floor(frame*fs);
    trp = NaN(length(trig1),1);
    for j = 1:length(trig1)
        ind = find(Trials(:,1)<trig1(j)+50,1,'last');
        if(~SL(i).righttrialsuccess(ind))
            continue;
        end
        ind = ceil(Trials(ind,1)/dwn);
        if ind+frame(2) < length(Filter2) 
            snip = Filter2(ind+frame(1):ind+frame(2));
            trp(j) = find(snip == max(snip)) + ind + frame(1);
        end
    end
    
    delays = trp-stim;
    
    
    for ch = Hemi
        delays = zeros(length(trig1),1);
        for j = 1:length(SL.trig1)
            % add 50 to trig1 since at 0 delay the trigger sometimes occurs before
            % the trial
            ind = find(Trials(:,1)<SL.trig1(j)+50,1,'last');
            snip = Filter(Trials(ind,1):Trials(ind,1)+window,ch);
            delays(j) = SL.trig1(j) - (Trials(ind,1) + find(snip == max(snip),1));
        end
    end
    
end
end