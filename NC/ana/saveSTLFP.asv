function saveSTLFP(fpath,fname,day,session_time,Channels,fs,bp)
% Summary
%   Calculates stLFP over the course of the entire experiment for each
%   spike and each recording channel. Saves packets to fpath\Packets.
%   Requires saveSleep.m and saveSortedSpikes to have run first.
%
% Inputs
%   fpath - path of the experiment
%   fname - path of the Matlab settings file
%   day - experiment day
%   session_time - length of the experiment
%   Channels - recorded channels
%   bp - bandpass for stLFP. Empty array for no filter
%
% RJY 06/22/2018

    % load channels with spikes
    if(~exist(fullfile(fpath,'SpkChannels.mat')))
        error('Run saveSortedSpikes.m first!');
    end
    load(fullfile(fpath,'SpkChannels'));

    % set time to bin stLFPs over
    window = 5*60; %5 minute window
    dt = 1*60; %1 minute shifts

    % set with of stLFP time to look at
    frame = 0.05;
    range = round(-frame*fs):round(frame*fs);

    % set sleep times for plotting
    if(~exist(fullfile(fpath,'Sleep.mat')))
        error('Run saveSleep.m first!');
    end
    load(fullfile(fpath,'Sleep'));
    marker = nan(1,length(sleep)); % which indices to use as sleep
    marker(sleep) = find(range==0); % plot at time 0 of stLFP, nan otherwise to leave blank
    sleepy = (1:length(marker))/sleepfs/dt; % put it into minutes
    shift = round(window/2*sleepfs); % plot for middle of each bin
    marker = marker(shift:end); sleepy = sleepy(shift:end);

    % loop through each channel with a spike
    packetPath = fullfile(fpath,'Packets');
    if(~exist(packetPath,'dir'))
        mkdir(packetPath)
    end
    for chn = 1:length(SpkChannels)

        % load the spike
        load(fullfile(fpath,['spk_',num2str(SpkChannels(chn))]))

        % set packet
        packet = [packetPath,'\',day,'_chn',num2str(SpkChannels(chn)),'_stLFP.ps'];
        if(exist(packet))
            delete(packet);
        end

        % loop through each recorded channel
        for ch = 1:length(Channels)

            start = 0;
            stLFP = [];

            % go through entire recording
            while start+window < session_time

                fprintf('TrigChn %d, Chn %d, Time %ds, %d%% done\n',SpkChannels(chn),Channels(ch),round(start),round(start/session_time*100));

                % set start and end times
                st = start-frame;
                if(st<0)
                    st = start;
                end
                wd = window+2*frame;
                if(wd>session_time)
                    wd = window;
                end

                % load data and bp filter
                Data = nc3data(Channels(ch), st, wd, fs, [], fname);
                if(~isempty(bp))
                    Data = bpfilt(Data,bp,fs,3);
                end

                % set triggers
                trig = spkTime-st;
                trig = trig(trig>=0 & trig<=wd);
                trig = round(trig*fs);

                inds = repmat(trig, length(range), 1) + repmat(range(:), 1, size(trig,2));
                inds(:,inds(1,:)<=0) = [];
                inds(:,inds(end,:)>length(Data)) = [];

                % get stLFP
                d = Data(inds); d = d - mean(d);
                stLFP = [stLFP,mean(d,2)];

                % increment loop
                start = start+dt;

            end

            % save data
            lfpPath = fullfile(fpath,'stLFP');
            if(~exist(lfpPath,'dir'))
                mkdir(lfpPath)
            end
            save(fullfile(lfpPath,['stLFP_',num2str(SpkChannels(chn)),'_',num2str(ch)]),...
                'stLFP','range','fs','window','dt','bp','chn','ch','-v7.3');
            
            % plot figure
            figure('visible','off'); imagesc(stLFP');
            xticks(1:((length(range)-1)/4):length(range))
            xticklabels(-frame*1000:(frame*2000/4):frame*1000);
            xlabel('Snippet Time (ms)');
            ylabel('Experiment Time (min)');
            box off; c = colorbar;
            ylabel(c,'Amplitude (uV)');
            title(['NC Channel ',num2str(Channels(ch))]);

            % plot sleep times
            hold on; 
            plot(marker,sleepy,'Color',[0,0,0],'linewidth',2);

            % print
            print('-painters','-fillpage',gcf, '-dpsc2', packet, '-append');
            close(gcf);

        end

        % convert to pdf
        callps2pdf(packet);

    end
end
