function saveSortedSpikes(fpath,fname,session_time,fs)
% Summary
%   Sorts channels using sortChannel.m and saves the spike channel,
%   parameters, and times to fpath. Must run saveSortingParams first.
%
% Inputs
%   fpath           path of the experiment
%   fname           path of the Matlab settings file
%   session_time    length of the experiment
%   fs              sampling rate
%
% RJY 06/22/2018

    % load parameters saved by saveSortingParams.m
    if(~exist(fullfile(fpath,'SortingParams.mat')))
        error('Run saveSortingParams.m first!');
    end
    load(fullfile(fpath,'SortingParams'));

    step = 3600; % look at an hour at a time
    times = step:step:session_time;
    times = [times,session_time];

    % for each detected channel
    for chn = 1:length(SpkChannels)

        spkTime = [];

        % loop through time bins
        Spike = [];
        for t = times

            % load data
            disp(['Channel ',num2str(SpkChannels(chn)),', Loading ',num2str(t-step),'s to ',num2str(t),'s']);
            Data = nc3data(SpkChannels(chn),t-step, step, fs, [], fname);

            % sort
            disp('Sorting...');
            [s,spk,~] = sortChannel(Data,Params{chn});
            Spike = [Spike,spk];
            spkTime = [spkTime,s/fs+t-step];

            disp([num2str(round(find(times==t)/length(times)*100*chn/length(SpkChannels))),'% Done']);
        end

        % save extracted spike times
        channel = SpkChannels(chn);
        p = Params{chn};
        save(fullfile(fpath,['spk_',num2str(channel)]),'channel','p','Spike','spkTime','-v7.3');
        
    end
    
    save(fullfile(fpath,'SpkChannels'),'SpkChannels'); % save channels with spikes for later use

end