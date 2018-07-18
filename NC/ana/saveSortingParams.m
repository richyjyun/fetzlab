function saveSortingParams(fpath,fname,Channels,fs,train_time)
% Summary
%   Goes through each channel and extracts spikes and parameters using
%   getSortingParams. Plots spike, noise, and autocorrelation and asks for
%   user input on whether to keep the spike or not. Saves parameters in
%   fpath as SortingParams.mat.
%
% Inputs
%   fpath       path of the experiment
%   fname       path of the Matlab settings file
%   Channels    recorded channels with possible spikes
%   fs          sampling rate
%   train_time  time used for training the parameters in seconds
%
% RJY 06/22/2018

    % load in the data
    Data = nc3data(Channels, 0, train_time, fs, [], fname);

    % get stim times and blank
    Events = nc3events(fname);
    Stim = floor(Events.stim{1}*fs)'; blank = floor(-0.015*fs):1:floor(0.015*fs);
    rm = repmat(Stim, length(blank), 1) + repmat(blank', 1, size(Stim,2));
    rm(:,rm(1,:)<=0) = [];
    rm(:,rm(end,:)>length(Data)) = length(Data);
    
    % determine if the detected spike is good or not
    Params = cell(0); SpkChannels = [];
    for chn = 1:length(Channels)
        [p,spk,ts,noise] = getSortingParams(Data(:,chn),fs,rm);

        % for plotting
        x = (1:size(spk,1))/fs*1000;

        % plot the detected spike
        figure('pos',[50 0 800 800]); subplot(2,2,1);
        plot(x,spk(:,unique(round(linspace(1,size(spk,2),1000)))),'Color',[0.5,0.5,0.5]);
        hold on; plot(x,mean(spk,2),'k','linewidth',2);
        xlim([x(1),x(end)]); xlabel('Time (ms)'); ylabel('\muV');
        title(['Channel ',num2str(Channels(chn)),' Detected Spike']);

        % plot the noise removed
        subplot(2,2,2);
        plot(x,noise(:,unique(round(linspace(1,size(noise,2),1000)))),'Color',[0.5,0.5,0.5]);
        hold on; plot(x,mean(noise,2),'k','linewidth',2);
        xlim([x(1),x(end)]); xlabel('Time (ms)'); ylabel('\muV');
        title('Removed Noise');

        % plot the auto correlation
        subplot(2,2,3);
        CrossCorr(ts/fs, 'ts2',ts/fs,'binsize', 0.002,'lag',[-0.1,0.1],'suppress_plot',0);
        yticks([]); xlabel('Time (ms)'); title(['Autocorrelation, ',num2str(size(spk,2)),' spikes']);

        % prompt
        subplot(2,2,4); axis off
        prompt = {'Type Y in console to keep', 'Type anything else to discard'};
        text(0.5,0.5,prompt,'horizontalalignment','center');

        % get user input
        in = input(['Channel ',num2str(Channels(chn)),', Keep the spike? [Y/N]'],'s');
        close(gcf);
        if strcmpi(in,'y')
            Params{length(Params)+1} = p;
            SpkChannels(end+1) = Channels(chn);
        end

    end

    % save the parameters
    save(fullfile(fpath,'SortingParams'),'Params','SpkChannels','train_time');

end