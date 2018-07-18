function [p,spk,ts,noise] = getSortingParams(Data,fs,rm)
% Summary
%   Determines sorting parameters using train data for further 
% sorting. Assumes only one cell in each channel.
%
% Inputs 
%   Data    vector of data
%   fs      sampling rate
%
% Outputs
%   p       parameters
%   spk     extracted spike waveforms
%   ts      time stamp of spikes
%   noise   threshold crossings determined to be noise
%
% Details
%   Obtains snippets of data using an arbitrary threshold. Performs 
% PCA and sets a distance from the centroid of the simple cluster.
% Determines distribution of trough and peak values and timings for two
% window discrimination. Use sortChannel.m with the parameters. 
%
%
% Returned parameter values:
% For SVD
%   p.Ur = Ur;
%   p.Sr = Sr;
%   p.cent = centroid;
%   p.dlim = dlim;
%
% For 2 window discrim. delays are in samples. 
%   p.win1del = round(nanmean(ttime));
%   p.win1lim = troughlims;
%   p.win2del = round(nanmean(ptime));
%   p.win2lim = peaklims;
%
% General
%   p.range = range;
%   p.fs = fs;
%   p.thresh = thresh;
%   p.bw = bw;
%
% RJY 06/22/2018

    %% Filter data
    % Bandwidth for bandwidth filter 
    bw = [300,3000];
    filt = bpfilt(Data,bw,fs,3);
    filt(rm) = 0;
    % filt = hardwareFilt(Data,bw,fs); % can do hardware-like filter

    %% Get snippets
    window = 0.001; % window to look around threshold crossing (s)
    window = round(window*fs);
    range = round(-window/2):round(1.3*window);

    thresh = -std(filt); % threshold for detecting spiks
    cross = filt < thresh; % threshold crossings
    spks = find(diff(cross) == 1)'; % times of crossings

    inds = repmat(spks, length(range), 1) + repmat(range(:), 1, size(spks,2));
    inds(:,inds(1,:)<=0) = [];
    inds(:,inds(end,:)>length(filt)) = [];

    snips = filt(inds); % Snippets. All analysis performed on this data
    rm = max(abs(snips))> 3*abs(thresh); % set amplitude upper limit
    snips(:,rm) = [];

    %% PCA and Clustering
    % SVD space
    [U,S,V] = svd(snips,'econ');

    lambda = (S.^2)/(length(S)-1); % variance 
    r = find(cumsum(diag(lambda)/sum(diag(lambda))) > 0.95,1); % get the rank that accounts for more than 95% of the variance

    % Low dimensional space
    Ur = U(:,1:r);
    Sr = S(1:r,1:r);
    Vr = V(:,1:r);

    % Centroid of the cluster in the low dimensional space
    centroid = mean(Vr);
    D = sqrt(sum((Vr-centroid).^2,2)); % distance from cluster to each point

    % Determine distance threshold
    q1 = quantile(D,0.25);
    q3 = quantile(D,0.75);
    dlim = q3+1.5*(q3-q1);
    bad = false(length(D),1);
    bad(D>dlim) = true; % set which ones are noise

    %% 2 Window Discrimination Parameters
    % Find values/time of peak/trough and time between the two
    tdiff = nan(size(snips,2),1);
    trough = nan(size(snips,2),1);
    peak = nan(size(snips,2),1);
    ptime = nan(size(snips,2),1);
    ttime = nan(size(snips,2),1);
    for i = 1:size(snips,2)
        [t,troughs] = findpeaks(-snips(round(window/2):end,i));
        [p,peaks] = findpeaks(snips(round(window/2):end,i));
        if(isempty(troughs) || isempty(peaks))
            continue;
        else
            ptime(i) = peaks(1);
            ttime(i) = troughs(1);
            tdiff(i) = peaks(1)-troughs(1);
            trough(i) = t(1);
            peak(i) = p(1);
        end
    end
    tdiff(tdiff<0) = nan;
    trough = -trough;
    trough(trough>0) = nan;
    peak(peak<0) = nan;

    % Determine threshold for each parameter
    q1 = quantile(tdiff,0.25);
    q3 = quantile(tdiff,0.75);
    tdifflims = [q1-1.5*(q3-q1),q3+1.5*(q3-q1)];

    q1 = quantile(trough,0.25);
    q3 = quantile(trough,0.75);
    troughlims = [q1-1.5*(q3-q1),q3+1.5*(q3-q1)];

    q1 = quantile(peak,0.25);
    q3 = quantile(peak,0.75);
    peaklims = [q1-1.5*(q3-q1),q3+1.5*(q3-q1)];

    % Define which ones lie outside the thresholds
    bad = (bad | isnan(tdiff) | tdiff < tdifflims(1) | tdiff > tdifflims(2)...
        |trough < troughlims(1) | trough > troughlims(2) | peak < peaklims(1) | peak > peaklims(2));

    %% Set function outputs
    spk = snips(:,~bad);
    ts = spks(~bad);
    noise = snips(:,bad);

    p = struct;

    % For SVD
    p.Ur = Ur;
    p.Sr = Sr;
    p.cent = centroid;
    p.dlim = dlim;

    % For 2 window discrim. delays are in samples. 
    p.win1del = round(nanmean(ttime));
    p.win1lim = troughlims;
    p.win2del = round(nanmean(ptime));
    p.win2lim = peaklims;

    % General
    p.range = range;
    p.fs = fs;
    p.thresh = thresh;
    p.bw = bw;

end