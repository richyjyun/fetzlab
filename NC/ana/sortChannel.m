function [spkTime,spk,noise,spkInd] = sortChannel(Data,p)
% Summary
%   Sorts the channel with the paramters. Assumes a single spike.
%
% Inputs
%   Data - vector of data
%   p - sorting parameters
%
% Outputs
%   spkTime - time of spikes (in samples)
%   spk - extracted spike waveforms
%   noise - threshold crossings determined to be noise
%
% Details
%   Sorts using PCA and 2 window discrimination using the parameters passed.
% See getSortingParams.m for parameter details
%
% RJY 06/22/2018

if any(size(Data)==1)
    %% Filter
    filt = bpfilt(Data,p.bw,p.fs,3);
    
    %% Get snippets
    cross = filt < p.thresh;
    cross = find(diff(cross) == 1)';
    cross(cross+p.win2del > length(filt)) = []; % check to make sure 2 window can be done
    
    inds = repmat(cross, length(p.range), 1) + repmat(p.range(:), 1, size(cross,2));
    inds(:,inds(1,:)<=0) = [];
    inds(:,inds(end,:)>length(filt)) = [];
    
    snips = filt(inds);
else
    snips = Data;
    cross = p.ts;
end

%% Transform into SVD space
Vr = snips'*p.Ur/p.Sr;

D = sqrt(sum((Vr-p.cent).^2,2));

bad = false(length(D),1);
bad(D>p.dlim) = true;

%% Do 2 window discrimination
win1 = snips(p.win1del,:) >= p.win1lim(1) & snips(p.win1del,:) <= p.win1lim(2);
win2 = snips(p.win2del,:) >= p.win2lim(1) & snips(p.win2del,:) <= p.win2lim(2);

bad(~(win1|win2)) = true;

%% Set outputs
spkTime = cross(~bad);
spkInd = find(~bad);
spk = snips(:,~bad);
noise = snips(:,bad);

end