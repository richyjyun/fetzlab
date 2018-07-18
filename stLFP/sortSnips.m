function [spk,ts,noise] = sortSnips(snips,ts)

snips = snips';

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
    [t,troughs] = findpeaks(double(-snips(:,i)));
    [p,peaks] = findpeaks(double(snips(:,i)));
    if(isempty(troughs) || isempty(peaks))
        continue;
    else
        pind = find(peaks==max(peaks));
        tind = find(troughs==max(troughs));
        ptime(i) = peaks(pind);
        ttime(i) = troughs(tind);
        tdiff(i) = peaks(pind)-troughs(tind);
        trough(i) = t(tind);
        peak(i) = p(pind);
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
ts = ts(~bad);
noise = snips(:,bad);

end