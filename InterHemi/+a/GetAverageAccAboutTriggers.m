function [AvAcc, ts] = GetAverageAccAboutTriggers(SL)
%
% pass this a single element from SL and it will calculate the average accs
% for the trigger times in that SL
%
% returns Nx2 array like [leftacc, rightacc]
%
% arb 050415

window = 0.8; % seconds

col = 1;

windowsamples = -.2*SL.fs:window*SL.fs;

accelfs = 1000;
if(isfield(SL,'accelfs'))
    accelfs = SL.accelfs;
end

inds = repmat(SL.trig1(:)*SL.fs/accelfs, 1, length(windowsamples)) + repmat(windowsamples, length(SL.trig1), 1);

inds = round(inds);

for accdat = {SL.accel_raw_l, SL.accel_raw_r} % for right and left hand accel

    inds(inds(:,end)>length(accdat{1}), :) = []; % delete rows that try to read out past end of recording

    inds(inds(:,1)<=0, :) = []; % delete rows that try to read out before start of recording

    tmpacc = accdat{1}(inds') - repmat(mean(accdat{1}(inds'),1), size(inds,2), 1); % subtract mean

    AvAcc(:,col) = mean(abs(tmpacc), 2);

    col=col+1;

end

ts = SL.fs^-1*(windowsamples);