function [accs,ts] = GetAverageAccRT(SL,window)

if ~isempty(SL.trig1)
    
    start = SL.trig1(1);
    finish = SL.trig1(end);
    preinds_l = SL.lefttrials(:,1)<start;
    condinds_l = SL.lefttrials(:,1)>=start & SL.lefttrials(:,1)<=finish;
    postinds_l = SL.lefttrials(:,1)>finish;
    preinds_r = SL.righttrials(:,1)<start;
    condinds_r = SL.righttrials(:,1)>=start & SL.righttrials(:,1)<=finish;
    postinds_r = SL.righttrials(:,1)>finish;
    
    pre = calculateAcc(SL,preinds_l,preinds_r,window);
    cond = calculateAcc(SL,condinds_l,condinds_r,window);
    post = calculateAcc(SL,postinds_l,postinds_r,window);
    
    accs(:,1) = pre(:,1);
    accs(:,2) = cond(:,1);
    accs(:,3) = post(:,1);
    accs(:,4) = pre(:,4);
    accs(:,5) = cond(:,4);
    accs(:,6) = post(:,4);
    
    
else
    
    all_l = 1:length(SL.lefttrials);
    all_r = 1:length(SL.righttrials);
    all = calculateAcc(SL,all_l,all_r,window);
    accs(:,1) = all(:,1);
    accs(:,4) = all(:,4);
    accs(:,6) = 0;
    
end

ts = linspace(window(1), window(2), size(accs,1));
end

% Slightly modified version of Andrew's code
function AvAcc = calculateAcc(SL,left,right,window)
col = 1;
RFilter = u.FilterAccForMovementTimes(SL.accel_raw_r, SL.fs, 'richardson');
LFilter = u.FilterAccForMovementTimes(SL.accel_raw_l, SL.fs, 'richardson');
for tstart = {SL.lefttrials(left,1)+SL.rts_l(left)*SL.fs/1000, SL.righttrials(right,1)+SL.rts_r(right)*SL.fs/1000}; % for left and right targets
    temp = tstart{1};
    tstart{1} = temp(~isnan(temp));
    
    windowsamples = round(window(1)*SL.fs:window(2)*SL.fs);
    inds = repmat((tstart{1}*SL.fs/1000), 1, length(windowsamples)) + repmat(windowsamples, length(tstart{1}), 1);
    inds = round(inds); 
    
    for accdat = {LFilter, RFilter} % for right and left hand accel
        inds(inds(:,end)>length(accdat{1}), :) = []; % delete rows that try to read out past end of recording
        inds(inds(:,1)<=0, :) = []; % delete rows that try to read out before start of recording
        tmpacc = accdat{1}(inds');% - repmat(mean(accdat{1}(inds'),1), size(inds,2), 1); % subtract mean
        AvAcc(:,col) = median(abs(tmpacc),2);
        col=col+1;
    end
end
end