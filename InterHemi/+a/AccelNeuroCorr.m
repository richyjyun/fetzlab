fname = '20170208_02';  % use aligned data, since matching accelerometer with neural

[accdata, trig1, trig2, lefttrials, righttrials, fs] = u.LoadTrain([fname,'.f32']);

accdown = [];
dwn = 10;
accdown(:,1) = decimate(accdata(:,1),dwn);
accdown(:,2) = decimate(accdata(:,2),dwn);
fsT = fs/dwn;
trig1 = floor(trig1./dwn);
lefttrials = floor(lefttrials./dwn);
righttrials = floor(righttrials./dwn);

RFilter = utils.FilterAccForMovementTimes(accdown(:,2), fsT, 'richardson');
LFilter = utils.FilterAccForMovementTimes(accdown(:,1), fsT, 'richardson');

% Load Neural data
[data, fsG, chnm, itrigch] = u.LoadGug(fname, dwn);

% Filter Neural data (band pass, try to get TRPs)
bpf = [10 30]; % chL/R transform: 1. bandpass filter (Hz) (around beta wave range)
[bbpf,abpf] = butter(1,bpf/(fs/2)); % 1st order bandpass (see daq_sapi_*)
test = filtfilt(bbpf,abpf,double(data));

correlations = [];
lags = [];
figure;
window = 0.500; %s
ind = fsG*window;
t = lefttrials;
a = test(:,12);
b = LFilter';
for i = 1:length(t)
    if t(i,2)-ind <= 0
        continue;
    end
    figure(1)
    hold on
    plot(a(t(i,2)-ind:t(i,2)))
    figure(2)
    hold on
    plot(b(t(i,2)-ind:t(i,2)))
    x = a(t(i,2)-ind:t(i,2));
    y = b(t(i,2)-ind:t(i,2));
    [c,l] = xcorr(x,y,'coeff');
    %     keyboard
    lags(end+1) = l(find((c) == max((c))));
    correlations(end+1) = corr(x,y);
end

figure
plot(lags)

% So far -
% 1. Right hemisphere - lefttrials, left hand movement - high correlation
% (sanity check)
% 2. Right hemisphere - righttrials, right hand movement - high correlation
% (unexpected for this to be so high)
% 3. Right hemisphere - righttrials, left hand movement - no change in
% correlation through conditioning phase. We were hoping to see changes
% (drop in correlation).
% 4. Right movement - righttrials, left hand movement - high correlation.
% Looks like he moves his opposite hand a bit for all trials. Same the
% other way around
% 5. Lefttrials - RFilter and LFilter - much higher correlation in
% conditioning period as expected.
% 6. lefttrials - left and right hemispheres (R3.5,L4.5) - relatively high correlation,
% but very consistent throughout. OPTIMAL LAG FOR MAXIMUM CORRELATION
% CHANGES DRASTICALLY DURING CONDITIONING PERIOD. (from ~100 to ~150) The
% neural response in the left hemisphere (stimulated hemisphere) changes in
% timing and amplitude during conditioning for left trials (response to
% stim?)
% 7. righttrials - left and right hemispheres (R3.5,L4.5)- extremely high correlation
% with close to zero offset throughout. Much better fit than 6. Left
% hemisphere affects right hemisphere more than vice versa?

% Maybe can't look at crosscorrelation, or need to limit the offset. Need a
% better way to check similiarities / existence of activity.
% TRY STORING THE LAGS AND SEE IF THEYRE CONSISTENT


% Split up into regions (pre/cond/post) using trig1
% 265, 806

