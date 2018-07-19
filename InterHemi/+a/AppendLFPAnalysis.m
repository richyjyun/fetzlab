function SL = AppendLFPAnalysis(SL, technique)
%
% plots a full page of details for a single experiment
%
% arb feb 20 2015

import OOARB.Utils.*

SL(1).LFP_chnm = [];

if ~exist('technique', 'var'), technique = 'envelope'; end 

for i = 1:numel(SL) % for all elements

    disp(['appending lfp ' num2str(100*i/numel(SL)) '% done'])
    
    if (SL(i).Bad), continue; end
    fname = [char(SL(i).Date),'_',char(SL(i).Session_Guger_Train(1))];
    
    if ~exist([fname,'.cfg']), continue; end
    dwn = 1;
    [dat, fs, chnm, trigch] = u.LoadGug(fname,dwn); %first four are left hem. last four are right hem

    if isempty(SL(1).LFP_chnm), SL(1).LFP_chnm = chnm; end
    
%     if fs~=SL(i).fs, error('somethin wrong'); end

    for filter = {'beta', 'gamma'}
    
        for ii = 1:size(dat, 2)
            
            %if strcmp(filter{1}, 'beta') % only do this once
                %SL(i,1).spectrum(:,ii) = ana.GetSpectrum(predat(:,ii), fs); % this goes from 0 to nyquist
                %SL(i,2).spectrum(:,ii) = ana.GetSpectrum(condat(:,ii), fs);
                %SL(i,3).spectrum(:,ii) = ana.GetSpectrum(posdat(:,ii), fs);
            %end

            eval(['[sig, avp] = InstBand(dat(:,ii), fs,''' technique ''',''' filter{1} ''');'])
            eval(['SL(i).av' filter{1} '(ii) = avp;'])
            eval(['[SL(i).TA_L' filter{1} '(:,ii), ts] = GetTriggeredAv(sig, SL(i).lefttrials(:,1), 1000, fs);']);
            eval(['SL(i).TA_R' filter{1} '(:,ii) = GetTriggeredAv(sig, SL(i).righttrials(:,1), 1000, fs);']);
            
            SL(i).TA_ts = [ts(1) ts(end)]; % append the range on the lags
            
        end
        
    end

end

function [TAv, ts] = GetTriggeredAv(sig, ind, dur, fs)
%
% for signal passed through sig, gathers windows starting at ind for dur
% samples and returns average
%

prewindowtime = .200; % sec to look before window

ind = 1000*(ind/1000); % convert to samples
dur = fs*(dur/1000);

windowinds = -prewindowtime*fs:dur;

inds = repmat(ind(:), 1, length(windowinds)) + repmat(windowinds(:)', length(ind), 1);

inds = round(inds);

inds(any(isnan(inds),2),:) = []; % delete nan

inds(inds(:,1)<1,:) = []; % delete rows that try to look before 0

inds(inds(:,end)>length(sig), :) = []; % delete rows that try to read out past end of recording

TAv = median(sig(inds'),2);

ts = windowinds/fs;



function [sig, avp] = InstBand(sig, fs, technique, band)

if ~exist('technique', 'var'), technique = 'envelope'; end

sig = sig-mean(sig); % subtract DC

if strcmp(band, 'beta')
    bandpass = [15 25];
elseif strcmp(band, 'gamma')
    bandpass = [45 120];
end

%avpnorm = mean(sig.^2);

% if strcmp(technique, 'envelope')
%     norm = mean(abs(hilbert(sig)));
% elseif strcmp(technique, 'rms')
%     norm = mean(abs(sig-mean(sig)));
% end

Wn_theta = [bandpass(1)/(fs/2) bandpass(2)/(fs/2)]; % normalized by the nyquist frequency
[btheta,atheta] = butter(3,Wn_theta);
sig = double(sig);
sig = filtfilt(btheta,atheta,sig);

avp = mean(sig.^2);%/avpnorm;

if strcmp(technique, 'envelope')
    sig = abs(hilbert(sig));
elseif strcmp(technique, 'rms')
    rms = 5; % low pass cutoff hz
    [brms,arms] = butter(1,rms/(fs/2),'low'); % lowpass (see daq_sapi_*)
    sig = max(filtfilt(brms,arms,abs(sig)),zeros(size(sig))); % filter with lowpass
else
    error('problem')
end

%sig = sig/norm; % fraction of average total power