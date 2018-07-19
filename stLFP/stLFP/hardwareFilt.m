function [Filter] = hardwareFilt(data,bp,fs)

% Filter data exactly as the TDT system filters it
% Biquad filter with 2nd order Butterworth. Need to do low and high pass separately instead of a bandpass. 
[z,p,k] = butter(2,bp(2)/(fs/2),'low');  % low pass filter
sos = zp2sos(z,p,k);
Filter = sosfilt(sos,data);
[z,p,k] = butter(2,bp(1)/(fs/2),'high'); % high pass filter
sos = zp2sos(z,p,k);
Filter = double(sosfilt(sos,Filter));

end