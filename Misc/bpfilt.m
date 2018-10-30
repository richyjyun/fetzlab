function A = bpfilt (data,freq_range,sr,order)

% Band-pass elliptical filtering of a vector
%
% Syntax:
% A = bpfilt (data, freq_range, sr, order);
%
% data          vector to be filtered
% freq_range    band range in Hz e.g. [30 50]
% sr            sampling frequency in Hz
% order         filter order (by default: 3)
%
% Written by Stavros Zanos, Spring 2006

data = double(data);
Wn = freq_range./(sr/2);
% [B,C] = butter(order,Wn);
[B,C] = ellip(order,.5,20,Wn);
A = filtfilt (B,C,data);