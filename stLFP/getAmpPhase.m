function [freq,amp,phase] = getAmpPhase(data,trialinds,fs,f1,f2,df)

% Loops through frequency range in 2hz periods and get amplitude
% and phase
freq = f1:df:f2; amp = []; phase = []; 
for f = 1:length(freq)-1
    filt = bpfilt(data',[freq(f),freq(f+1)],fs,3)';
    d = filt;
    d = d(floor(trialinds));
    h = hilbert(d);
    temp = abs(h);
    amp(f,:) = temp(ceil(size(temp,1)/2),:);
    temp = abs(hilbert(filt));
    temp = angle(h);
    phase(f,:) = temp(ceil(size(temp,1)/2),:);
end

end