function [f,Y] = Spectra(data,fs)

Y = fft(data,[],1);
f = fs*(0:(size(data,1)/2))/size(data,1);
Y = Y(1:length(f),:);
Y = (abs(Y).^2)/size(data,1);

end