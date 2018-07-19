% From Mark Kramer (Boston University)
% http://math.bu.edu/people/mak/sfn/tutorial.pdf

v0 = NSnips(:,1) .* hann(length(NSnips(:,1)))';             %Multiply data by Hann taper.
pow = (abs(fft(v0)).^2) / length(v0);     %Compute the spectrum.
pow = 10*log10(pow/max(pow));             %Convert to decibels.
pow = pow(1:length(NSnips(:,1))/2+1);              %Ignore negative frequencies.

f0 = 1920;                                %Determine the sampling frequency.
fNQ = f0 / 2;                             %Determine the Nyquist frequency.
df = fNQ/length(pow);                           %Determine the frequency resolution.
faxis = (0:df:fNQ);                       %Construct frequency axis.

figure
plot(faxis(1:length(pow)), pow);
xlim([0 100])
xlabel('Frequency [Hz]')
ylabel('Power [dB]')