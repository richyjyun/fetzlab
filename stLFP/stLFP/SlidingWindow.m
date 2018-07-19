trig = (Snips.ts(Snips.chan == 83 & Snips.sortcode == 1)' - T1)*fs;

spectra = [];

start = 1; window = 500; step = window/5;
j = 25;
while start + window <= length(trig)
    
        trig2 = trig(start:start+window); start = start+step;
        trialinds = repmat(trig2, length(range), 1) + repmat(range(:), 1, size(trig2,2));
        trialinds(:,floor(trialinds(1,:))<=0) = [];
        trialinds(:,floor(trialinds(end,:))>length(LFPs.data)) = [];
        
        
        d = LFPs.data(j,:);
        d = d(floor(trialinds));
        d = d-mean(d);
        
        Y = fft(d);
        f = fs*(0:(size(d,1)/2))/size(d,1);
        Y = Y(1:length(f),:);
        
        Y = mean(abs(Y),2);
        spectra(end+1,:) = Y(f>10 & f<50);
        
end

figure;
imagesc(spectra); 
xticks(0:10:40)
xticklabels({'10','20','30','40','50'})
        