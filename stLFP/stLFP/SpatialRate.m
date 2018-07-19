function [spkrate, avg] = SpatialRate(Mani,spk,fs,window,dt,plt)

% 1 is RU (radial (up) is positive) and 2 is FE (flexion (left) is positive)
Mani.data(2,:) = -Mani.data(2,:); % convert for left/right

FE_lim = [min(Mani.data(2,:)),max(Mani.data(2,:))];
RU_lim = [min(Mani.data(1,:)),max(Mani.data(1,:))];

dloc = 0.1;
window = round(window*fs); dt = round(dt*fs);

FE_bins = FE_lim(1):dloc:FE_lim(2);
RU_bins = RU_lim(1):dloc:RU_lim(2);

spkrate = cell(length(FE_bins),length(RU_bins));

w1 = 1; w2 = w1+window-1;
while w2 < length(Mani.data)  
    FE = mean(Mani.data(2,w1:w2)); RU = mean(Mani.data(1,w1:w2));

    if(FE < FE_lim(1) || FE > FE_lim(2)...
            || RU < RU_lim(1) || RU > RU_lim(2))
        w1 = w1+dt; w2 = w2+dt;
        continue;
    end
    
    FE_ind = find(FE >= FE_bins,1,'last');
    RU_ind = find(RU >= RU_bins,1,'last'); 

    spkcnt = sum(spk >= w1 & spk < w2);
    
    spkrate{FE_ind,RU_ind} = [spkrate{FE_ind,RU_ind},spkcnt];
    
    w1 = w1+dt; w2 = w2+dt;
end

avg = cellfun(@mean,spkrate); % get average # spikes
avg = avg/window; % get absolute firing rate in Hz
if(plt)
    imagesc(avg); set(gca,'YDir','normal'); axis off; colorbar;
end

end