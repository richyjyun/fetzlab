function [Theta,PLV] = getPercentiles(amp,phase,prcnt)

Theta = zeros(length(prcnt)-1,size(phase,1)); 
PLV = zeros(length(prcnt)-1,size(phase,1));
for p = 1:length(prcnt)-1
    for ang = 1:size(phase,1)
        p1 = prctile(amp(ang,:),prcnt(p)); p2 = prctile(amp(ang,:),prcnt(p+1));
        inds = amp(ang,:) >= p1 & amp(ang,:) < p2;
        
        % mean of circular values. kind of undoing the "angle"
        % function when saving phase
        r = sum(exp(1j*phase(ang,inds)));
        Theta(p,ang) = angle(r);
        PLV(p,ang) = abs(r)/sum(inds);
    end
end

end