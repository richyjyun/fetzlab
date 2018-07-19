function [Stim,idx,stimT] = getStimTrials(trig,Trials,err)
Stim = nan(1,length(trig)); idx = nan(1,length(trig)); stimT = [];
for j = 1:length(trig)
%     ind = find(Trials(:,1) > (trig(j)+50),1)-1;
%     
%     if isempty(ind)
%         continue;
%     end
%     
%     if(abs(Trials(ind,2)-trig(j)) > abs(Trials(ind+1,1)-trig(j)))
%         ind = ind+1;
%     end

    ind = find(Trials(:,1) < (trig(j)+50) & Trials(:,2) > trig(j)-err);
    if(isempty(ind))
        continue;
    end

    Stim(j) = ind(1);
    idx(j) = j;
    stimT(end+1) = trig(j) - Trials(ind(1),1);
end
Stim(isnan(Stim)) = []; 
Stim = unique(Stim);
idx(isnan(idx)) = [];
end
