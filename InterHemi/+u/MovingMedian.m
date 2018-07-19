% data is given in two rows, with the first being times and second being
% reaction times. window is defined in seconds. returns vals, with first
% row the time stamps and the second with median values
function vals = MovingMedian(data,window)

limit = data(1,max(find(~isnan(data(1,:)))));
medians = nan(2,ceil(limit/window));

for i = 0:window:limit
    snip = [];
    if(i+window > limit)
        bound = limit;
    else
        bound = i+window;
    end
    snip = data(2,(data(1,:)>i & data(1,:)<=bound));
    m = nanmedian(snip);
    if(isempty(snip) || isnan(m))
        continue
    else
        medians(2,i/window + 1) = nanmedian(snip); 
    end
    medians(1,i/window + 1) = i;
end

trials = medians(1,:);
times = medians(2,:);
vals(1,:) = trials;%(~isnan(trials));
vals(2,:) = times;%(~isnan(times));

end