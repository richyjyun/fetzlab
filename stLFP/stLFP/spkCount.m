% Print number of spikes per channel/sortcode combination
% returns cell/sortcode pairs that are above count
% RJY 08/16/2018

function [cells,sortcodes] = spkCount(Snip,count)

if(~exist('count','var'))
    count = 5000;
end

cells = []; sortcodes = [];
chns = unique(Snip.chan);
for c = 1:length(chns)
    codes = Snip.sortcode(Snip.chan == chns(c));
    codes = unique(codes);
    codes(codes==0 | codes==31) = [];
    for sc = 1:length(codes)
        spks = sum(Snip.chan==chns(c) & Snip.sortcode == codes(sc));
        fprintf('%d %d: %d\n',chns(c),codes(sc),spks);
        if(spks >= count)
            cells(end+1) = chns(c);
            sortcodes(end+1) = codes(sc);
        end
    end
end

end