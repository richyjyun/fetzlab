function snips = meanSubtract(snips,params)
    snips = snips-mean(snips);
    snips = snips./repmat(std(snips),size(snips,1),1); % divide by standard deviation to normalize
%     snips = rmlinesc(snips,params,[],[],60);
end