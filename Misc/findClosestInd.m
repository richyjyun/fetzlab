function i = findClosestInd(array,value)
    [~,i] = min(abs(array-value));
end