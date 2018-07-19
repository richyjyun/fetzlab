function SL = SortSL(SL,field)
    if(isfield(SL,'trig1'))
        temp = extractfield(SL,field);
        [~,ind] = sort(temp);
        SL = SL(ind);
    else
        return;
    end
end