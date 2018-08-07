function [avg,plv] = meanPhase(data,dim)
    add = sum(exp(1j*data),dim);
    avg = angle(add);
    plv = abs(add);
end