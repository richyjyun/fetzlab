function [C] = PlotERPDiff(SL)
C = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]);

for i = 1:length(SL)
    IpsiR = median(SL(i).IpsiNeural(1:16,:));
    IpsiL = median(SL(i).IpsiNeural([17:24,27:30],:));
    ContR = median(SL(i).ContNeural(1:16,:));
    ContL = median(SL(i).ContNeural([17:24,27:30],:));

    
    
end

end