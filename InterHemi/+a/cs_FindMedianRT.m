mediansum = 0;
trialsum = 0;
varsum = 0;

for i = 1:3
    
    switch i
        case 1
            load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgor.mat')
            constring = 'nostim';
        case 2
            load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKato.mat')
            constring = 'NaN';
        case 3
            load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat')
            constring = 'Control';
    end
    
    for j = 1:length(SL)
        if(strcmp(SL(j).Condition, constring))
            varsum = varsum + nanstd(SL(j).rts_l)^2/sum(~isnan(SL(j).rts_l));
            varsum = varsum + nanstd(SL(j).rts_r)^2/sum(~isnan(SL(j).rts_r));
            mediansum = mediansum + nanmedian(SL(j).rts_l)*sum(~isnan(SL(j).rts_l));
            mediansum = mediansum + nanmedian(SL(j).rts_r)*sum(~isnan(SL(j).rts_r));
            trialsum = trialsum + sum(~isnan(SL(j).rts_l));
            trialsum = trialsum + sum(~isnan(SL(j).rts_r));
        end
    end
    
end

medianRT = mediansum/trialsum;
semRT = sqrt(varsum);