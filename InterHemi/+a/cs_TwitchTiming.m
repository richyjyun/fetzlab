load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat')
load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKato.mat')
load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgor.mat')

for i = 1:length(SL)
    if(isempty(SL(i).trig1) || ~isempty(SL(i).Bad))
        continue;
    end
    if(strcmp(SL(i).StimHemi,'R'))
        Trials = SL(i).lefttrials;
        Accel = SL(i).accel_raw_l;
    elseif(strcmp(SL(i).StimHemi,'L'))
        Trials = SL(i).righttrials;
        Accel = SL(i).accel_raw_r;
    else
        continue;
    end
    window = 300; %100ms
    window = floor(window*SL(i).fs/1000);
    
    inds = -floor(window/4):1:window;
    
    trialinds = repmat(SL(i).trig1'.*SL(i).fs./1000 + str2num(SL(i).Stim_Delay)*SL(i).fs./1000, length(inds), 1) + repmat(inds(:), 1, size(SL(i).trig1,1));
    
    Snips = Accel(floor(trialinds));
    
    figure;
    plot(inds,Snips)
    title(char(SL(i).Date))
    
end