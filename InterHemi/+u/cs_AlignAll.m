clear; close all;

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat')

for i = 1:length(SL)
    
    D = SL(i).Date;
    S = SL(i).Session_Guger_Train;
    Session = [char(D),'_',char(S(2))];
    disp(['Session ',Session])
    
    if(~isempty(SL(i).Bad) || exist([Session,'.f32']))
        continue;
    end
    
    if(~isempty(SL(i).trig2))
        u.trainalign3(Session,2);
    elseif(~isempty(SL(i).trig1))
        u.trainalign3(Session,1);
    end
    
end