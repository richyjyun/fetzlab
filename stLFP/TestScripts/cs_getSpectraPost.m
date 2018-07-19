clear;

% days = {'Spanky-180122-130528';'Spanky-180129-131027';'Spanky-180204-155626';...
%     'Spanky-180208-133323';'Spanky-180208-135737';'Spanky-180208-144817';...
%     'Spanky-180213-130955';'Spanky-180216-165346';'Spanky-180228-142824'};
%     
% t = {[3,5];[3,4,5];[3,5,7];[3];[1];[1];[3,5,7];[3];[3]};

days = {'Spanky-180122-130528';'Spanky-180129-131027';'Spanky-180204-155626';...
    'Spanky-180208-133323';...
    'Spanky-180213-130955';'Spanky-180216-165346';'Spanky-180228-142824'};
    
t = {[3,5];[3,4,5];[3,5,7];[3];[3,5,7];[3];[3]};


tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';

for i = 1:length(days)
    blockname = days{i};
    trigChns = 32;
    stimChn = 32;
    epochs = t{i};
    
    TT = TDT2mat([tankpath,blockname],'TYPE',2);
    Dscm = TT.epocs.Dscm;
    [val,ind] = findpeaks(Dscm.data);
    ind = ind(val>1000); val = val(val>1000);
    
    for j = 1:length(epochs)
        tests = epochs(j); times = [ind(tests)-val(tests),ind(tests)];
        times = Dscm.onset(times);
        
        SaveAmpPhase(tankpath,blockname,trigChns,1,stimChn,times,j);
    end
end
