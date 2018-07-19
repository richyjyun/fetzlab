load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbiFinal.mat');
load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbiNeuro.mat');


for S = 1
    switch S
        case 1
            SL = UbiSL;
        case 2
            SL = IgorSL;
        case 3
            SL = KatoSL;
    end
    for i = 1:length(SL)
        if(length(SL(i).Condition)<6 || ~strcmp(SL(i).Condition(1:6),'Contra') || isempty(SL(i).trig1) ||~isempty(SL(i).Bad) ...
                || SL(i).Condition(end)=='R')
            continue;
        end
        
        %determine contra and ipsi trials
        Contra = []; Ipsi = []; ContraRT = []; IpsiRT = [];
        if(strcmp(SL(i).StimHemi,'L'))
            Contra = SL(i).righttrials;
            Ipsi = SL(i).lefttrials;
            ContraRT = SL(i).rts_r;
            IpsiRT = SL(i).rts_l;
        elseif(strcmp(SL(i).StimHemi,'R'))
            Contra = SL(i).lefttrials;
            Ipsi = SL(i).righttrials;
            ContraRT = SL(i).rts_l;
            IpsiRT = SL(i).rts_r;
        end
        
        %determine which trials have stimulation
        if(length(SL(i).Condition) >=6 && strcmp(SL(i).Condition(1:6),'Contra'))
            StimTrials = Contra;
            Trials = [Contra(:,1);Ipsi(:,1)];
            RT = [ContraRT;IpsiRT];
            Label = char([ones(length(ContraRT),1)*67;ones(length(IpsiRT),1)*73]);  %uint16 code for C and I
%             ContraExp(prevInd+i) = 1;
        else%if(strcmp(SL(i).Condition(1:4),'Ipsi'))
            continue;
            StimTrials = Ipsi;
            Trials = [Ipsi(:,1);Contra(:,1)];
            RT = [IpsiRT;ContraRT];
            Label = char([ones(length(IpsiRT),1)*73;ones(length(ContraRT),1)*67]);  %uint16 code for C and I
            IpsiExp(prevInd+i) = 1;
        end
        
        StimInd = zeros(length(StimTrials),1);
        StimLag = StimInd;
        trigInd = [];
        for j = 1:length(SL(i).trig1)
            norm = StimTrials(:,1)-50-SL(i).trig1(j);
            ind = find(norm>0,1)-1;
            if(~isempty(ind))
                StimInd(ind) = 1;
                trigInd(end+1) = j;
                % Replace with subtracting RP timing 
                StimLag(ind) = StimTrials(ind,1)+RT(ind)-200-SL(i).trig1(j);
            end
        end
        StimInd = [StimInd;zeros(length(Trials)-length(StimInd),1)];
        StimLag = [StimLag;zeros(length(Trials)-length(StimLag),1)];
        
        %set up order of all trials
        [Trials,Order] = sort(Trials);
        StimInd = StimInd(Order);
        StimLag = StimLag(Order);
        Label = Label(Order);
        RT = RT(Order);
        StimHist = -1*ones(length(StimInd),1);
        for j = 1:length(StimInd)
            ind = find(StimInd(1:j),1,'last');
            if(~isempty(ind))
                StimHist(j) = j-ind;
            end
        end
         
        ContraInd = Label=='C' & StimHist==0;
        Contra = RT(ContraInd);
        ContraLag = StimLag(ContraInd);

        IpsiInd = Label=='I' & StimHist==1;
        Ipsi = RT(IpsiInd);
        IpsiInd = find(IpsiInd)-1;
        IpsiLag = StimLag(IpsiInd);
        
        ContraPre = StimHist(strfind(Label','C'));
        ContraPre = ContraRT(ContraPre==-1);
        IpsiPre = StimHist(strfind(Label','I'));
        IpsiPre = IpsiRT(IpsiPre==-1);
        
        
        
    end
end



