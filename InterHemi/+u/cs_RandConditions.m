Total = [];
Conditions = {'Ipsi_0','Ipsi_2','Ipsi_4','Ipsi_6','Ipsi_R','Contra_6','Contra_R'};
rep = [3,1,2,3,3,1,1];
for i = 1:length(Conditions)
    Total = [Total; repmat(Conditions(i), rep(i),1)];
end

Order = Total(randperm(sum(rep)));