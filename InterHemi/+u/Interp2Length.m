function D = Interp2Length(d,len)
    int = floor(len/(length(d)-2));
    D = zeros(1,len);
    if int == 0
        D = d;
        return;
    end
    for i =  1:length(d)-1
        ind = (i-1)*int+1;
        D(ind:ind+int) = linspace(d(i),d(i+1),int+1);
    end
    if(length(D)>len)
        D = D(1:len);
        return;
    else
        ext = d(end-1)-d(end);
        D(ind+int:end) = linspace(d(end),d(end)-ext,len-(ind+int-1));
    end
end