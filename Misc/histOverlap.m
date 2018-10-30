function N = histOverlap(data,edges)
N = zeros(1,size(edges,2));
for i = 1:(size(edges,2))
    N(i) = sum(data > edges(1,i) & data < edges(2,i));
end

end