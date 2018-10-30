function N = discretizeOverlap(data,edges)
    N = cell(1,size(edges,2));
    for i = 1:(size(edges,2))
        N{i} = find(data > edges(1,i) & data < edges(2,i));
    end
end