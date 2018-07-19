function ind = findfirstineachcolumn(tf)
%
% for logical array tf, finds row number for each column corresponding to
% first 1

N = size(tf,2);

ind = nan(1, N);

% Modified to check from the back and find when it becomes zero to make sure random movement in the beginning does not change results.
% RJY
for j = 1:N
    prev = 0;
    for i = size(tf,1):-1:1
        if (~tf(i,j) && prev)
            ind(j) = i+1;
            break
        end
        prev = tf(i,j);
    end
end


% for j = 1:N
%     for i = 1:size(tf,1)
%         if tf(i,j)
%             ind(j) = i;
%             break
%         end
%     end
% end

