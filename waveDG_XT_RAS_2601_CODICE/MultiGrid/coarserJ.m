function [J] = coarserJ(J,n)

% This script eliminates half the elements n-times
% n = number of eliminations (0 returns the same matrix)

if n == 0
    return
else
for j = 1:n
  i = 12;
while(i+12 <= size(J,1))
    
J((i+1):(i+12),:) = [];
i = i+12;
end
end

end


end