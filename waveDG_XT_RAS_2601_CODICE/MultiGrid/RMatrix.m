function [R] = RMatrix(n,k)

% Creates a restriction matrix for FEM matrices with 2nd order polynomials
% These are made of 6x6 block matrices, in pairs of 2 (uh,wh)
% n = matrix dimension
% k = number of halvings of the elements (0 returns an identity)

R = sparse(eye(n));

if k == 0
    return
else
    
for j = 1:k % loop over n. of halvings
  i = 12;
  a = 0;
while(a == 0)  
     
try
R((i+1):(i+12),:) = [];   % <--- This is horrible but it works 
i = i+12;
catch
    a = 1; % end of the current halving
end

end

end


end