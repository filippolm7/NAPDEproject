function [R] = newRMatrix(n,k)

% Creates a restriction matrix for FEM matrices with 2nd order polynomials
% These are made of 6x6 block matrices, in pairs of 2 (uh,wh)
% n = matrix dimension
% k = number of halvings of the elements (0 returns an identity)

e = ones(n,1);

R = spdiags([0.5*e e 0.5*e],[-12,0,12],n,n);

if k == 0
    return
else
    
for j = 1:k % loop over n. of halvings
  i = 0;
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