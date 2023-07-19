function [uj,wj] = zoop(uj_wj)
% Extracts the solutions vector for plotting purposes

uj = [];
wj = [];
i = 1;
j = 1;

while(i< length(uj_wj))
    
    if rem(j,2) == 0
            wj = [wj;uj_wj(i:(i+5),1)];
    else
           uj = [uj;uj_wj(i:(i+5),1)]; 
    end
    i = i+6;
    j = j + 1;
    
end

end