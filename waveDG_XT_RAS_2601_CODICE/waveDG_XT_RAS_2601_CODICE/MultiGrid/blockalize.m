function [Aj,bj] = blockalize(A,b)
    
    n = size(A,1);
    J = createJ(n);
    Aj = J*A*(J');
    bj = J*b;
    
end