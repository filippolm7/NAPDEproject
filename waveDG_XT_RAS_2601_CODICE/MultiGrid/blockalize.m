function [Aj,bj] = blockalize(A,b,n)
    % This function creates the block matrices from A and b, performing an
    % elimination of half the rows (n-times)

    k = size(A,1);
    J = createJ(k);
    
    
    J =  coarserJ(J,n);
    
    % Blockalizes the system matrices
    Aj = J*A*(J');
    bj = J*b;
    
end