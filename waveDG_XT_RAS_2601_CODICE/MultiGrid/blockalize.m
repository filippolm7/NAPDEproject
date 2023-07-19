function [Aj,bj] = blockalize(A,b)

    % Code to pass from a (uh1,uh2,...,uhm,wh1,wh2,...,whm) system into a 
    % (uh1,wh1,uh2,wh2,...,uhm,whm) system

    k = size(A,1);
    J = createJ(k);

    % Blockalizes the system matrices
    Aj = J*A*(J');
    bj = J*b;
    
end