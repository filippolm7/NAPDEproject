function u = VCycle(A,f,u,l)

% Performs one loop of V-Cycle algorithm using Jacobi
% with the following parameters:

nu1 = 1; nu2 = 1;
w = 2/3;

if l == 1
    u = A\f;
else
    for i = 1:nu1
        u = u + w*(f-A*u)./diag(A);
    end
    r = f-A*u;
    R = RMatrix(size(A,1),1);
    Ah = R*A*R';
    e = VCycle(Ah,R*f,zeros(size(Ah,1),1),l-1);
    u = u + R'*e;
    
    for i = 1:nu2
        u = u + w*(f-A*u)./diag(A);
    end
end