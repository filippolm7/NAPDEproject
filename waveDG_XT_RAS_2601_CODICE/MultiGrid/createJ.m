function [J]= createJ(n)

% Code to create the matrix necessary to make A a block matrix, i.e. to
% pass from a (uh1,uh2,...,uhm,wh1,wh2,...,whm) system into a
% (uh1,wh1,uh2,wh2,...,uhm,whm) system

% Made considering FEM matrices with 2nd order polynomials, made of 6x6
% blocks

if rem(n,2) ~= 0
    disp("Error, input must be an even number")
    return
end

m = n/2;

j1 = 1;
j2 = 1;
J = sparse(n,n);
i = 1;
it = 1;

while(i <= n) % righe 
   
      if rem(it,2) == 0
            for k = 0:5                
                J(i+k,j1+m+k) = 1;            
            end
            j1 = j1+6;
      else
          
          for k = 0:5              
             J(i+k,j2+k) = 1;  
          end
          j2 = j2+6;
      end
      it = it+1;
      i = i+6;
end
    
end


