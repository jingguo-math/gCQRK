function [S1] = idstf(rhs,J,s)

   S1=[];
   for i=1:s
   S=transpose(idst(transpose(idst(reshape(rhs(:,i),J-1,J-1)))));
   S1=[S1 reshape(S,(J-1)^2,1)];
   end
end