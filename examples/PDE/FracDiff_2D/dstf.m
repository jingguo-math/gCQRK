function [S1] = dstf(rhs,J,s)

   S1=[];
   for i=1:s
   S=transpose(dst(transpose(dst(reshape(rhs(:,i),J-1,J-1)))));
   S1=[S1 reshape(S,(J-1)^2,1)];
   end
end