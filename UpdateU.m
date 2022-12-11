function [U] = UpdateU(alpha,lambda,Ls,c)
       N=size(Ls,1);
       temp=lambda*Ls+alpha*eye(N);
       [P,Lambda,Q] = svd(temp);
       Lambda=inv(Lambda);
       Lambda(1:c,1:c)=0;
       U=alpha*Q*Lambda*P';    
end

