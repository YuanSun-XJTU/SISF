function [W] = SISF(X,alpha,gamma,c,k20,m,lambda,A)
% X d*n data matrix
    Iter=20;
    eps_loss=0.000001;
    objvalue=[];
    N=size(X,2);
    d=size(X,1);
    S=rand(N);
    S=S./sum(S);
    S=(S+S')/2;
    Ls=diag(sum(S,2)) - S;
    W0 = orth(rand(d,m));
    for iter=1:Iter 
        %Update U
        U=UpdateU(alpha,lambda,Ls,c);
        V=X*(lambda*U'*Ls*U+alpha*(eye(N)-U)'*(eye(N)-U))*X';
        V=(V + V') / 2;
        eigvalue = eigs(V,1,'LA');
        V=eigvalue*eye(d)-V;
        %Update W
        W = IPU(V, m, k20,d,W0);
        %Update F
        F=U*X'*W; 
        obj=trace(-A'*S)+gamma*norm(S,'fro')^2+lambda*trace(F'*Ls*F)+alpha*norm(X'*W-F,'fro')^2;
        objvalue=[objvalue;obj]; 
        %Update S Ls 
        S=UpdateS(-A,F,gamma,lambda);
        Ls=diag(sum(S,2)) - S;
        if iter>1
            if abs(objvalue(iter)-objvalue(iter-1))/objvalue(iter-1)<eps_loss
                 break;
            end
        end
    end
end

