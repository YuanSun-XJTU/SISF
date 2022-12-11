function [S] = UpdateS(D,F,gamma,lambda)
    N=size(D,1);
    d_hat = L2_distance_1(F',F');
    
    S= zeros(N);
    for i=1:N

        idxa0 = 1:N;
        dfi = d_hat(i,idxa0);
        dxi = D(i,idxa0);
        ad = -(dxi+lambda*dfi)/(2*gamma);
        S(i,idxa0) = EProjSimplex_new(ad);
    end 
    S=(S+S')/2;   
end

