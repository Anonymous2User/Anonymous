function [D,P]=sinkhorn(a,b,M, lambda)
K=exp(-lambda*M);
K(K<1e-100)=1e-100;
U=K.*M;

I=(a>0);
if ~all(I)
    K=K(I,:);
    U=U(I,:);
    a=a(I);
end

ainvK=bsxfun(@rdivide,K,a);
i=0;
u=ones(size(a,1),size(b,2))/size(a,1);

% iterations
while i<1000
        u=1./(ainvK*(b./(K'*u)));
    i=i+1;

    if mod(i,20)==1 || i==1000     
        v=b./(K'*u);
        u=1./(ainvK*v);
                        
        % check whether to stop
        stop=norm(sum(abs(v.*(K'*u)-b)),inf);
        if stop< .5e-2|| isnan(stop) % norm of all || . ||_1 differences between the marginal of the current solution with the actual marginals.
            break;
        end      
        i=i+1;
    end
end

D=sum(u.*(U*v));
%P=diag(u)*K*diag(v);
P1=bsxfun(@times,K,u);  
P =bsxfun(@times,P1,v');