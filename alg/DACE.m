function [H,Ri,Pi,alpha] = DACE(Yi,n_a,eta)
opts.record = 0;
opts.mxitr  = 200;
opts.xtol = 1e-5;
opts.gtol = 1e-5;
opts.ftol = 1e-8;
opts.tau = 1e-3;

options=[];
options.Display='none';

m=length(Yi);
[n,c]=size(Yi{1});
maxiter=20;

alpha=ones(1,m)./m;
[H,Ri,Pi]=initialize(Yi,n_a,eta);
for i=1:m-1
    PY{i}=Pi{i}*Yi{i+1};
end
a=ones(n-n_a,1);
b=a;
dd=zeros(1,m);

for iter=1:maxiter
    %V
    A1=Yi{1}*Ri{1};
    a1=2*(alpha(1)^2).*sum((H-A1).^2,2);
    for i=2:m
        A1=PY{i-1}*Ri{i};
        a1=a1 + 2*(alpha(i)^2).*sum((H-A1).^2,2);
    end
    if iter==1
        A1_sort=sort(a1,'ascend');
        gamma=A1_sort(floor(n*0.1));
    end
    v=gamma./a1;
    v=min(v,1);

    %H
    A2=(alpha(1)^2).*Yi{1}*Ri{1};
    for i=2:m
       A2=A2+(alpha(i)^2).*PY{i-1}*Ri{i};
    end

    A=bsxfun(@times,A2,v.^2);
    A=A./max(sum(alpha.^2),eps);
    [H]= OptStiefelGBB(H, @solveH_UCE_sp2, opts, v.^2,A);

    %Ri
    VH=bsxfun(@times,H,v.^2);
    B=Yi{1}'*VH;
    [U2,~,V2]=svd(B,'econ');
    Ri{1}=U2*V2';
    for i=2:m 
        B=PY{i-1}'*VH;
        [U2,~,V2]=svd(B,'econ');
        Ri{i}=U2*V2';
    end
    
    %Pi
    VH = bsxfun(@times,H,v.^2);
    for i=1:m-1
        M=-VH*Ri{i+1}';
        M2=M(n_a+1:n,:);
        Y2=Yi{i+1}';
        Y2=Y2(:,n_a+1:n);
        MY=M2*Y2;
        [~,P]=sinkhorn(a,b,MY, eta);
        Pi{i}=[eye(n_a),zeros(n_a,n-n_a);zeros(n-n_a,n_a),P];
        PY{i}=Pi{i}*Yi{i+1};
    end

    %alpha
    D = zeros(1,m);
    A3 = Yi{1}*Ri{1};
    B3 = bsxfun(@times,v,H-A3);
    D(1) = 1/sum(sum(B3.^2));
    
    for i = 1:m-1
        A3 = PY{i}*Ri{i+1};
        B3 = bsxfun(@times,v,H-A3);
        D(i+1) = 1/sum(sum(B3.^2));
    end
    alpha =D./sum(D);

    gamma=gamma*1.1;
end

