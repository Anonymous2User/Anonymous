function [H,Ri,alpha] = UCE_full(Yi)
% function [H,Ri,Pi,alpha,acc_iter,nmi_iter] = UCE_sp4(Yi,n_a,eta,y_real)
%UCE 此处显示有关此函数的摘要
%   此处显示详细说明

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

%H=Yi{1};
alpha=ones(1,m)./m;
[H,Ri,Pi]=initialize_full(Yi);
for i=1:m-1
    % Pi{i}=eye(n,n);
    PY{i}=Pi{i}*Yi{i+1};
end
% [Y,R]=mydiscretisation(H,1);
% YR=Y*R';
% for i=1:m-1
%     Pi{i}=eye(n);
%     Ri{i}=eye(c);
%     PY{i}=Yi{i+1};
% end
% Ri{m}=eye(c);
% Y=H;
% R=eye(c);
% YR=Yi{1};
% a=ones(n-n_a,1);
% b=a;%%%
dd=zeros(1,m);
% objall=[];
acc_iter = zeros(1,maxiter);
nmi_iter = acc_iter;
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

   % objall=[objall obj(Yi,H,Ri,Pi,v,alpha,Y,R,lambda,gamma)];
   % objall=[objall obj(Yi,H,Ri,Pi,v,alpha,gamma)];
    %H
    A2=(alpha(1)^2).*Yi{1}*Ri{1};
    for i=2:m
       A2=A2+(alpha(i)^2).*PY{i-1}*Ri{i};
    end

    A=bsxfun(@times,A2,v.^2);
    % A=A+lambda.*YR;

    A=A./max(sum(alpha.^2),eps);
    %sum_alpha=sum(alpha.^2);
    [H]= OptStiefelGBB(H, @solveH_UCE_sp2, opts, v.^2,A);

    % objall=[objall obj(Yi,H,Ri,Pi,v,alpha,Y,R,lambda,gamma)];
% objall=[objall obj(Yi,H,Ri,Pi,v,alpha,gamma)];
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
    
    % objall=[objall obj(Yi,H,Ri,Pi,v,alpha,Y,R,lambda,gamma)];
% objall=[objall obj(Yi,H,Ri,Pi,v,alpha,gamma)];
    %Pi
    % VH = bsxfun(@times,H,v.^2);
    % for i=1:m-1
    %     M=-VH*Ri{i+1}';
    %     M2=M(n_a+1:n,:);
    %     Y2=Yi{i+1}';
    %     Y2=Y2(:,n_a+1:n);
    %     MY=M2*Y2;
        % VPY=bsxfun(@times, PY{i},v.^2);
        % % RR=Ri{1}*Ri{i+1}';
        % tmp2=(alpha(i+1)^2).*VPY;
        % tmp2=tmp2(n_a+1:n,:);
        % MY=MY+tmp2*Y2;

        % [~,P]=sinkhorn(a,b,MY, eta);
        % P = ipot_WD(a,b,MY,beta);
        % G=ones(n-n_a,n-n_a);
        % P = greenkhorn(MY,a,b,sum(a));
    % 
    %     Pi{i}=[eye(n_a),zeros(n_a,n-n_a);zeros(n-n_a,n_a),P];
    %     PY{i}=Pi{i}*Yi{i+1};
    % end
    % objall=[objall obj(Yi,H,Ri,Pi,v,alpha,Y,R,lambda,gamma)];
% objall=[objall obj(Yi,H,Ri,Pi,v,alpha,gamma)];
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

% objall=[objall obj(Yi,H,Ri,Pi,v,alpha,Y,R,lambda,gamma)];
 % objall=[objall obj(Yi,H,Ri,Pi,v,alpha,gamma)];
    %Y,R
    % [Y,R]=mydiscretisation(H,0,R);
    % YR=Y*R';
% objall=[objall obj(Yi,H,Ri,Pi,v,alpha,Y,R,lambda,gamma)];
    gamma=gamma*1.1;
    % objall=[objall obj(Yi,H,Ri,Pi,v,alpha,gamma)];
%    objall=[objall obj(Yi,H,Ri,Pi,alpha,YR,lambda)];


% [Y,~]=mydiscretisation(H,1);
% [~,y] = max(Y,[],2);
% res=ClusteringMeasure(y_real,y);
% acc_iter(iter) = res(1);
% nmi_iter(iter) = res(2);

end
% figure('visible','on');
% plot(objall)
%objall=reshape(objall,[6,20]);
%objall=objall'
%plot(objall)


% for i=2:20
%     PPP=Pi{i-1};
%     [~,iidx]=max(PPP,[],2);
%     iidx=iidx(n_a+1:n);
%     PPP=zeros(3,3);
%     iidx=iidx-n_a;
%     for jj=1:3
%         PPP(jj,iidx(jj))=1;
%     end
%     idxit=idxi{i};
%     idxx=idxit(203:-1:201)';
%     [PPP*idxx,[203;202;201]]
% end

% for i=2:2
%     PPP=Pi{i-1};
%     [~,iidx]=max(PPP,[],2);
%     iidx=iidx(n_a+1:n);
%     PPP=zeros(4,4);
%     iidx=iidx-n_a;
%     for jj=1:4
%         PPP(jj,iidx(jj))=1;
%     end
%     idxit=idxi{i};
%     idxx=idxit(9:12)';
%     [PPP*idxx,[9;10;11;12]]
% end
