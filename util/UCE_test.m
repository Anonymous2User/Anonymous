function [Y,H,Ri,Pi,alpha] = UCE_sp2(Yi,n_a,lambda,eta,idxi)
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
[H,Ri,Pi]=initialize(Yi,n_a,eta);
for i=1:m-1
    PY{i}=Pi{i}*Yi{i+1};
end
[Y,R]=mydiscretisation(H,1);
YR=Y*R';
% for i=1:m-1
%     Pi{i}=eye(n);
%     Ri{i}=eye(c);
%     PY{i}=Yi{i+1};
% end
% Ri{m}=eye(c);
% Y=H;
% R=eye(c);
% YR=Yi{1};
a=ones(n-n_a,1);
b=a;
dd=zeros(1,m);
objall=[];
for iter=1:maxiter
    %V
    A1=alpha(1).*Yi{1}*Ri{1};
    for i=2:m
        A1=A1+alpha(i)*PY{i-1}*Ri{i};
    end
    a1=2.*sum((H-A1).^2,2);
    if iter==1
        A1_sort=sort(a1,'ascend');
        gamma=A1_sort(floor(n*0.1));
    end
    v=gamma./a1;
    v=min(v,1);

  %  objall=[objall obj(Yi,H,Ri,Pi,v,alpha,Y,R,lambda,gamma)];
    %H
    %A=alpha(1).*Yi{1}*Ri{1};
    %for i=2:m
    %    A=A+alpha(i).*PY{i-1}*Ri{i};
    %end
    A=bsxfun(@times,A1,v.^2);
    A=A+lambda.*YR;

    %sum_alpha=sum(alpha.^2);
    [H]= OptStiefelGBB(H, @solveH_UCE_sp2, opts, v.^2,A);

 %    objall=[objall obj(Yi,H,Ri,Pi,v,alpha,Y,R,lambda,gamma)];
%    objall=[objall obj(Yi,H,Ri,Pi,alpha,YR,lambda)];
    %Ri
    
    VH=bsxfun(@times,H,v.^2);
    % for i=1:m
    %     if i==1
    %         B=Yi{1}'*VH;
    %         for j=2:m
    %             VPY=bsxfun(@times, PY{j-1},v.^2);
    %             B=B-alpha(j).*Yi{1}'*VPY*Ri{j};
    %         end
    %         [U2,~,V2]=svd(B,'econ');
    %         Ri{1}=U2*V2';
    %     else
    %         VPY=bsxfun(@times, Yi{1},v.^2);
    %         B=PY{i-1}'*VH-alpha(1).*PY{i-1}'*VPY*Ri{1};
    %         for j=2:m
    %             if j==i
    %                 continue;
    %             else
    %                 VPY=bsxfun(@times, PY{j-1},v.^2);
    %                 B=B-alpha(j).*PY{i-1}'*VPY*Ri{j};
    %             end
    %         end
    %         [U2,~,V2]=svd(B,'econ');
    %         Ri{i}=U2*V2';
    %     end
    % end
    
%     objall=[objall obj(Yi,H,Ri,Pi,v,alpha,Y,R,lambda,gamma)];
%    objall=[objall obj(Yi,H,Ri,Pi,alpha,YR,lambda)];
    %Pi
    for i=1:m-1
        M=-VH*Ri{i+1}';
        M2=M(n_a+1:n,:);
        Y2=Yi{i+1}';
        Y2=Y2(:,n_a+1:n);
        MY=M2*Y2;
        VPY=bsxfun(@times, Yi{1},v.^2);
        RR=Ri{1}*Ri{i+1}';
        tmp2=alpha(1).*VPY*RR;
        tmp2=tmp2(n_a+1:n,:);
        MY=MY+tmp2*Y2;

        for j=1:m-1
            if j==i
                continue;
            else
                VPY=bsxfun(@times,PY{j},v.^2);
                RR=Ri{j+1}*Ri{i+1}';
                tmp3=alpha(j+1).*VPY*RR;
                tmp3=tmp3(n_a+1:n,:);
                MY=MY+tmp3*Y2;
            end
        end

        [~,P]=sinkhorn(a,b,MY, eta);
        Pi{i}=[eye(n_a),zeros(n_a,n-n_a);zeros(n-n_a,n_a),P];
        PY{i}=Pi{i}*Yi{i+1};
    end
 %    objall=[objall obj(Yi,H,Ri,Pi,v,alpha,Y,R,lambda,gamma)];
 %   objall=[objall obj(Yi,H,Ri,Pi,alpha,YR,lambda)];
    %alpha
    Q=zeros(m,m);
    RY=Yi{1}*Ri{1};
    VRY=bsxfun(@times,RY,v.^2);
    Q(1,1)=sum(sum(RY.*VRY));
    for i=2:m
        VRY=PY{i-1}*Ri{i};
        VRY=bsxfun(@times,VRY,v.^2);
        Q(1,i)=sum(sum(RY.*VRY));
        Q(i,1)=Q(1,i);
    end
    for i=2:m
        for j=i:m
            VRY1=PY{i-1}*Ri{i};
            VRY2=PY{j-1}*Ri{j};
            VRY1=bsxfun(@times,VRY1,v.^2);
            Q(i,j)=sum(sum(VRY1.*VRY2));
            Q(j,i)=Q(i,j);
        end
    end
    f=zeros(m,1);
    f(1)=-sum(sum((VH'*Yi{1}).*Ri{1}'));
    for i=2:m
        PYR=PY{i-1}*Ri{i};
        f(i)=-sum(sum(VH.*PYR));
    end

    alpha=quadprog(Q,f,[],[],ones(1,m),1,zeros(m,1),ones(m,1),alpha,options);


% objall=[objall obj(Yi,H,Ri,Pi,v,alpha,Y,R,lambda,gamma)];
 %   objall=[objall obj(Yi,H,Ri,Pi,alpha,YR,lambda)];
    %Y,R
    [Y,R]=mydiscretisation(H,0,R);
    YR=Y*R';
% objall=[objall obj(Yi,H,Ri,Pi,v,alpha,Y,R,lambda,gamma)];
    gamma=gamma*1.1;
%    objall=[objall obj(Yi,H,Ri,Pi,alpha,YR,lambda)];
end
%plot(objall)
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
