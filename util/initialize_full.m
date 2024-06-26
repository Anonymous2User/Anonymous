function [H,Ri,Pi2] = initialize_full(Yi)
%INITIALIZE 此处显示有关此函数的摘要
%   此处显示详细说明
m=length(Yi);
[n,c]=size(Yi{1});

% a=ones(n-n_a,1);
% b=a;%%%

% for i=1:m
%     YY=Yi{i};
%     Yi_a{i}=YY(1:n_a,:);
%     Yi_u{i}=YY(n_a+1:n,:);
% end

KC=zeros(n,n);
for i=1:m
    KC=KC+Yi{i}*Yi{i}';
    % KC=KC+Yi_a{i}*Yi_a{i}';
end
KC=KC./m;
H = mykernelkmeans(KC,c);

maxiter=10;

for j=1:m
    Ri{j}=eye(c);
end

for i=1:maxiter
    
    H1=zeros(n,c);
    for j=1:m
        A=Yi{j}'*H;
        for k=1:m
            if k==j
                continue;
            else
                A=A-1/m.*Yi{j}'*Yi{k}*Ri{k};
            end
        end
        [U2,~,V2]=svd(A,'econ');
        Ri{j}=U2*V2';
        H1=H1+Yi{j}*Ri{j};
    end
    H1=H1./m;

    [U1,~,V1]=svd(H1,'econ');
    H=U1*V1';
end


% H_u=Yi_u{1}*Ri{1};
% 
% for j=1:m
%     YR{j}=Yi_u{j}*Ri{j};
% end
for j=1:m-1
    Pi2{j}=eye(n);
end
% for i=1:maxiter
%     H2=1./m.*Yi_u{1}*Ri{1};
%     A=H_u-H2;
% 
%     for j=1:m-1
%         B=-A*YR{j+1}';
%         for k=1:m-1
%             if k==j
%                 continue;
%             else
%                 B=B+1./m.*Pi{k}*YR{k+1}*YR{j+1}';
%             end
%         end
%         % [~,P]=sinkhorn(a,b,B, eta);
%         % P = ipot_WD(a,b,B,beta);
%         % G=ones(n-n_a,n-n_a);
%         % P = greenkhorn(B,a,b,sum(a));
%         % Pi{j}=P;
%         H2=H2+1./m.*Pi{j}*YR{j+1};
%     end
%     [U3,~,V3]=svd(H2,'econ');
%     H_u=U3*V3';
% end
    
% [U4,~,V4]=svd(H_all,'econ');
% H=U4*V4';

% for j=1:m-1
%     Pi2{j}=[eye(n_a),zeros(n_a,n-n_a);zeros(n-n_a,n_a),Pi{j}];
% end

end

