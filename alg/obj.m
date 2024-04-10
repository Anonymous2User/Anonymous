function obj1 = obj(Yi,H,Ri,Pi,v,alpha,gamma)
%OBJ 此处显示有关此函数的摘要
%   此处显示详细说明
A=Yi{1}*Ri{1};
B=bsxfun(@times,v,H-A);
obj1=(alpha(1)^2)*sum(sum(B.^2));
m=length(alpha);
for i=2:m
    A=Pi{i-1}*Yi{i}*Ri{i};
    B=bsxfun(@times,v,H-A);
    obj1=obj1+(alpha(i)^2)*sum(sum(B.^2));
end

% obj1=obj1+lambda*sum(sum((Y-H*R).^2))-gamma.*sum(v);
obj1=obj1-gamma.*sum(v);


% obj1=alpha(1).^2.*sum(sum((H-Yi{1}).^2));
% m=length(Yi);
% for i=1:m-1
%     obj1=obj1+alpha(i+1).^2.*sum(sum((H-Pi{i}*Yi{i+1}*Ri{i}).^2));
% end
% obj1=obj1+lambda.*sum(sum((Y-H*R).^2));
end

