function [Y2_new] = hungarian_align(Y1,Y2,n_a,if_pca)
%HUNGARIAN_ALIGN 此处显示有关此函数的摘要
%   此处显示详细说明
[n,c1]=size(Y1);
[n,c2]=size(Y2);
if if_pca==1
    Y11=Y1(1:n_a,:);
    Y21=Y2(1:n_a,:);
    coeff1=pca(Y11);
    coeff2=pca(Y21);
    Y12=Y1(n_a+1:n,:);
    Y12=bsxfun(@minus,Y12,mean(Y12,1));
    Y22=Y2(n_a+1:n,:);
    Y22=bsxfun(@minus,Y22,mean(Y22,1));
    X1=Y12*coeff1(:,1:2);
    X2=Y22*coeff2(:,1:2);
else
    X1=Y1(n_a+1:n,:);
    X2=Y2(n_a+1:n,:);
end
options.KernelType='Gaussian';
options.t=optSigma2(X1,X2);
K = constructKernel(X1,X2,options);
if n-n_a>5
    k=5;
else
    k=n-n_a-1;
end
G=knn_new(K,k);
[C,T]=hungarian(-G);
Y22_ori=Y2(n_a+1:n,:);
Y21_ori=Y2(1:n_a,:);
Y22_new = zeros(size(Y22_ori));

for i=n_a+1:n
    Y22_new(i-n_a,:)=Y22_ori(C(i-n_a),:);
end
Y2_new=[Y21_ori;Y22_new];

