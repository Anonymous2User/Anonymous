A=diag([1,2,3,4]);
B=[0,1,0,0;1,0,0,0;0,0,0,1;0,0,1,0];
M=-B*A';

a=ones(4 ,1);
b=a;
% lambda=10
% [D,P]=sinkhorn(a,b,M, lambda)
beta = 2;
[P,loss] = ipot_WD(a,b,M,beta)