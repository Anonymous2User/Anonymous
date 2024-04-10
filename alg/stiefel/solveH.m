function [ O,G ] = solveH( H,v,A,sum_alpha )
%SOLVEYSTAR �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

B=bsxfun(@times,H,v);
C=sum(sum(H.*B)).*sum_alpha;
O=C-2.*sum(sum(H.*A));


G=2.*sum_alpha.*B-2.*A;


end

