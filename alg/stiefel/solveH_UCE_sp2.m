function [ O,G ] = solveH_UCE_sp2( H,v,A )
%SOLVEYSTAR �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

B=bsxfun(@times,H,v);
C=sum(sum(H.*B));
O=C-2.*sum(sum(H.*A));


G=2.*B-2.*A;


end

