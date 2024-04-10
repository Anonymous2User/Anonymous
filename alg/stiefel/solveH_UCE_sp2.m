function [ O,G ] = solveH_UCE_sp2( H,v,A )
%SOLVEYSTAR 此处显示有关此函数的摘要
%   此处显示详细说明

B=bsxfun(@times,H,v);
C=sum(sum(H.*B));
O=C-2.*sum(sum(H.*A));


G=2.*B-2.*A;


end

