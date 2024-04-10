function [ O,G ] = solveH( H,v,A,sum_alpha )
%SOLVEYSTAR 此处显示有关此函数的摘要
%   此处显示详细说明

B=bsxfun(@times,H,v);
C=sum(sum(H.*B)).*sum_alpha;
O=C-2.*sum(sum(H.*A));


G=2.*sum_alpha.*B-2.*A;


end

