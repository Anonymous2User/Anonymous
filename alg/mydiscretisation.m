function [EigenvectorsDiscrete,R]=mydiscretisation(EigenVectors,flag,R)


[n,k]=size(EigenVectors);

if flag==1
    vm = max(sqrt(sum(EigenVectors.*EigenVectors,2)),eps);
    EigenVectors = EigenVectors./repmat(vm,1,k);
end

if nargin<3
    R=zeros(k);
    R(:,1)=EigenVectors(1+round(rand(1)*(n-1)),:)';
    c=zeros(n,1);
    for j=2:k
        c=c+abs(EigenVectors*R(:,j-1));
        [minimum,i]=min(c);
        R(:,j)=EigenVectors(i,:)';
    end
end

lastObjectiveValue=0;
exitLoop=0;
nbIterationsDiscretisation = 0;
nbIterationsDiscretisationMax = 20;%voir
while exitLoop== 0 
    nbIterationsDiscretisation = nbIterationsDiscretisation + 1 ;   
    EigenvectorsDiscrete = discretisationEigenVectorData(EigenVectors*R);
    [U,S,V] = svd(EigenvectorsDiscrete'*EigenVectors,0);    
    NcutValue=2*(n-trace(S));
    
    if abs(NcutValue-lastObjectiveValue) < eps | nbIterationsDiscretisation > nbIterationsDiscretisationMax
        exitLoop=1;
    else
        lastObjectiveValue = NcutValue;
        R=V*U';
    end
end