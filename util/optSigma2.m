function sigma = optSigma(X,Y)
    N = size(X,1);
    dist = EuDist2(X,Y); 
    dist = reshape(dist,1,N*N);
    sigma = median(dist);
     if sigma==0
         dist2=dist(find(dist>0));
         sigma=median(dist2);
     end