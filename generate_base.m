clear;
path='./datasets/';
DIRS=dir([path,'*.mat']); 
nn=length(DIRS);
m=10;

%path2='./graphs/';
for ii=107:107
 %   try
    name=DIRS(ii).name
    load([path name]);
    name1=name(1:end-4);
    n=length(y);
    c=length(unique(y));
    idx=randperm(n);
    [aa,bb]=size(y);
    if aa==1
        y=y';
    end
    X=double(X(idx,:));
    y=y(idx,:);
    for iter=1:100
        ypred=litekmeans(X,c);
        res=ClusteringMeasure(y,ypred);
        acc(iter)=res(1);
        nmi(iter)=res(2);
        YY=zeros(n,c);
        for t=1:n
            YY(t,ypred(t))=1;
        end
        Yi{iter}=YY;
    end
    save(['./base_full/' name],'X','Yi','y','acc','nmi');
    mkdir(['./base/' name]);
    for ratio= 1:9
        n_a=floor(n*0.1*ratio);
        for iter=1:100
            if mod(iter,10)==1
                Xi{iter}=X;
                yi{iter}=y;
                idxi{iter}=1:n;
                ypred=litekmeans(Xi{iter},c);
                res=ClusteringMeasure(yi{iter},ypred);
                acc(iter)=res(1);
                nmi(iter)=res(2);
                YY=zeros(n,c);
                for t=1:n
                    YY(t,ypred(t))=1;
                end
                Yi{iter}=YY;
            else
                idx2=randperm(n-n_a);
                idx2=idx2+n_a;
                idx2=[1:n_a,idx2];
                Xi{iter}=X(idx2,:);
                yi{iter}=y(idx2,:);
                idxi{iter}=idx2;
                ypred=litekmeans(Xi{iter},c);
                res=ClusteringMeasure(yi{iter},ypred);
                acc(iter)=res(1);
                nmi(iter)=res(2);
                YY=zeros(n,c);
                for t=1:n
                    YY(t,ypred(t))=1;
                end
                Yi{iter}=YY;
            end
        end
        
        save(['./base/' name '/' num2str(ratio) '.mat'],'Yi','yi','idxi','acc','nmi','n_a');
    end
end
