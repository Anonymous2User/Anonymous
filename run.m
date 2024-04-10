clear
rand('state',0);

path='./base/';
DIRS=dir([path]);
nn=length(DIRS);
m=10;
for i=3
    name=DIRS(i).name
    DIRS2=dir([path name '/*.mat']);
    nn2=length(DIRS2);
        our_acc=zeros(nn2,m);
        our_time=zeros(nn2,m);
        our_nmi=zeros(nn2,m);
        for j=1:nn2
            load([path name '/' num2str(j)]);
            nSmp=length(yi{1});
            for iter=1:m
                for k=1:m
                    idx=(iter-1)*10+k;
                    Yi_input{k}=Yi{idx};
                end
                eta= 100;

                tic
                [H,Ri,Pi,alpha] = DACE(Yi_input,n_a,eta);
                tt=toc
                [Y,~]=mydiscretisation(H,1);
                [~,ypred]=max(Y,[],2);
                res=ClusteringMeasure(yi{1},ypred);
                our_acc(j,iter)=res(1);
                our_nmi(j,iter)=res(2);
                our_time(j,iter)=tt;
            end
        end
        acc =mean(our_acc');
        nmi =mean(our_nmi'); 
        time=mean(our_time');
        save(['./result_our/' name '.mat'], 'acc','nmi','time');
end