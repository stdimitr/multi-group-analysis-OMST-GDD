

%% add toolbox of threshold schemes

addpath(genpath('C:\Users\mpnsd\Desktop\SOFTWARE\threshold_schemes\threshold_schemes\threshold_schemes'))

%% create a number of random full-weighted brain networks

no=10;
rois=90;
toy=rand(10,rois,rois);

%% transform them to undirected
for k=1:rois
    for l=(k+1):rois
        toy(:,l,k)=toy(:,k,l);
    end
end

thresholded=zeros(no,rois,rois);

%%% topologically filtered of brain networks with OMST

for k=1:no
    [nCIJtree CIJtree mdeg  globalcosteffmax costmax E]=threshold_omst_gce_wu(squeeze(toy(k,:,:)),0);
    thresholded(k,:,:)=CIJtree;
end

%%distance matrix with graph diffusion distance metric

dist=zeros(no,no);

for k=1:no
    for l=(k+1):no
        A1=squeeze(thresholded(k,:,:));
        A2=squeeze(thresholded(l,:,:));
        [gdd,t,t_upperbound]=compute_gdd(A1,A2);
        dist(k,l)=gdd;
        dist(l,k)=dist(k,l);
    end
end

%% dist can be projected to a 2D space with e.g. multi-dimensional scaling



