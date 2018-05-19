

%% add toolbox of threshold schemes

addpath(genpath('C:\Users\mpnsd\Desktop\SOFTWARE\threshold_schemes\threshold_schemes\threshold_schemes'))

%% create a number of random full-weighted brain networks

no=10;
rois=90;
toy=rand(10,rois,rois);

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

%% take the sum of rows
sum1=sum(dist);

%normalize and keep the coefficients of the linear combination of
%individualized brain networks
coef=sum1./sum(sum1);

%%integrated brain network
integr=zeros(90,90);

for k=1:no
    integr=integr + coef(k).*squeeze(thresholded(k,:,:));
end

%% normalize the integrated brain network
integr=integr./max(integr);

%% topological filtering of integrated brain network with OMST
 [nCIJtree CIJtree mdeg  globalcosteffmax costmax E]=threshold_omst_gce_wu(integr,1);
 
 figure(2),imagesc(coef) ; colorbar
           title('Coefficients of the linear combination of brain networks')



