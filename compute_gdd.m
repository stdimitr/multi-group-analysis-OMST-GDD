% compute_gdd
%
% function [gdd_star,t_star,t_upperbound]=compute_gdd(A1,A2)
%
% Computes the graph diffusion distance between two connectivity
% matrices
%
% Inputs : 
% A1, A2 - symmetric NxN matrices representing the edge weights
%
% Outputs
%
% gdd_star - the graph diffusion distance
% t_star - the diffusion time corresponding to the observed maximum
% t_upperbound - largest diffusion time probed during optimization
%
% 
%

%-DISCLAIMER---------------------------------------------------------------
% This code has been implemented based on the GlobalSIP'13 conference
% paper. You can use this source code for non commercial research and
% educational purposes only, without licensing fees and is provided without
% guarantee or warrantee expressed or implied. You cannot repost this file
% without prior written permission from the authors. If you use this code
% please cite the following paper:

% David K Hammond, Yaniv Gur, Chris R Johnson, "Graph Diffusion Distance:
% A Difference Measure for Weighted Graphs Based on the Graph Laplacian
% Exponential Kernel", in proceedings of the 2013 IEEE Global Conference on
% Signal and Information Processing (GlobalSIP), December 3-5, 2013,
% Austin, Texas, USA.

% Code author: Dimitriadis Stavros I
%https://github.com/stdimitr/multi-group-analysis-OMST-GDD
% https://www.researchgate.net/profile/Stavros_Dimitriadis

function [gdd_star,t_star,t_upperbound]=compute_gdd(A1,A2)
    L1=full(graph_laplacian(A1));
    L2=full(graph_laplacian(A2));
    [V1,D1]=eig(L1);
    [V2,D2]=eig(L2);
    
    [gdd_star,t_star,t_upperbound]=compute_gdd_pd(V1,D1,V2,D2);
end

% Internal functions
% _pd : prediagonalized : spectral decomposition of A1, A2 already computed
function [gdd_star,t_star,t_upperbound]=compute_gdd_pd(V1,D1,V2,D2)
    alleigval=[diag(D1);diag(D2)];
    nonzeroeigval=alleigval(alleigval>1e-14); % hard coded threshold here
    nonzeroeigval=sort(nonzeroeigval);
    t_upperbound=1/nonzeroeigval(1); % reciprocal of minimum nonzero eigenvalue
    
    m_dmf=@(t) -1*gdd_xi_t(V1,D1,V2,D2,t) ;
    %options=optimset('Display','iter','TolX',1e-4);
    options=optimset('Display','off','TolX',1e-4);
    
    [x,fval,exitflag]=fminbnd(m_dmf,0,t_upperbound,options);
    gdd_star=sqrt(-fval);
    t_star=x;
    
end
% gdd_xi_t : squared graph diffusion distance, at scale t
%
% Computes frobenius norm of difference of laplacian exponential diffusion kernels,
% at specified timepoints
%
% Vn,Dn : eigenvector/eigenvalue decomposition of Ln, so Vn*Dn*Vn'=Ln, for
% n=1,2 
% t - array of times to compute gdd function
% Outputs :
% E - same size as t, contains differences of frobenius norms
function E=gdd_xi_t(V1,D1,V2,D2,t)
    E=zeros(size(t));
    
    for kt=1:numel(t)
        
% $$$ tmp=expm(-t*L1)-expm(-t*L2);    
% $$$ tmp=V1*diag(exp(-t(kt)*diag(D1)))*V1' - ...
% $$$     V2*diag(exp(-t(kt)*diag(D2)))*V2';  
% code below is equivalent to code above, but faster as
% it avoids matrix-matrix multiply for multiplying by
% diagonal matrix
        
        ed1=exp(-t(kt)*diag(D1));
        ed2=exp(-t(kt)*diag(D2));
        tmp=V1*bsxfun(@times,ed1(:),V1')-...
            V2*bsxfun(@times,ed2(:),V2');
        
        E(kt)=sum(tmp(:).^2);
    end
    
end

function L = graph_laplacian(A)
    L=diag(sum(A))-A;
end
