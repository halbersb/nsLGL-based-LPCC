function [dag,Observed,Latent,SOC] = runLPCC(clusters,clusters_unrounded,pcc_type,data,labels,major_clusters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find dag with latent exogenous and colliders.
% exogenous and colliders are indistinguishable yet
%
% input:
% [clusters]           - (matrix) integers: number of clusters -X- number of observed
% [clusters_unrounded] - (matrix) real: number of clusters -X- number of observed
% [pcc_type]           - (scalar) indicator for how to run pairwise cluster comparison (PCC)
% [data]               - (matrix) the original data set
% [labels]             - (vector) the cluster label for each sample in the data
% [major_clusters]     - (vector) indeces of the major clusters
%
% output:
% [dag]                - (matrix) the learned dag
% [Observed]           - (structure array) observed variables
% [Latent]             - (structure array) latents variables
% [SOC]                - (array) list of suspected observed colliders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% create major clusters matrix
if pcc_type>1
    m_clusters=clusters_unrounded(major_clusters,:);
    [mcc,~]=createPcc(m_clusters,pcc_type,data,labels,major_clusters); %create pcc for major clusters only (mcc stands for major clusters comparison)
else
    m_clusters=clusters(major_clusters,:);
    [mcc,~]=createPcc(m_clusters,pcc_type); %create pcc for major clusters only (mcc stands for major clusters comparison)
end
[Observed,Latent,~]=findPossibleLatent(mcc); %based on MSO - the returned set of latents includes both exogenous and colliders
[emcc]=extendPCCMatrix(mcc,Latent); %emcc stands for extended major clusters comparison
%disp(emcc); %print emcc
SOC=z_findSuspectedObservedColliders(Observed,Latent,emcc); %SOC stands for suspected observed colliders
[Latent]=findColliderLatents(Latent,emcc); %separate the latent set into latent exogenous and latent colliders
dag=createDag(Observed,Latent);