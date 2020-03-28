function [Observed,Latent_final,NUM_PATHS,PATH_Len] = latentSplit2_pr3(clusters,clusters_unrounded,try_big,major_cluster,l,Observed,Latent,l_ch,pcc_type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find split for current latent based on pairwise comparison of
% first order minor cluster (FMC) to major clusters
%
% input:
% [clusters]           - (matrix) the clusters centers
% [clusters_unrounded] - (matrix) real: number of clusters -X- number of observed
% [try_big]            - (vector) subset of minor clusters (sms)
% [major_cluster]      - (scalar) current major cluster to compare to
% [l]                  - (scalar) current latent
% [Observed]           - (structure array) observed variables
% [Latent]             - (structure array) latents variables
% [l_ch]               - (vector) list of current latent children
% [pcc_type]           - [scalar] indicator for how to run pairwise cluster comparison
%
% output:
% [Observed]           - (structure array) updated observed variables
% [Latent_final]       - (structure array) new latents variables
% [NUM_PATHS]          - (scalar)
% [PATH_Len]           - (scalar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
if pcc_type>1, clusters=clusters_unrounded; end %if pcc_type is not regular use the unrounded clusters
FMC=clusters(try_big,:); %FMC = first order minor clusters
m_clusters=clusters(major_cluster,:); %current major cluster (to compare to sms)
Latent_final=[];
if isempty(FMC)
    Latent_final=struct([]);
    Obsreved=struct([]);
    NUM_PATHS=0;
    PATH_Len=0;
    return
end

%step 1: create pairwise cluster compartion between all clusters in FMC and all cluster in major clusters 
[PCCS,~]=createPccMM(FMC,m_clusters,pcc_type); %was: m_clusters(2,:) - why 2??? 
PCCS_M=PCCS(:,Latent(l).OCH); %take only the children of Latent(l) into consideration (columns selection)
PCCS_M=unique(PCCS_M,'rows'); %dan added: remove duplicate rows

%step 2: filter out PCC that have exctly one '0' or exctly one '1'
PCCS_2S=[];
for m=1:size(PCCS_M,1)
    %if only one observed changes value filter it out:
    %because latent must have at least two observed
    if (sum(PCCS_M(m,:))==length(Latent(l).OCH)-1) || (sum(PCCS_M(m,:))==1) %good engouh for real problem??
        %only one element in the row is "0" or only one is "1"
        continue;
    end
    PCCS_2S=[PCCS_2S;(PCCS_M(m,:))];
end
if isempty(PCCS_2S) %if PCCS_2S is empty - break
    Latent_final=struct([]);
    Obsreved=struct([]);
    NUM_PATHS=0;
    PATH_Len=0;
    return;
end

%step 3: split current latent
[Observed,Latent_final,~]=findPossibleLatentSplit3(PCCS_2S,l,Observed,l_ch);
if isempty(Latent_final) %if no possible splits were found - break
    Latent_final=struct([]);
    Obsreved=struct([]);
    NUM_PATHS=0;
    PATH_Len=0;
    return;
end

%step 4: create latent tag (L_tag) which is a copy of Latent_final that indicates on number of paths
%Latent_final(i) will be added if at least in one row of PCCS_2S the
%latent's LF_CH_IN are all '1's and they are the only '1's in that row
%also updated '.CH_2S'
L_tag=[];
for i=1:length(Latent_final)
    LF_CH_IN=[];
    for c=1:length(Latent_final(i).CH)
       ind=find(l_ch==Latent_final(i).CH(c));
       LF_CH_IN=[LF_CH_IN ind];
    end
    Latent_final(i).CH_2S=LF_CH_IN;
    for j=1:size(PCCS_2S,1)
        %if only LF_CH_IN are '1's in the row
        if sum(PCCS_2S(j,:))==length(Latent_final(i).CH) && sum(PCCS_2S(j,Latent_final(i).CH_2S))==length(Latent_final(i).CH)
            flag=1;
            for x=1:length(L_tag)
                if L_tag(x).PCC_index==Latent_final(i).PCC_index
                    flag=0;
                    break;
                end
            end
            if flag==1
                L_tag=[L_tag Latent_final(i)];
                continue;
            end
        end
    end
end
if isempty(L_tag) %if no L_tag were found - break
    Latent_final=struct([]);
    Obsreved=struct([]);
    NUM_PATHS=0;
    PATH_Len=0;
    return;
end

%step 5: update Latent_final and calculate NUM_PATHS
PCCS_2S_tag={};
for x=1:length(L_tag)
    PCCS_2S_tag{x}=[];
    for j=1:size(PCCS_2S,1)
        if sum(PCCS_2S(j,L_tag(x).CH_2S))==0 %length(L_tag(x).CH)
            PCCS_2S_tag{x}=[PCCS_2S_tag{x};PCCS_2S(j,:)];
        end
    end
    PCCS_2S_tag{x}=unique(PCCS_2S_tag{x},'rows'); %remove duplications
    for i=1:length(Latent_final)
        Latent_final(i).K{x}=0;
        Latent_final(i).W{x}=[];
        for j=1:size(PCCS_2S_tag{x},1)
            if sum(PCCS_2S_tag{x}(j,Latent_final(i).CH_2S))==length(Latent_final(i).CH)
                Latent_final(i).K{x}=Latent_final(i).K{x}+1;
                Latent_final(i).W{x}=[Latent_final(i).W{x} j];
            end
        end
    end
    PATH_Len{x}=0;
    for i=1:length(Latent_final)
        Latent_final(i).PATH{x}=[];
        for j=1:length(Latent_final)
            flag=1;
            for w=1:length(Latent_final(i).W{x})
                if ~(sum(PCCS_2S_tag{x}(Latent_final(i).W{x}(w),Latent_final(i).CH_2S))==length(Latent_final(i).CH) && sum(PCCS_2S_tag{x}(Latent_final(i).W{x}(w),Latent_final(j).CH_2S))==length(Latent_final(j).CH))
                    flag=0;
                end
            end
            if (Latent_final(i).K{x}==Latent_final(j).K{x}-1) && flag==1
                Latent_final(i).PATH{x}=[Latent_final(i).PATH{x} Latent_final(j)];
                PATH_Len{x}=PATH_Len{x}+1;
            end
        end
    end
end
%NUM_PATHS>2 :Diverging connection (3 branches or more)
%NUM_PATHS=2 :Serial connection or two branches diverging
%NUM_PATHS=1 :Two latents only 
NUM_PATHS=length(L_tag);