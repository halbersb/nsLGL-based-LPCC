function [pdag,pdag_c,DAG,Observed,Latent,LL_score,BIC_score,clusters,clusters_size,major_clusters,clusters_unrounded,init_clusters,SM_250,labels,SOC] = LPCC(data,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main for LPCC for Desktop or SGE (cluster) - learn latent BN model
%
% input:
% [data]                 - can be string\matrix\structure
% [min_k,max_k]          - (scalar) inputs for automatic SOM model, upper and lower bounds for Davis Bouldin index to be evaluated
% [observed_cardinality] - (vector) the cardinality of each of the observed variables
% [C]                    - (matrix) clusters centers - if provided by the user
% [CS]                   - (vector) clusters sizes - if provided by the user
% [CUR]                  - (matrix) unrounded clusters centers - if provided by the user
% [pcc_type]             - [scalar] indicator for how to run pairwise cluster comparison (PCC)
% [flag]                 - [scalar] indicator to use EM (or just run LPCC stage 1)
% [flag2]                - [scalar] indicator to calculate BIC (flag must also be 1)
% [data2]                - same as data but for temporal data it includes all sequence rather than just 2-TBN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% if data is not a matrix
if ischar(data), data=load(data); end %if 'data' is a string with data path
if isstruct(data), data=data.data; end %if 'data' is a structure

%read input
for i=1:2:length(varargin)
    switch varargin{i},
        case 'max_k', max_k=varargin{i+1};
        case 'min_k', min_k=varargin{i+1};
        case 'observed_cardinality', observed_cardinality=varargin{i+1};
        case 'C', C=varargin{i+1};
        case 'CS', CS=varargin{i+1};
        case 'CUR', CUR=varargin{i+1};
        case 'server_ind', server_ind=varargin{i+1};
        case 'pcc_type', pcc_type=varargin{i+1};
        case 'flag', flag=varargin{i+1};
        case 'flag2', flag2=varargin{i+1};
        case 'data2', data2=varargin{i+1};
        otherwise, error(['invalid argument name: ' varargin{i}]);       
    end
end

%input validation
if nargin<5, error('not enough input arguments'); end
if exist('max_k','var') && ~isnumeric(max_k), error('max_k must be integer'); end
if exist('min_k','var') && ~isnumeric(min_k), error('min_k must be integer'); end
if exist('C','var') && size(C,2)~=size(data,2), error('number of variables in clusters is different from number of variables in data'); end
if exist('CS','var') && length(CS)~=size(C,1), error('the sizes of cluster vector is different from the number of clusters'); end
if exist('CUR','var') && (size(CUR,2)~=size(data,2) || size(CUR,1)~=size(C,1)), error('number of variables in unrounded clusters is different from number of variables in data or the number of cluster is different from C'); end
if exist('server_ind','var') && ~isnumeric(server_ind), error('server_ind must be 1 or 0'); end
if exist('pcc_type','var') && ~isnumeric(pcc_type), error('pcc_type must be 1, 2, or 3'); end
if exist('flag','var') && ~isnumeric(flag), error('flag must be 0 or 1'); end
if exist('flag2','var') && ~isnumeric(flag2), error('flag2 must be 0 or 1'); end
if ~exist('observed_cardinality','var'), observed_cardinality=max(data); end %assumeing each variable states are in the form 1,2,3...
if ~exist('server_ind','var'), server_ind=0; end
if ~exist('pcc_type','var'), pcc_type=1; end
if ~exist('flag','var'), flag=1; end
if ~exist('flag2','var'), flag2=1; end
if ~exist('data2','var'), data2=data; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find clusters with SOM if it was not provided as input %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clusters_unrounded={};init_clusters={};SM_250={};labels={};clusters_size={};
if exist('max_k','var') && exist('min_k','var') %if 'C' and 'CS' are also provided by user then ignore them
    min_k=max_k; %to average the centroids we need to assume that they are the same k size.
    %run automatic SOM
    for i=1:10
        [~,~,clusters_unrounded{i},init_clusters{i},SM_250{i},~]=runAutoSOM(data,max_k,min_k,server_ind,200);
        [~,clusters_size{i},clusters_unrounded{i},labels{i}]=z_reCalcCentroids(data,clusters_unrounded{i}); %Dan added: recalculate clusters not based on som unit (but directly from data)
    end
    try
        [clusters,clusters_size,clusters_unrounded]=z_avgClusters(clusters_size,clusters_unrounded); %Dan added: take 10 SOM results smooth it to a singe cluster martrix 
    catch
        try mkdir('clusters'); end
        save(['clusters/' 'temp' '.mat']); %save workspace for debuging
    end
else
    if ~exist('C','var') || ~exist('CS','var'), error('if max_k and min_k are not supply then C and CS must be given as input'); end
    %take from input
    clusters=C; clusters_size=CS;
    if exist('CUR','var'), clusters_unrounded=CUR; else clusters_unrounded=C; end %CUR is an optional input
end
[clusters,clusters_size,clusters_unrounded]=z_findDistinctClusters(clusters,clusters_size,clusters_unrounded,1); %Dan added: merge clusters if due to rounding they are equal

%%%%%%%%%%%%%%%%%%%
% initializations %
%%%%%%%%%%%%%%%%%%%

major_clusters{1}=find(clusters_size>mean(clusters_size)); %find initial major clusters
if length(major_clusters{1})<2, major_clusters{1}=sort([major_clusters{1} find(clusters_size==max(clusters_size([1:major_clusters{1}-1,major_clusters{1}+1:end])))]); end %if there is only one major cluster then add the next one as well (by its size)
MAX_ITER=10;
iter=1;
DAG={};
SOC={};

try
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find major clusters, latent exogenous and latent colliders %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while iter<MAX_ITER %iterations to find major cluster set
        [dag,Observed,Latent,SOC{iter}]=runLPCC(clusters,clusters_unrounded,pcc_type,data,labels,major_clusters{iter});
        if iter>1 && isequal(dag,DAG{iter-1}) %first stop condition: if the found dag is equal to the one found in the previous iteration
            break;
        end
        latents=[]; %array of latent indices
        observeds=[]; %array of observed indices
        for i=1:length(Latent)
            latents(i)=Latent(i).PCC_index;
        end
        for i=1:length(Observed)
            observeds(i)=Observed(i).PCC_index;
        end
        EXO_L=[]; %array of exogenous latent indices
        t=0;
        for i=1:length(Latent)
            if isempty(Latent(i).PA) %if the latent has no parents
                t=t+1;
                EXO_L(t)=Latent(i).PCC_index; %add the latent index to exogenous latents array
            end
        end
        [latents_cardinality,l_ch,ch_v]=findLatentCardinality(dag,latents,clusters,clusters_unrounded,major_clusters{iter},pcc_type); %find the cardinality of each latent based on the distinct configurations of its children. also, for each state in the cardinality find the corresponding children values - ch_v
        new_majors=sort(selectNewMajors(EXO_L,Observed,clusters,clusters_size,clusters_unrounded,l_ch,ch_v,[observed_cardinality,latents_cardinality],pcc_type)); %update the major cluster set

        if length(new_majors)<2 || isequal(major_clusters{iter},new_majors) %second stop condition: if a new major set was not found
            break;
        else %if a new major set was found
            major_clusters{iter+1}=new_majors;
        end
        DAG{iter}=dag;
        iter=iter+1;
        if iter==MAX_ITER, warning('iterations have reached the maximum value'); end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % end of stage 1 of LPCC % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %learn the bnet parameters (which are used to find first order minor clusters - 1-MC)
    cases_inc=createIncompleteData(data2,latents,observeds,ch_v,l_ch,pcc_type); %fill in as much as possible of the latents values - only possible to fill latent sample when the sample's latent children configuration is equal to one of the latent's ch_v
    for i=1:length(latents), latents_cardinality(i)=length(unique(cell2mat(cases_inc(latents(i),:)))); end %validate latents cardinality
    ns=[observed_cardinality,latents_cardinality];
    if flag
        [bnet,~,~,~]=learnParametersEM(dag,ns,cases_inc,length(latents));
        %to view the learned parameters, we use a little Matlab hackery 
        CPT=cell(1,length(ns));
        for i=1:length(ns)
            s=struct(bnet.CPD{i}); %violate object privacy
            CPT{i}=s.CPT;
        end
        %find max-minor and max-major for each latent
        [Max_minors,Max_majors,Max_priors]=findMaxPr(CPT,Observed,Latent,EXO_L,l_ch);
    else
        Max_minors=[];Max_majors=[];Max_priors=[];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find new non-collider latents - stage 2 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    max_PCC_index=Latent(end).PCC_index; %max_PCC_index is the index of the last latent (i.e., the last latent that was added)
    Latent_final2=[]; %init the final list of latents
    for l=1:length(Latent)
        Latent_final{l}=[]; %temp list of latents
        if length(Latent(l).OCH)>3 %at least 4 children (otherwise no split is possible) 
            [clust_new,~,clusters_size_new]=findNewClusters(clusters,clusters_size,clusters_unrounded,Observed,Latent,l_ch,ch_v,l,pcc_type); %marginalize the other latents, i.e., all the children of the other latents which have major values
            [minors,minors_s]=unUsedClusters(major_clusters{end},clust_new,clusters_size_new); %remove the major clusters from 'clust_new' since 1-MCs are smaller than the smallest major cluster 
            [try_big,~]=findBigClusters(minors,minors_s,l,size(data,1),Latent,Max_minors,Max_majors,Max_priors); %take only the big minor clusters (those that are bigger than 2-MCs)
            for m=1:length(major_clusters{end}) %Dan added: a loop over all major cluster
                %TODO consider loop over try_big each time with all major clustsers: for m=1:length(try_big) and: ...try_big(m),major_clusters{end},...
                [Observed,Latent_final{l},NUM_PATHS,~]=latentSplit2_pr3(clusters,clusters_unrounded,try_big,major_clusters{end}(m),l,Observed,Latent,l_ch{l},pcc_type);
                if length(Latent_final{l})>1 %if a split happened
                    break;
                end
            end
            if length(Latent_final{l})>1 %if a split happened
                %Dan added: NUM_PATHS==1
                if NUM_PATHS==1 %only two latents with an undirected edge between them
                    Latent_final{l}(1).PCC_index=Latent(l).PCC_index; %take the original (split) latent PCC_index
                    Latent_final{l}(1).name=['L',num2str(Latent_final{l}(1).PCC_index)];
                    max_PCC_index=max_PCC_index+1; %increase max_PCC_index by 1
                    Latent_final{l}(2).PCC_index=max_PCC_index; %new latent will be added to the end in terms of PCC_index
                    Latent_final{l}(2).name=['L',num2str(Latent_final{l}(2).PCC_index)];
                    %add an undirected edge between them
                    Latent_final{l}(1).UD=Latent_final{l}(2).PCC_index;
                    Latent_final{l}(2).UD=Latent_final{l}(1).PCC_index;
                end
                if NUM_PATHS==2 %serial connection or two branches diverging - undirected edges
                    for lf=1:length(Latent_final{l})
                        if (Latent_final{l}(lf).K{1})==0 %this is the EX of this serial connection and it replaces the original ex in the final graph
                            previous_PCC_index=Latent_final{l}(lf).PCC_index;
                            Latent_final{l}(lf).PCC_index=Latent(l).PCC_index; %take the original latent PCC_index
                            Latent_final{l}(lf).name=['L',num2str(Latent_final{l}(lf).PCC_index)];
                            for lf2=1:length(Latent_final{l})
                                if length(Latent_final{l}(lf2).PATH{1})>0
                                    if Latent_final{l}(lf2).PATH{1}.PCC_index==previous_PCC_index
                                        Latent_final{l}(lf2).PATH{1}=Latent_final{l}(lf);
                                    end
                                end
                            end
                        end
                    end
                    for lf=1:length(Latent_final{l})
                        if (Latent_final{l}(lf).K{1})>0
                            previous_PCC_index=Latent_final{l}(lf).PCC_index;
                            max_PCC_index=max_PCC_index+1;
                            Latent_final{l}(lf).PCC_index=max_PCC_index; %take the last added latent index
                            Latent_final{l}(lf).name=['L',num2str(Latent_final{l}(lf).PCC_index)];
                            for lf2=1:length(Latent_final{l})
                                if length(Latent_final{l}(lf2).PATH{1})>0
                                    if Latent_final{l}(lf2).PATH{1}.PCC_index==previous_PCC_index
                                        Latent_final{l}(lf2).PATH{1}=Latent_final{l}(lf); 
                                    end
                                end
                            end
                        end
                    end
                    %add the undirected edges
                    for lf=1:length(Latent_final{l})
                        for lf2=1:length(Latent_final{l})
                            if length(Latent_final{l}(lf2).PATH{1})>0
                                if Latent_final{l}(lf2).PATH{1}.PCC_index==Latent_final{l}(lf).PCC_index
                                    Latent_final{l}(lf2).UD=[Latent_final{l}(lf2).UD Latent_final{l}(lf).PCC_index];
                                    Latent_final{l}(lf).UD=[Latent_final{l}(lf).UD Latent_final{l}(lf2).PCC_index];
                                end
                            end
                        end
                    end
                end
                if NUM_PATHS>2 %diverging connection with 3 branches or more - with directed edges     
                    for x=1:NUM_PATHS
                        for lf=1:length(Latent_final{l})
                            if (Latent_final{l}(lf).K{x})==0 %this is the EX of this serial edge and replaces the original ex in the final graph
                                previous_PCC_index=Latent_final{l}(lf).PCC_index;
                                Latent_final{l}(lf).PCC_index=Latent(l).PCC_index; %take the original latent PCC_index
                                Latent_final{l}(lf).name=['L',num2str(Latent_final{l}(lf).PCC_index)];
                                for lf2=1:length(Latent_final{l})
                                    if length(Latent_final{l}(lf2).PATH{x})>0
                                        if Latent_final{l}(lf2).PATH{x}.PCC_index==previous_PCC_index
                                            Latent_final{l}(lf2).PATH{x}=Latent_final{l}(lf);
                                        end
                                    end
                                end
                            end
                        end
                        for lf=1:length(Latent_final{l})
                            if (Latent_final{l}(lf).K{x})>0
                                previous_PCC_index=Latent_final{l}(lf).PCC_index;
                                max_PCC_index=max_PCC_index+1;
                                Latent_final{l}(lf).PCC_index=max_PCC_index;
                                Latent_final{l}(lf).name=['L',num2str(Latent_final{l}(lf).PCC_index)];
                                for lf2=1:length(Latent_final{l})
                                    if length(Latent_final{l}(lf2).PATH{x})>0
                                        if Latent_final{l}(lf2).PATH{x}.PCC_index==previous_PCC_index
                                            Latent_final{l}(lf2).PATH{x}=Latent_final{l}(lf); 
                                        end
                                    end
                                end
                            end
                        end
                        %add the directed edges
                        for lf=1:length(Latent_final{l})
                            for lf2=1:length(Latent_final{l})
                                if length(Latent_final{l}(lf2).PATH{x})>0
                                    if Latent_final{l}(lf2).PATH{x}.PCC_index==Latent_final{l}(lf).PCC_index
                                        Latent_final{l}(lf2).D=[Latent_final{l}(lf2).D Latent_final{l}(lf).PCC_index];
                                    end
                                end
                            end
                        end 
                    end
                end
                %add all latents split from Latent(l) to the final latents list
                for lf=1:length(Latent_final{l})
                    Latent_final2=[Latent_final2 Latent_final{l}(lf)];
                end
            else
                %copy_Latent=Latent(l);
                %if ~isempty(Latent_final2), copy_Latent.PCC_index=Latent_final2(end).PCC_index+1; end %Dan added: update latent index
                Latent_final2=[Latent_final2 Latent(l)]; %if length(Latent_final{l})<1, then no split occurred although latent has more than 3 children
            end
        else
            %copy_Latent=Latent(l);
            %if ~isempty(Latent_final2), copy_Latent.PCC_index=Latent_final2(end).PCC_index+1; end %Dan added: update latent index
            Latent_final2=[Latent_final2 Latent(l)]; %if length(Latent(l).OCH)<=3, then there was not split. just add the curr latent as is
        end
    end
    pdag=createDag(Observed,Latent_final2); %add the new observed-latents edges after the split
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % post processing steps % 
    %%%%%%%%%%%%%%%%%%%%%%%%%

    %step 1: preserve v-structures of latents if there were any before the splits
    for l1=1:length(Latent) 
       for l2=1:length(Latent) 
           if dag(Latent(l1).PCC_index,Latent(l2).PCC_index)==1
               pdag(Latent(l1).PCC_index,Latent(l2).PCC_index)=1;
           end
       end
    end
    Latent=Latent_final2;
    %step 2: add bi-directed edges if necessary
    un_directed=[];
    for lf=1:length(Latent) 
        if ~isempty(Latent(lf).UD) 
            for ud=1:length(Latent(lf).UD)
                pdag(Latent(lf).PCC_index,Latent(lf).UD(ud))=2;
                pdag(Latent(lf).UD(ud),Latent(lf).PCC_index)=2;
                un_directed=[un_directed;sort([Latent(lf).PCC_index,Latent(lf).UD(ud)])];
            end
        end
    end
    un_directed=unique(un_directed,'rows'); %remove duplications
    if flag && flag2
        %step 3: calculate the new parameters and log-likelihood
        for i=1:length(Latent), latents(i)=Latent(i).PCC_index; end %update latents vector
        [latents_cardinality,l_ch,ch_v]=findLatentCardinality(pdag,latents,clusters,clusters_unrounded,major_clusters{end},pcc_type);
        cases_inc=createIncompleteData(data,latents,observeds,ch_v,l_ch,pcc_type);
        for i=1:length(latents), latents_cardinality(i)=length(unique(cell2mat(cases_inc(latents(i),:)))); end %validate latents cardinality
        ns=[observed_cardinality,latents_cardinality];
        %step 4: connect unconnected latents (based on PC algorithm)
        dag=z_covertToDag(pdag,un_directed);
        %run EM only if there was a split and EM flag is up
        if ~isequal(dag,bnet.dag)
            [~,cases_comp,LL_trace,BIC_score]=learnParametersEM(dag,ns,cases_inc,length(latents));
            LL_score=LL_trace(end);
            [pdag_c,~]=z_connectLatents(pdag,Latent,latents,cases_comp(sort(latents),:)); %PC algorithm. pdag_c is the pdag which was completed based on the PC
        else
            pdag_c=[];LL_score=[];BIC_score=[];
        end
    else
        pdag_c=[];LL_score=[];BIC_score=[];
    end
    %step 5: plot the graph including names of latents (relevant only for PC and not server)
    if server_ind
        names={};
        for i=1:size(data,2), names=[names,['X' int2str(i)]]; end
        for j=i+1:size(pdag,2), names=[names,['L' int2str(j-i)]]; end
        if sum(sum(pdag))>0
            bio=biograph(pdag,names);
            for k=i+1:j, bio.nodes(k).shape='ellipse'; end
            view(bio);
        end
    end
catch 
    if ~exist('pdag','var')
        pdag=[];LL_score=[];BIC_score=[];
    end
    if ~exist('pdag_c','var')
        pdag_c=[];LL_score=[];BIC_score=[];
    end
end