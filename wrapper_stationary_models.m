function wrapper_stationary_models(G,C,E,S,alg,flag,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run wrapper
%
% input:
% G             - [scalar] indicator which graph to run
% C             - [scalar] number of changes
% E             - [scalar] number of slices between changes - epoch
% S             - [scalar] the required data size to be sampled
% alg           - [scalar] indicator which algorithm to run: 
%                          1-SEN_DBN, 2-stationary LPCC
% flag          - [scalar] 1 - use prelearn local graphs
%
% optional:
% ns            - [scalar] node size for all nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read input
for i=1:2:length(varargin)
    switch varargin{i},
        case 'ns', ns=varargin{i+1};
        otherwise, error(['invalid argument name: ' varargin{i}]);       
    end
end
if ~exist('ns','var'), ns=2; end
data_perm=5;
display(['G=' int2str(G) ' , C=' int2str(C) ' , E=' int2str(E) ' , S=' int2str(S) ' , Algo=' int2str(alg) ' , ns=' int2str(ns) ' , D_M=' int2str(1)])


sign='\';
save_path='after_databases\';

switch G
    case 1
        file_name=['G1_C' int2str(C) '_E' int2str(E) '_O3_L1_Size' int2str(S) '_P9_N' int2str(ns) '_DM' int2str(1)];
        num_slices=E*(C+1);
        num_var_slice=4;
        num_obs_slice=3;
        if ns==2, k=6; end %k for k-means
        if ns==3, k=8; end %k for k-means
        node_sizes=repmat(ns,1,num_slices*num_var_slice);
    case 2
        file_name=['G2_C' int2str(C) '_E' int2str(E) '_O6_L2_Size' int2str(S) '_P9_N' int2str(ns) '_DM' int2str(D_M)];
        num_slices=E*(C+1);
        num_var_slice=8;
        num_obs_slice=6;
        if ns==2, k=10; end %k for k-means
        if ns==3, k=16; end %k for k-means
        node_sizes=repmat(ns,1,num_slices*num_var_slice);
    otherwise, error('invalid experiment');
end

for p=1:data_perm
    
    load(['databases' sign 'G' int2str(G) sign file_name sign file_name '_' int2str(p) '.mat']); 
    data=data(1:0.9*size(data,1),:); %leave the last 10% for data imputation (testing/evaluation)
    TBN=[];final_2TBN=[];init_2TBN=[];pdags=[];
    
    switch alg,
        case 1
            %run SEM-DBN
            if p==1, save_path=[save_path 'After_SEM_DBN' sign]; end
            if ~flag
                [TBN,init_2TBN,~,~]=learn_struct_nsdbn_EM(mk_random_2TBN(num_var_slice,num_obs_slice,1,1,G*E*S*p),data,node_sizes,num_var_slice-num_obs_slice,10,1);
            else
                load([save_path file_name sign file_name '_' int2str(p) '.mat']);
            end
        case 2
            %run LPCC per slice
            if p==1, save_path=[save_path 'After_LPCC' sign]; end
            if ~flag
                for t=1:num_slices-1
                    try
                        [pdags{t},~,~,~,~,~,~,~,~,~,~,~,~,~,~]=LPCC(data(:,(((t-1)*num_obs_slice+1):(t+1)*num_obs_slice)),'min_k',k,'max_k',k,'pcc_type',1,'flag',1,'flag2',0);
                    catch
                        if ~exist([save_path 'errorsLPCC' sign file_name],'dir'), mkdir([save_path 'errorsLPCC' sign file_name]); end
                        save([save_path 'errorsLPCC' sign file_name sign file_name '.mat']); %save workspace for debuging
                    end
                end
            else
                load([save_path file_name sign file_name '_' int2str(p) '.mat']); %load pdags (local graphs)
            end
            try final_2TBN=combine_local_struct(pdags,num_obs_slice*2); end 
        otherwise, error('unrecognize algorithm input');
    end
    
    if ~exist([save_path file_name],'dir'), mkdir([save_path file_name]); end
    switch alg,
        case 1
            save([save_path file_name sign file_name '_' int2str(p) '.mat'],'TBN','init_2TBN','pdags');
        case 2
            save([save_path file_name sign file_name '_' int2str(p) '.mat'],'final_2TBN','init_2TBN','pdags');
    end

end