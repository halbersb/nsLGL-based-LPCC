function wrapper_LPCC_based_nsLGL(G,C,E,S,alg,window,jump,flag,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run wrapper
%
% input:
% G             - [scalar] indicator which graph to run
% C             - [scalar] number of changes
% E             - [scalar] number of slices between changes  - epoch
% S             - [scalar] number of samples
% alg           - [scalar] indicator which algorithm to run: 
%                          1-sliding window, 2-weighted version 
% window        - [scalar] number of slices for LGL-LPCC
% jump          - [scalar] number of slices to be jumped between windows
% flag          - [scalar] 1 - use prelearn local graphs
%
% optional:
% ns            - [scalar] node size for all nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read input
for i=1:2:length(varargin)
    switch varargin{i},
        case 'ns', ns=varargin{i+1};
        otherwise, error(['invalid argument name: ' varargin{i}]);       
    end
end
if ~exist('ns','var'), ns=2; end
data_perm=5;
display(['G=' int2str(G) ' , C=' int2str(C) ' , E=' int2str(E) ' , S=' int2str(S) ' , Algo=' int2str(alg) ' , W=' int2str(window) ' , ns=' int2str(ns)])

%add folder to path
sign='\';
save_path='after_databases\';

switch G
    case 1
        file_name=['G1_C' int2str(C) '_E' int2str(E) '_O3_L1_Size' int2str(S) '_P9_N' int2str(ns) '_DM' int2str(1)];
        num_slices=E*(C+1);
        num_obs_slice=3;
        if ns==2, k=6; end %k for k-means
        if ns==3, k=8; end %k for k-means
    case 2
        file_name=['G2_C' int2str(C) '_E' int2str(E) '_O6_L2_Size' int2str(S) '_P9_N' int2str(ns) '_DM' int2str(1)];
        num_slices=E*(C+1);
        num_obs_slice=6;
        if ns==2, k=10; end %k for k-means
        if ns==3, k=16; end %k for k-means
    otherwise, error('invalid experiment');
end

for p=1:data_perm

    load(['databases' sign 'G' int2str(G) sign file_name sign file_name '_' int2str(p) '.mat']); 
    data=data(1:0.9*size(data,1),:); %leave the last 10% for data imputation (testing/evaluation)

    switch alg,
        case 1
            %sliding window
            save_name=[file_name '_W' int2str(window) '_J' int2str(jump)];
            if p==1, save_path=[save_path 'After_Sliding' sign]; end
            if ~flag
                [TBN,pdags]=sliding_window(data,num_slices,num_obs_slice,2,1,k,save_path,save_name,sign,1);
            else
                load([save_path save_name sign save_name '_' int2str(p) '.mat']) %load pdags (local graphs)
                [TBN,pdags]=sliding_window(data,num_slices,num_obs_slice,2,1,k,save_path,save_name,sign,1,'pdagsPr',pdags);
            end
        case 2
            %weighted algorithm (jump always=1) no window (i.e., window=2)
            save_name=file_name;
            if p==1, save_path=[save_path 'After_Weights' sign]; end
            if ~flag
                [~,pdags]=sliding_window(data,num_slices,num_obs_slice,2,1,k,save_path,save_name,sign,0);
            else
                load([save_path save_name sign save_name '_' int2str(p) '.mat']) %load pdags (local graphs)
                [~,pdags]=sliding_window(data,num_slices,num_obs_slice,2,1,k,save_path,save_name,sign,0,'pdagsPr',pdags);
            end
            try for t=1:length(pdags), pdags{t}=pdags{t}{1}; end; end
            for t=1:length(pdags)
                TBN{t}=combine_local_struct_forWeights(pdags,num_obs_slice*2,t); %here t is t_star
            end
        otherwise, error('unrecognize algorithm input');
    end

    if ~exist([save_path save_name],'dir'), mkdir([save_path save_name]); end
    save([save_path save_name sign save_name '_' int2str(p) '.mat'],'TBN','pdags');
    
end