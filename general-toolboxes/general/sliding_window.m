function [TBN,pdags] = sliding_window(data,num_slices,num_obs_slice,window,jump,k,save_path,save_name,sign,flag,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run sliding_window algorithm
%
% input:
% data          - [matrix] data set
% num_slices    - [scalar] number of time-slices
% num_obs_slice - [scalar] number of observed variable per slice
% window        - [scalar] window size for the local learning
% jump          - [scalar] jumps between slices per window size 
% k             - [scalar] k for k-means of LPCC
% save_path     - [scalar] path to save the results
% save_name     - [scalar] mane of the saved result file
% sign          - [char] '\' or '/'
% flag          - [scalar] indicator if the function was called from
%                          sliding or weights algorithm
%
% optional:
% pdagsPr       - [cell array] local graphs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read input
for i=1:2:length(varargin)
    switch varargin{i},
        case 'pdagsPr', pdagsPr=varargin{i+1};
        otherwise, error(['invalid argument name: ' varargin{i}]);       
    end
end

pdags=cell(1);TBN=cell(1);
i=1;
for ind=1:jump:(num_slices-window+1)
    if ~exist('pdagsPr','var')
        for t=ind:(window+ind-2) %always 1 when called from 'AlgoWeights' (with window=2)
            try
                curr_data=data(:,(((t-1)*num_obs_slice+1):(t+1)*num_obs_slice));
                curr_data(any(isnan(curr_data),2),:)=[]; %remove rows with NaN
                if size(curr_data,1)<200, continue; end %break the loop if data is too small
                [pdags{i}{t-ind+1},~,~,~,~,~,~,~,~,~,~,~,~,~,~]=LPCC(curr_data,'min_k',k,'max_k',k,'pcc_type',1,'flag',1,'flag2',0);
                clear curr_data
            catch
                if ~exist([save_path 'errorsLPCC' sign save_name],'dir'), mkdir([save_path 'errorsLPCC' sign save_name]); end
                save([save_path 'errorsLPCC' sign save_name sign save_name '.mat']); %save workspace for debuging
            end
        end
    else
        pdags=pdagsPr;
    end
    if flag %a condition to avoid 'AlgoWeights'
        try
            TBN{ind}=combine_local_struct(pdags{i},num_obs_slice*2);
        catch
            if ~exist([save_path 'errorsLG' sign save_name],'dir'), mkdir([save_path 'errorsLG' sign save_name]); end
            save([save_path 'errorsLG' sign save_name sign save_name '.mat']); %save workspace for debuging
        end
    end
    i=i+1;
end