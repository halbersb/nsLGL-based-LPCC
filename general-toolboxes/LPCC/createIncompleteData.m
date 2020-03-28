function [cases_inc] = createIncompleteData(data,latents,observeds,ch_v,l_ch,pcc_type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extend data to include latent values - incomplete data
%
% create the incomplete data by extending the data to include some of the 
% values of the latent variables, filled in by the corresponding value 
% of their children in each case
%
% input:
% [data]                 - (matrix) data to be used [rows=samples,columns=vars]
% [latents]              - (vector) list of latents indices
% [observeds]            - (vector) list of observed indices
% [observed_cardinality] - (vector) cardinality of each observed variable
% [cardinality]          - (vector) list with the cardinality of each latent variable
% [ch_v]                 - (array)  array of matrices with the latent children value configurations
% [l_ch]                 - (array)  array of vectors with the latent children
% [pcc_type]             - [scalar] indicator for how to run pairwise cluster comparison
%
% output:
% [cases_inc]            - (cell matrix) incomplete data [rows=vars,columns=samples]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
ncases_inc=size(data,1); %number of row in the database
cases_inc=cell(length(observeds)+length(latents),ncases_inc); %allocate cell matrix
cases_inc(1:length(observeds),:)=num2cell(data'); %fill observed data

for f=1:length(latents)
    for k=1:ncases_inc
        cases_inc_ch=zeros(1,length(l_ch{f}));
        for ch=1:length(l_ch{f})
            cases_inc_ch(ch)=cell2num(cases_inc(l_ch{f}(ch),k)); %vector of children of the current latent
        end
        %loop over all configurations\joint values of centers of the current latents children (ch_v) and search for a match to 'cases_inc_ch'
        %for synthetic data all observed in the vector will probably have the same value since we set max(p) where i=j in the matrix (i in current observed and j in latent variable)
        for i=1:size(ch_v{f},1)
            %if the relevant observed variables in the current data sample are equal to the i unique configuration (in ch_v) then fill the sample with i in 'f' colums
            %otherwise the sample will remains empty in the 'f' column - this is because the subset vector of observed is due to minor effects
            if isequal(createPcc([cases_inc_ch;ch_v{f}(i,:)],pcc_type),zeros(1,length(l_ch{f})))
                cases_inc{latents(f),k}=i; %i represent the i state of the latent (out of cardinality of size(ch_v{f},1))
                break;
            end
        end
    end
end