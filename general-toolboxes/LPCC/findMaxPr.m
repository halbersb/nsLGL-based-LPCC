function [Max_minors,Max_majors,Max_priors] = findMaxPr(CPT,Observed,Latent,EXO_L,l_ch) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the max probabilty for minor per Obs for each latent
% and the max probabilty for major per Obs for each latent
% and the max prioir for each latent
%
% input:
% [CPT]        - (cell) conditional probabilities tables
% [Observed]   - (structure array) observed variables
% [Latent]     - (structure array) latents variables
% [EXO_L]      - (vector) indices of exogenous latent variables
% [l_ch]       - (array) array of vectors each contains children of latent
%
% output:
% [Max_minors] - (array) array of vectors with max minors for each latent
% [Max_majors] - (array) array of vectors with max majors for each latent
% [Max_priors] - (vector) the maximum prior for each latent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Max_priors=[];
Max_minors={};
Max_majors={};

%find the maximum prior value for each latent
for i=1:length(EXO_L)
    max_prob=0;
    for k=1:size(CPT{EXO_L(i)},1) %pass over the row of the exo latent in the CPT
        if CPT{EXO_L(i)}(k)>max_prob
            max_prob=CPT{EXO_L(i)}(k);
        end
    end
    Max_priors(i)=max_prob;
end

%create max minor and major vectors per latent
for i=1:length(Latent)
    Max_minors_ch=[];
    Max_majors_ch=[];
    for ch=1:length(l_ch{i})
        if l_ch{i}(ch)<=length(Observed) %if it is an observed child
            min_p=[];
            max_p=[];
            for k=1:size(CPT{l_ch{i}(ch)},1) %pass over the row of the child in the CPT
                min_p=[min_p min(CPT{l_ch{i}(ch)}(k,:))]; %add min of the row
                max_p=[max_p max(CPT{l_ch{i}(ch)}(k,:))]; %add max ofthe row
            end
            if ~isempty(min_p), max_minor_ch=max(min_p); else max_minor_ch=0; end % take te max out of the minor (not min!)
            if ~isempty(max_p), max_major_ch=max(max_p); else max_major_ch=0; end
        end
        Max_minors_ch=[Max_minors_ch max_minor_ch];
        Max_majors_ch=[Max_majors_ch max_major_ch];    
    end
    Max_minors{i}=sort(Max_minors_ch,'descend'); %sort the probabiliteis for the curr latent
    Max_majors{i}=sort(Max_majors_ch,'descend');
end