function [bnet,cases_comp,LLtrace,bic] = learnParametersEM(dag,ns,cases_inc,n_lat,flag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the EM algorithm to learn the parameters of the bnet using the
% incomplete data created
%
% input:
% [dag]        - (matrix) the found dag
% [ns]         - (vector) nodes dize
% [cases_inc]  - (cell matrix) rows are variables and columns are samples
% [n_lat]      - (scalar) number of latent variables
%
% optional:
% [flag]       - [scalar] mode. if flag exists try to fill all variables
%
% output:
% [bnet]       - (object) a Bayesian network structure
% [cases_comp] - (cell matrix) the filled-in matrix [rows=vars,columns=samples]
% [LLtrace]    - (vector) all loglikelihood srcores according to max_iter
% [bic[        - (scalar) the final BIC score
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5, flag=0; end

%% create bnt
bnet=mk_bnet(dag,ns);
for i=1:length(ns)
    bnet.CPD{i}=tabular_CPD(bnet,i);
end
engine=jtree_inf_engine(bnet);
max_iter=10;
[bnet,LLtrace,engine2]=learn_params_em(engine,cases_inc,max_iter);
cases_comp=cases_inc;

%fill in empty cells in the data according to em parameters
n_obs=length(ns)-n_lat;
for i=1:size(cases_comp,2)
    [engine3,~]=enter_evidence(engine2,cases_comp(:,i));
    %loop over all latents in the sample
    if flag, n_lat=n_obs+n_lat;n_obs=0; end %this row was added on 09-18 for AISTAT exp. to fill-in only latent variables
    for j=1:n_lat
        %if it's empty then fill it in
        if isempty(cell2mat(cases_comp(n_obs+j,i)))
            m=marginal_nodes(engine3,n_obs+j);
            cases_comp{n_obs+j,i}=find(m.T==max(m.T),1);
        end
    end
end
%calculate bic score (cases_inc are now complete)
d{1}=dag;
try
    bic=score_dags(cases_comp,ns,d,'discrete',1:length(dag),'scoring_fn','bic');
catch
    bic=[];
end