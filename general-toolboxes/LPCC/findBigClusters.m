function [sms,max_th] = findBigClusters(c,c_s,l,ncases,Latent,Max_minors,Max_majors,Max_priors)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find subset of minors cluster (defenition is non-major) for which their
% size is greater than some threshold
% they are the FMC (first order minor clusters)
%
% input:
% [c]          - (vector) indices of minor clusters (FMC)
% [c_s]        - (vector) list of clusters size
% [l]          - (scalar) the curr latent index
% [ncases]     - (scalar) number of cases\samples in the data
% [Latent]     - (structure array) latents variables
% [Max_minors] - (array) array of vectors with max minors for each latent
% [Max_majors] - (array) array of vectors with max majors for each latent
% [Max_priors] - (vector) the maximum prior for each latent
%
% output:
% [sms]        - (vector) list of found clusters (under size condition)
%                         sms stands for sub minor set
% [max_th]     - (scalar) threshold based on p1-6 and ncases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
sms=[];
max_th=0;
lb=0.5;

if ~isempty(Max_minors)
    p0=prod(Max_priors);
    p1=prod(Max_minors{l}(1:2)); %exactly 2 observed children are minors
    p2=prod(Max_majors{l}(1:length(Latent(l).CH)-2)); %all the rest are majors
    p3=1;
    for k=1:length(Latent) %all children of all other latents are major
        if k~=l
            for r=1:length(Latent(k).OCH)
                p3=p3*Max_majors{k}(r);
            end
        end
    end
    %calculate the symetrical case (obtaining the same observed configutration)
    %but with a different value of the exogenous - two majors and the rest are minors
    p4=prod(Max_majors{l}(1:2)); %exactly 2 observed children are majors
    p5=prod(Max_minors{l}(1:length(Latent(l).CH)-2)); %all the rest are minors
    p6=p3;
    max_th=(p0*p1*p2*p3*ncases) + (p0*p4*p5*p6*ncases);
    t=0;
    for i=1:length(c_s)
        if c_s(i)>=max_th %take only minor clusters for which their size is greater than max_th (max threshold)
            t=t+1;
            sms(t)=c(i);
        end
    end
else
    sms=c;
    lb=0.6;
end

%% Dan added: in any case don't take lower 50% in terms of size of clusters
if length(c_s)>1 && length(sms)>1
    c_s=c_s(ismember(c,sms));
    temp=sort(c_s);
    sms(ismember(c_s,temp(1:ceil(length(temp)*lb))))=[];
end