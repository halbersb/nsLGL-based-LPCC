function D = naneucdist(XI,XJ,weights)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function measure Euclidean distance while ignoring coordinates with NaNs
%
% input:
% [XI]      - (vector) the i'th sample
% [XJ]      - (matrix) the second sample (can be a matrix for multiple comparsions)
%
% optional:
% [weights] - weights for the features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    if nargin==2
        weights=ones(1,size(XI,2)); %all weights are 1
    end

    if size(XI,1)>1, error('XI is a single sample'); end
    n=size(XI,2);
    D=zeros(size(XJ,1),1);
    for i=1:size(XJ,1)
        Xj=XJ(i,:);
        sqdx=(XI-Xj).^2;
        sqdx=times(sqdx,weights);
        nstar=sum(~isnan(sqdx),2); %number of pairs that do not contain NaNs
        nstar(nstar==0)=NaN; %to return NaN if all pairs include NaNs
        Dsquared=nansum(sqdx,2).*n./nstar; %correction for missing coordinates
        D(i)=sqrt(Dsquared);
    end
end