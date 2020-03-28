function [exv] = calcEXV(ns)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find all conbination of assigned values based on node sizes
% this is a recursive function
%
% input;
% [ns]  - (vector) nodes size
%
% output:
% [exv] - (matrix) all combinations
%                  columns are nodes rows are combinations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialization
exv=[];
t=0;

if isempty(ns), return; end %stop condition for the recursive function
if size(ns,2)>1
    exv_temp=calcEXV(ns(2:size(ns,2))); %recursive call
    for i=1:ns(1)
        for j=1:size(exv_temp,1)
            t=t+1;
            exv(t,:)=[i exv_temp(j,:)];
        end
    end
else
    for i=1:ns(1)
        t=t+1;
        exv(t,:)=i;
    end
end