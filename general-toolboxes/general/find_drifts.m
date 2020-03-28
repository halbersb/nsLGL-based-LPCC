function [drifts] = find_drifts(TBN,E,C)

%find temporal Cs
drifts=[];
for i=2:length(TBN)
    if ~isequal(TBN{i},TBN{i-1})
        drifts=[drifts i];
    end
end
%remove double drifts around Cs
to_remove=[];
for i=1:E:(C*E+1)
    if i>1 && ismember(i,drifts) && ismember(i-1,drifts)
        to_remove=[to_remove find(drifts==i-1)];
    end
end
drifts(to_remove)=[];