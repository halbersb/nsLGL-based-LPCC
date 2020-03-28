function [score] = measure_drifts(vec1,vec2,window,flag)

%count hits: how many out of vec1 (true) exist in vect (pred)
%flag: an indicator for precision or recall
if window<4,window=1; else window=2; end
score=0;
if ~isempty(vec2)
    for i=1:length(vec1)
        if ismember(vec1(i),vec2)
            score=score+1;
        else %try to increase according to window
            for j=1:window
                if ismember(vec1(i)-j,vec2) || ismember(vec1(i)+j,vec2)
                    score=score+1;
                    break;
                end
            end
        end
    end
end

if flag
    %recall
    score=score/length(vec1); %max score == length(vec1)
else
    %precision
    if ~isempty(vec2)
        score=score/length(vec2); %max score == length(vec2)
    else
        score=0;
    end
end