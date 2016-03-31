function [indivs] = mating(indivs,encounter,handling)

total_indivs = size(indivs,1);
if total_indivs == 0
    error('population went extinct, mating')
end

indivs = sortrows(indivs,[4,1]); % sort by natal patch and id

sires = zeros(1,size(indivs,1));
for m = unique(indivs(:,4))'
    thisPop = find(indivs(:,4) == m); % row number in indivs = global id
    
    if length(thisPop)==1 % alone
        sires(thisPop) = -9999; % no mating
        
    elseif length(thisPop) > 1
        thisPop_M = thisPop(find(indivs(thisPop,2)==0)); % find males, global id
        thisPop_F = thisPop(find(indivs(thisPop,2)==1)); % find females
        
        P_mate = 1-exp(-encounter*length(thisPop_M)/(1+handling*encounter*length(thisPop_F))); %%%%%%
        mating_F = thisPop_F(rand(size(thisPop_F)) < P_mate); % global id
        
        if (isempty(thisPop_M) || isempty(thisPop_F) || isempty(mating_F))% only one sex present or no mating females
            sires(thisPop) = -9999; % no mating
        elseif length(thisPop_M) == 1 && length(mating_F) <= 3 %%%%%%
            sires(mating_F) = thisPop_M;
        elseif length(thisPop_M) == 1 && length(mating_F) > 3 %%%%%%
            mating_F = randsample(mating_F,3,false);%%%%%%
            sires(mating_F) = thisPop_M;%%%%%%
        elseif length(thisPop_M) > 1 %> 1 %%%%%%
            attempts = randsample(thisPop_M,length(mating_F),true); % sample with replacement
           counts = histc(attempts,unique(thisPop_M)); %%%%%%
            for i = thisPop_M(counts > 3)'
                sel_tired_M = randsample(mating_F(attempts == i),counts(thisPop_M == i)-3);
                for j = 1:length(sel_tired_M)
                    attempts(mating_F == sel_tired_M(j)) = [];
                    mating_F(mating_F == sel_tired_M(j)) = [];
                end
            end %%%%%%
            sires(mating_F) = attempts;
        end
        
        sires(setdiff(thisPop,mating_F)) = -9999; % males or not mating females
    else
        error('error in mating')
    end
end

indivs(:,5) = sires;
end