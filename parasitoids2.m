function [indivs] = parasitoids2(indivs,npatch,fit_dd_slope,fit_dd_min,fit_pos_dd_min,fit_pos_dd_slope)

% remove parasitized individuals

if size(indivs,1) == 0
    error('population went extinct, parasitoids')
end

for m = 1:npatch % negative dd
    thisPop = find(indivs(:,3)==m); % function of the natal patch
    if ~isempty(thisPop)
        %prob_surv_paras = (exp(-fit_dd_slope*size(thisPop,1))*(1-fit_dd_min)+fit_dd_min);
        
        
        if size(thisPop,1) <= fit_dd_min
        prob_surv_paras = (-fit_dd_slope/fit_dd_min)*size(thisPop,1)+1;
        elseif size(thisPop,1) > fit_dd_min
            prob_surv_paras = 1-fit_dd_slope;
        end
        
        who_survive = rand(size(thisPop,1),1)<prob_surv_paras;
        indivs(thisPop(who_survive==0),6) = 0;
    end
end 

for m = 1:npatch % positive dd % 1-A*exp(-B*N), monotonically increasing curve
    thisPop = find(indivs(:,4)==m); % function of the destination patch
    if ~isempty(thisPop)
        prob_surv_paras = (1-fit_pos_dd_min*exp(-fit_pos_dd_slope*(size(thisPop,1))));
        
        who_survive = rand(size(thisPop,1),1)<prob_surv_paras;
        indivs(thisPop(who_survive==0),6) = 0;
    end
end