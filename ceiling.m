function [indivs] = ceiling(indivs,npatch,capacities,metapopCV,cap_scaling,min_vbl)

% % Otso version
% k = log(capacities);
% %SD = metapopSD*k.^cap_scaling; %metapopSD*k(0.1857); %metapopSD*1.5./(metapopSD*k.^1.5);  %sqrt((metapopSD^2)/sum(k.^2)).*k; % 1.5,
% RR = randn(1);
% Xt = RR*metapopSD+1; %*SD+1;
% new_capacities = k.*Xt;%log(capacities.*Xt); %

% original version
k = log(capacities);
SD_i = metapopCV*k.^cap_scaling; %metapopSD*k(0.1857); %metapopSD*1.5./(metapopSD*k.^1.5);  %sqrt((metapopSD^2)/sum(k.^2)).*k; % 1.5,
RR = randn(1);
Xt = RR*SD_i+1; % SD = CV_patch_i
new_capacities = k.*Xt;%log(capacities.*Xt); %

% % yet another version -- does not really work (variability does not look
% right)
% %k = log(capacities);
% SD = metapopCV*capacities.^cap_scaling; %metapopSD*k(0.1857); %metapopSD*1.5./(metapopSD*k.^1.5);  %sqrt((metapopSD^2)/sum(k.^2)).*k; % 1.5,
% RR = randn(1);
% Xt = RR*SD+1;
% new_capacities = log(capacities.*Xt);%log(capacities.*Xt); %

% % O2 suggestion
% Xt = randn(1)*metapopCV+log(ratio);
% new_capacities = log(exp(Xt)*patches(:,1).^cap_scaling);

for m = unique(indivs(:,4))'%1:npatch
    if m <= npatch
        thisPop = find(indivs(:,4)==m); % global id
        
        sum_nKids = sum(indivs(thisPop,7));
        if (sum_nKids > exp(new_capacities(m)))
            survive = exp(new_capacities(m))/sum_nKids;
            
%         if rand(1) > survive % scramble competition
%             indivs(thisPop,7) = 0;
%         end
            
            % binomial survival
            ret = zeros(size(thisPop,1),1);
            
            for i = find(indivs(thisPop,7)>0)'%1:size(thisPop,1)
                ret(i) = sum(rand(1,indivs(thisPop(i),7))<survive);
            end
            indivs(thisPop,7) = ret;
        end
        
        if rand(1) < min_vbl/(sum(indivs(thisPop,7))) % small populations fail to spin nests
            indivs(thisPop,7) = 0;
        end
        
    end
end
