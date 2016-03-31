function indivs = fitnessD(indivs,npatch,r0,eta,cutoff)

if size(indivs,1) == 0
    error('population went extinct, fitnessD')
end

fit = zeros(size(indivs,1),1);

for m = unique(indivs(:,3))'%1:npatch 3rd column = old patch
    if m <= npatch
    thisPop = find(indivs(:,3)==m); % global id
    if ~isempty(thisPop)
        
rkThis = r0;
        
        w1 = randn(1);
        rates = exp(w1*eta + rkThis); % expected fitness
        %rates = exp(w1*eta + rkThis)*exp(-fit_dd_slope*size(thisPop,1))*(1-fit_dd_min)+fit_dd_min;
        %rates = exp(w1*eta + rkThis)*(exp(-fit_dd_slope*size(thisPop,1))*(1-fit_dd_min)+fit_dd_min)*(1-fit_pos_dd_min*exp(-fit_pos_dd_slope*(size(thisPop,1))));
        %rates = exp(w1*eta + rkThis)*size(thisPop,1)^-0.1;
        %rates = exp(rkThis)*gamrnd(eta,1/eta)*(1-0.9*exp(-0.007*(size(thisPop,1))));
        %rates = exp(w1*eta + rkThis)*exp(-0.005*size(thisPop,1))*(1-0.9*exp(-0.007*(size(thisPop,1))));
        
%         ETA = eta-0*sqrt(patches(m,1)); % eta decreases with patch area
%         if ETA < 0
%             ETA = 0;
%         end
%         rates = exp(w1*ETA + rkThis)*exp(-fit_dd_slope*size(thisPop,1))*(1-fit_dd_min)+fit_dd_min;

% rates = 0;
% for i = 1:ceil(patches(m,1)/0.5)
%     w1 = randn(1);
%     rates = rates + exp(w1*eta + rkThis);
% end
% rates = (rates/i)*exp(-fit_dd_slope*size(thisPop,1))*(1-fit_dd_min)+fit_dd_min;
        
        while sum(rates>cutoff)>0
            w1 = randn(1); 
            
            rates=exp(w1*eta + rkThis); 
            %rates=exp(w1*eta + rkThis)*exp(-fit_dd_slope*size(thisPop,1))*(1-fit_dd_min)+fit_dd_min;
            %rates = exp(w1*eta + rkThis)*(exp(-fit_dd_slope*size(thisPop,1))*(1-fit_dd_min)+fit_dd_min)*(1-fit_pos_dd_min*exp(-fit_pos_dd_slope*(size(thisPop,1))));
            %rates = exp(w1*eta + rkThis)*size(thisPop,1)^-0.1;
            %rates = exp(rkThis)*gamrnd(eta,1/eta)*(1-0.9*exp(-0.007*(size(thisPop,1))));
            %rates = exp(w1*eta + rkThis)*exp(-0.005*size(thisPop,1))*(1-0.9*exp(-0.007*(size(thisPop,1))));
            
            %rates = exp(w1*ETA + rkThis)*exp(-fit_dd_slope*size(thisPop,1))*(1-fit_dd_min)+fit_dd_min;
            
% rates = 0;
% for i = 1:ceil(patches(m,1)/0.5)
%     w1 = randn(1);
%     rates = rates + exp(w1*eta + rkThis);
% end
% rates = (rates/i)*exp(-fit_dd_slope*size(thisPop,1))*(1-fit_dd_min)+fit_dd_min;

        end
                
        fit(thisPop) = rates; %*size(thisPop,1)^-0.25;
    end
    end
end
livingdead = (indivs(:,4) > npatch); % these have died during migration, but may be someone's sires; genotypes are needed
fit(livingdead) = 0; %figure(1);histogram(fit(indivs(:,2)==1),41)
fit(indivs(:,5) == -9999) = 0; % males or not mating females
indivs(:,6) = fit; 