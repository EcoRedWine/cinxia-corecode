function [indivs,flag] = inbreeding(indivs,k,inbreed,t,families)

flag = 0;

if size(indivs,1) == 0
    disp('population went extinct, diploid reproduction')
    flag = 1;
    return
end

% reproduction
All_F = find(indivs(:,2)==1); % find females, global id

if isempty(All_F)
    flag = 1;
    size(All_F)
    disp('no females')
    %nKids = [];
    return
end

% correct fitness based on relatedness of pair
if inbreed == 1 % 1 = inbreeding ON, 0 = inbreeding OFF
    if t > 1
        for i = 1:length(All_F)
            if indivs(All_F(i),5) ~= -9999 % mating females (males and not mating females have -9999)
                parents = [families(All_F(i),2:3),families(indivs(All_F(i),5),2:3)];
                
                if length(parents) < 4
                    error('there are not enough parents')
                end
                
                if length(unique(parents)) == 2
                    indivs(All_F(i),6) = indivs(All_F(i),6)*k; % 30% reduction in fitness in the first round of sib-sib mating
                end
                
            end
        end
    end
end

indivs(All_F,7) = poissrnd(indivs(All_F,6)); % nKids realized fecundity

%rand(1)

indivs(indivs(:,5)==-9999,7) = 0;  % males have no fitness and unsired females have no fitness