function indivs = herd_migration(indivs,wij,kap,npatch,oldPatch)
if size(indivs,1) == 0
    error('population went extinct, herd_migration')
end

newPatches = zeros(size(indivs(:,4)));

for m = unique(indivs(:,4))' %1:npatch
    if m <= npatch
    thisPop = find(indivs(:,4)==m); % global id
    nNow = length(thisPop);
    if nNow > 0 %occupied patch
        
        wij_hat = sample_dirichlet(kap*wij(:,m),1);   % realized dispersal probabilities
        destin = randsample((npatch+1),nNow,true,wij_hat);
        newPatches(thisPop) = destin;
        
    end
    end
end
indivs(:,4) = newPatches;
indivs(:,3) = oldPatch;
