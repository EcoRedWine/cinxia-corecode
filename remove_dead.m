function [indivs,families] = remove_dead(indivs,npatch,families)

if size(indivs,1) == 0
    error('population went extinct, remove dead')
end

dead = find(indivs(:,4) > npatch); % global id
if ~isempty(dead)
    indivs(dead,:) = [];
    families(dead,:) = [];
end 