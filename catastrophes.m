function [indivs,flag,families] = catastrophes(indivs,catast,families)

flag = 0;

A = unique(indivs(:,4)); 
B = (rand(size(A,1),1)< catast); % 1 = to be wiped out
catast_patch = A(B);

% AA = hist(indivs(:,4),A)';
% BB = logical((AA > 5 & AA < 15).*(rand(size(AA))<0.3));
% catast_patch = A(BB);

% A = unique(indivs(:,4)); 
% A = intersect(A,find(patches(A,1)<0.5));
% B = (rand(size(A,1),1)< catast);
% catast_patch = A(B);

for i = 1:length(catast_patch)
   families(indivs(:,4) == catast_patch(i),:) = [];
   indivs(indivs(:,4) == catast_patch(i),:) = [];   
end

if size(indivs,1) == 0
    disp('population went extinct, catastrophe')
    flag = 1;
    return
end
