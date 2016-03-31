function [opThis] = add_specs(opThis,indivs,studyPatches,t,npatch,SNP1,SNP2,n_L2,seedi)

if size(indivs,1) == 0
    error('population went extinct, add_specs')
end

if ~isnan(SNP1)
    newOP = zeros(length(studyPatches),6);
    gk = indivs(:,SNP1)+indivs(:,SNP2); % 2 = homozyg recessive
    elim = ((sum(gk)/(2*size(indivs,1)) < 0.2) + (sum(gk)/(2*size(indivs,1)) > 0.8) == 1); % eliminate minor alleles
else
    newOP = zeros(length(studyPatches),4);
end

j = 0;

for i = studyPatches
    j = j + 1;
    newOP(j,1) = t;
    newOP(j,2) = i;
    thisPop = find(indivs(:,4)==i);
    n1 = length(thisPop);
    newOP(j,3) = n1;
    
    if n1 > 0 % age
        if (t == 1 || isnan(opThis((t-2)*npatch+j,4)))
            newOP(j,4) = 1;
        else
            newOP(j,4) = opThis((t-2)*npatch+j,4)+1;
        end
    else
        newOP(j,4) = NaN;
    end
    
    % summarize genetics
    if ~isnan(SNP1)
        if ~isempty(thisPop)
            gk = indivs(thisPop,SNP1)+indivs(thisPop,SNP2); % 2 = homozyg recessive
           
            heterozygotes = sum(sum(gk == 1)); % number of heterozygous LOCI in the pooled loci over the population
            newOP(j,5) = heterozygotes/(size(thisPop,1)*(n_L2));
            
            gk(:,elim)=[]; % eliminate minor alleles
            heterozygotes = sum(sum(gk == 1)); % number of heterozygous LOCI in the pooled loci over the population, excluding loci with minor allele frequencies < 0.2
            newOP(j,6) = heterozygotes/(size(thisPop,1)*(n_L2-sum(elim)));
            
        else
            newOP(j,5) = NaN;
            newOP(j,6) = NaN;
        end
    end
    
% s = RandStream('mt19937ar','Seed',seedi); %randi(2^20) or 1 (=seedi) for fixed seed
% RandStream.setGlobalStream(s)


end

opThis = [opThis;newOP];