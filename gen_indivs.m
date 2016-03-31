function indivs = gen_indivs(Omega,npatch,nColInd,n_L2,allele_freq)

total_indivs = sum(Omega(:,1));

if total_indivs > 0
    
    indivs = zeros(total_indivs,nColInd);
    indivs(:,1) = 1:total_indivs; % individual id
    
    for m = 1:npatch
        indivs(sum(Omega(1:(m-1),1))+1:sum(Omega(1:m,1)),4) = m; % which patch
    end
    
    if n_L2 > 0
        
        uniq_patches = unique(indivs(:,4))';
        for i = uniq_patches
            thisPop = find(indivs(:,4) == i);
            temp = rand(length(thisPop),n_L2);
            
            
            %%%%%temp = ones(length(thisPop),n_L2);
            
            
            A = (temp < repmat(allele_freq,length(thisPop),1).^2); % 11
            B = (temp >= repmat(allele_freq,length(thisPop),1).^2 & temp < repmat(allele_freq,length(thisPop),1).^2+repmat(allele_freq,length(thisPop),1).*(1-repmat(allele_freq,length(thisPop),1))); % 10
            C = (temp >= repmat(allele_freq,length(thisPop),1).^2+repmat(allele_freq,length(thisPop),1).*(1-repmat(allele_freq,length(thisPop),1)) & temp < repmat(allele_freq,length(thisPop),1).^2+2*repmat(allele_freq,length(thisPop),1).*(1-repmat(allele_freq,length(thisPop),1))); % 01
            
            indivs(thisPop,8:7+2*n_L2)=[A,A]+[B,zeros(size(B))]+[zeros(size(C)),C];
            
            %         test = indivs(thisPop,5+2*n_L1:4+2*n_L1+n_L2)+indivs(thisPop,5+2*n_L1+n_L2:4+2*n_L1+2*n_L2);
            %         plot(i,mean(sum(test == 2)/length(thisPop)),'.'); hold on
            %         plot(i,mean(sum(test == 1)/length(thisPop)),'ro')
            %
        end
    end
    
% s = RandStream('mt19937ar','Seed',seedi); %randi(2^20) or 1 (=seedi) for fixed seed
% RandStream.setGlobalStream(s)

    % sex determination (second column)
    indivs(:,2) = (rand(total_indivs,1) > 0.5); % 0 = male, 1 = female
    
else
    error('Error, gen_indivs.m')
end
