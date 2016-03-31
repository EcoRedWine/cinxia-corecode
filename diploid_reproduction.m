function [indivs,flag,families,families_1older] = diploid_reproduction(indivs,nColInd,max_id,t,families,families_1older,n_L2,achiasma,n_chrom,SNP1,SNP2)

flag = 0;

if size(indivs,1) == 0
    disp('population went extinct, diploid reproduction')
    flag = 1;
    return
end

All_F = find(indivs(:,2)==1); % find females, global id

if isempty(All_F)
    flag = 1;
    size(All_F)
    disp('no females')
    return
end

% No of individuals
nin1 = size(indivs,1);

offSpring = zeros(nansum(indivs(:,7)),nColInd); %% can sometimes get an error due to fit = Inf or Nan CHECK
families_1older = families; % save the one generation back
families = zeros(nansum(indivs(:,7)),3); % new families vector is made here (families only keeps one parents-offspring generation)
tf = 0;

% reproduction
for i = 1:nin1
    nNo = indivs(i,7);
    
    if nNo > 0 && indivs(i,5) ~= -9999 % the latter condition is really needed if females are unmated
        tc = (tf+1):(tf+nNo);
        offSpring(tc,4) = indivs(i,4); % where their mothers reproduced
        offSpring(tc,3) = indivs(i,3); % where their mothers came from
        
        %%% only ecology version
        if n_L2 == 0 % no genetics
            
            %offSpring(tc,4) = indivs(i,4); % where their mothers reproduced
            families(tc,2:3) = repmat([indivs(i,1),indivs(indivs(i,5),1)],indivs(i,7),1); % Mom and Dad id
            
            %%% genetics version
        elseif n_L2 > 0
            
% s = RandStream('mt19937ar','Seed',seedi); %randi(2^20) or 1 (=seedi) for fixed seed
% RandStream.setGlobalStream(s)

            dam = indivs(i,:);
            sire =  indivs(indivs(i,5),:);
            
            for k = tc
                % SNPs
                %if n_L2 > 0
                % dams
                if achiasma == 0 %0 = yes recombination in females
                    sw = randn(1,n_L2)<0.5;
                elseif achiasma == 1 % 1 = no recombination in females
                    sw = (rand(1,n_chrom)<0.5); % 30 chromosomes are independent
                    
                    %%%%%sw = (ones(1,n_chrom)<0.5);
                    
                    
                    sw = repmat(sw,n_L2/n_chrom,1);
                    sw = reshape(sw,1,n_L2);
                end
                
                offSpring(k,SNP1(sw==1)) = dam(SNP1(sw==1));
                offSpring(k,SNP1(sw==0)) = dam(SNP2(sw==0));
                
                % sires
                sw = rand(1,n_L2);
                
                %%%%%sw = ones(1,n_L2);
                
                
                sw(1:(n_L2/n_chrom):n_L2) = sw(1:(n_L2/n_chrom):n_L2)<0.5; % first loci on each chrom
                for j = 2:(n_L2/n_chrom):n_L2
                    sw(j:j+n_L2/n_chrom-2) = (sw(j:j+n_L2/n_chrom-2)<0.0625); % second to the 9th loci
                end
                
                for j = 2:(n_L2/n_chrom):n_L2
                    for jj = j:j+n_L2/n_chrom-2
                        if sw(jj) == 0 % not switch
                            sw(jj) = sw(jj-1);
                        else %sw(j) == 1 % do switch
                            sw(jj) = abs(sw(jj-1)-1);
                        end
                    end
                end
                
                offSpring(k,SNP2(sw==1)) = sire(SNP1(sw==1));
                offSpring(k,SNP2(sw==0)) = sire(SNP2(sw==0));
    
            families(k,2:3) = [dam(1),sire(1)]; % Mom and Dad id
            end
        end
        %%%%%%%%%
    end
    tf = tf + nNo;
end

% 
% s = RandStream('mt19937ar','Seed',seedi); %randi(2^20) or 1 (=seedi) for fixed seed
% RandStream.setGlobalStream(s)
%rand(1)

% sex determination and assign id
offSpring(:,2) = rand(size(offSpring,1),1) > 0.5; % 0 = male, 1 = female
offSpring(:,1) = ((max_id+1):(max_id+size(offSpring,1)))';
families(:,1) = offSpring(:,1);

if size(offSpring,1) == 0
    disp(['No offspring, diploid_reproduction/mutation. time = ',num2str(t)])
    flag = 1;
    return
end

indivs = offSpring;
