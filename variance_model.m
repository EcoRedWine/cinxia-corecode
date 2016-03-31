function [opThis,indivs,migrants,nKids_out] = variance_model(outPutData,nMC,T,iter,seedi,tSkip,Omega,patches,params,allele_freq,repid,save_migs,save_indivs,save_indivs_start,census_time,transSIN_mat,indivs)
%% set up
s = RandStream('mt19937ar','Seed',seedi); %randi(2^20) or 1 (=seedi) for fixed seed
RandStream.setGlobalStream(s)

burnin = T-iter;
%nColInd = 7;

%% parameters
[r0,gamma,environSD,locSD,eps,migDeath,alfa,kap,n_L1,n_L2,ratio,migr_opt,MR1,MR2,K_area_scaling,eta,zeta_em,...
    encounter,handling,metapopCV,achiasma,n_chrom,k,cutoff,inbreed,zeta_im,catast,fit_dd_slope,fit_dd_min,fit_pos_dd_min,fit_pos_dd_slope,cap_scaling,min_vbl] = fix_params(params);

% set up the landscape
npatch = size(patches,1); % individual patches already given in input
studyPatches = 1:npatch;

wij = calculate_distanceD(patches,npatch,alfa,migDeath,eps,zeta_em,zeta_im,transSIN_mat); % immigration probabilities
capacities = ratio*(patches(:,1)).^K_area_scaling;  % carrying capacities
patches = [patches;zeros(1,4)]; % adding a dummy patch for keeping track of dead sires

% indices for genetic materials (SNP only)
nColInd = 7+n_L2*2;
if n_L2 > 0 % if inbreeding is modeled
    SNP1 = 8:7+n_L2;
    SNP2 = 8+n_L2:nColInd;
else
    SNP1 = NaN;
    SNP2 = NaN;
end

% placeholders
opThis = [];
families = []; % need this for t = 0 (after that disploid_reproduction will generate it)
families_1older = [];
migrants = [];
%counts = zeros(100,1);
nKids_out = [];

% check input consistency
if size(Omega,1) ~= npatch
    error('dim of Omega and npatch not match, gen_indivs')
end
if sum(Omega(:,1)) == 0
    error('Population does not exist, initialization of indivs')
end

% initialize individuals if indiv == []
if isempty(indivs)
    indivs = gen_indivs(Omega,npatch,nColInd,n_L2,allele_freq);
end

for t = 1:T
    
    %rand(1)
    
    max_id = max(indivs(:,1));
    
    [indivs] = mating(indivs,encounter,handling);
    
    % Migration
    oldPatch = indivs(:,4);
    [indivs] = herd_migration(indivs,wij,kap,npatch,oldPatch);
    
%     if (save_migs == 1 && t >= save_indivs_start)
%        migrants = [migrants; t*ones(size(indivs,1),1), indivs(:,3:4)];  
%        migrants(migrants(:,3)>npatch,:) = [];
%     end
    
    % Fitness
    [indivs] = fitnessD(indivs,npatch,r0,eta,cutoff);

if (fit_dd_slope > 0)%|| fit_dd_min > 0 || fit_pos_dd_min > 0 || fit_pos_dd_slope > 0)    
   [indivs] = parasitoids2(indivs,npatch,fit_dd_slope,fit_dd_min,fit_pos_dd_min,fit_pos_dd_slope);
end

%     % set fitness of females colonizing patches with 0 or 1 other female equal to 0 (they
%     % "leave" the patch becasue there are not enough females
%     n = hist(indivs(indivs(:,2)==1,4),1:npatch+1);
%     few_females = find(n(1:npatch)>0 & n(1:npatch) < 3);
%     disp(['number of lonely females after migr = ' num2str(sum(n(few_females))) ' in ' num2str(length(few_females)) ' patches'])
%     disp(['occupancy after migr = ', num2str(size(unique(indivs(:,4)),1)/391)])
%     for i = few_females
%         indivs(indivs(:,4)==i & indivs(:,2)==1,6) = 0;
%     end
%[indivs] = reroute_females(indivs,wij,kap,npatch);
    
    % Inbreeding
    [indivs,flag] = inbreeding(indivs,k,inbreed,t,families);
    
    %rand(1)
    
    if flag == 1
        break
    end
    
    % population ceiling
    [indivs] = ceiling(indivs,npatch,capacities,metapopCV,cap_scaling,min_vbl);
    nKids_out = [nKids_out;indivs(indivs(:,2)==1,7)];
    if (save_migs == 1 && t >= save_indivs_start)
        migrants = [migrants; t*ones(size(indivs,1),1), indivs(:,[3:4,7])];
        migrants(migrants(:,3)>npatch,:) = [];
        migrants(migrants(:,4)<=0,:) = [];
    end

    % Reproduction 
    [indivs,flag,families,families_1older] = diploid_reproduction(indivs,nColInd,max_id,t,families,families_1older,n_L2,achiasma,n_chrom,SNP1,SNP2);

    if flag == 1
        break
    end
    
    % remove the individuals in the death patch
    [indivs,families] = remove_dead(indivs,npatch,families); % can also input/output families
    
    
    % random catastrophes wiping out populations
    [indivs,flag,families] = catastrophes(indivs,catast,families);
    
    if flag == 1
        break
    end 
    
%     % counting the number of offspring per successful mothers
%     [counts] = unique_mothers(families,counts)
    
    % after selection Census
    if (t > burnin  && mod(t,tSkip) == 0 && census_time == 3)
        opThis = add_specs(opThis,indivs,studyPatches,t,npatch,SNP1,SNP2,n_L2,seedi);
    end
    
     display(['t =' num2str(t)])
     display(['number of indivs = ' num2str(size(indivs,1))])
%     display(['occupancy = ' num2str(size(unique(indivs(:,4)),1)/391)])

% output opThis %
if (outPutData == 2 && mod(t,20) == 0)
    %[offsp_counts,xbin] = hist(counts,0:max(counts));
    output_filename = sprintf('opThis_nMC%1i_T%1i_t%1i_aland%1i',nMC,T,t,repid);
    patches = patches(1:npatch,:);
    save(output_filename,'opThis','params') % can also save 'families'
end

if (save_indivs == 1 && t >= save_indivs_start)
    cd ..
    output_filename = sprintf('indivs_nMC%1i_T%1i_t%1i_aland%1i',nMC,T,t,repid);
    save(output_filename,'indivs','t','families') 
    cd('corecode')
end

end % Time loop closes

%if ~isnan(opThis)
    opThis(opThis(:,3)==0,:) = []; % drop rows with no pop
%end

sprintf(['t = ',num2str(t),', simulation ended'])

% output opThis %
if outPutData == 1
    cd ..
    %[offsp_counts,xbin] = hist(counts,0:max(counts));
    output_filename = sprintf('opThis_nMC%1i_T%1i_aland%1i',nMC,T,repid);
    patches = patches(1:npatch,:);
    save(output_filename,'opThis','params','Omega','patches','nKids_out','migrants') % can also save 'families'
end
