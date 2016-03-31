function [r0,gamma,environSD,locSD,eps,migDeath,alfa,kap,n_L1,n_L2,ratio,migr_opt,MR1,MR2,K_area_scaling,eta,zeta_em,encounter,handling,metapopCV,achiasma,n_chrom,k,cutoff,inbreed,zeta_im,catast,fit_dd_slope,fit_dd_min,fit_pos_dd_min,fit_pos_dd_slope,cap_scaling,min_vbl] = fix_params(params)

r0 = params(1); % r0 > 0
gamma = params(2);
var_addit = params(3); % additive variance
h2 = 0.5; %heritability
environSD = sqrt(var_addit*(1-h2)/h2); % SD of residual effect???
locSD = params(4); 
eps = params(5); % migration probability
migDeath= params(6);
alfa = params(7); % small = long range migration
kap = params(8);
n_L1 = params(9);
n_L2 = params(10); % #SNP/inbreeding loci (when used as neutral SNPs, set s = 0)
ratio = params(11); % patch area to carrying cap ratio
migr_opt = params(12); % migration option
MR1 = params(13); % mutation rate for autosomes
MR2 = params(14); % mutation rate for inbreeding loci
K_area_scaling = params(15); % patch area scaling for carrying cap
eta = params(16); % called nu in R code, also act as sigma_r0
zeta_em = params(17); % area scaling parameter for emigraion
zeta_im = params(18); % area scaling parameter for immigraion
encounter = params(19); % encounter rate with females
handling = params(20); % handling time for mating
metapopCV = params(21); % 0.15, SD of metapopulation size
achiasma = params(22); % %0 = yes recombination in females, 1 = no recomb in females
n_chrom = params(23); % number of chromosomes
k = params(24); % inbreeding penalty
cutoff = params(25); % max female fitness
inbreed = params(26); % 0 = inbreeding off
catast = params(27); % rate of random destruction
fit_dd_slope = params(28); % rate or fitness decay with pop dens
fit_dd_min = params(29); % minimum of the function
fit_pos_dd_min  = params(30); % minimum of the pos dd function
fit_pos_dd_slope = params(31); % slope of the pos dd function
cap_scaling = params(32); % the scaling exponent of carrying cap for ceiling.m
min_vbl = params(33); % minimum viable butterfly number in spring (min #butts expected from one successful nest)

if mod(n_L2, n_chrom) ~= 0
    error('n_L2 needs to be a multiple of n_chrom')
end

