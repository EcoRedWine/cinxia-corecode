function [wij] = calculate_distanceD(patches,npatch,alfa,migDeath,eps,zeta_em,zeta_im,transSIN_mat)

% dispersal probability matrix
[X,Y] = meshgrid(patches(:,2),patches(:,3));
area_dist = repmat(patches(:,1).^zeta_im,1,npatch).*exp(-alfa*sqrt((X-X').^2+(Y'-Y).^2).*transSIN_mat); % proportional to disp prob from patch j to i
area_dist = area_dist.*(1-eye(size(area_dist)));

% stay put is a function of area (Wahlberg et al 2002 Oecologia)
area_modified_eps = eps*patches(:,1).^-zeta_em;
area_modified_eps(area_modified_eps > 1) = 1;

% % migDeath (Eqn 2 from Heino & Hanski 2001 Am Nat)--Does not make sense
% % to me. Why did they sum Sj? Probably to calculate survival for each
% % patch (The better the patches are connected, the more emigrate and suffer
% % from mortality?
% Sj = sum(area_dist,1);
% migSurv = Sj.^2./(20+Sj.^2); % survival probability: with lambda = 30, mean is about 0.78
% phi_jk = (area_dist./repmat(Sj,npatch,1)).*repmat(migSurv,npatch,1);
% wij = repmat(area_modified_eps',npatch,1).*phi_jk; % basically replace (1-migDeath) with migSurv
% wij(logical(eye(size(wij)))) = (1-area_modified_eps);
% wij = [wij;(area_modified_eps)'.*(1-migSurv).*ones(1,npatch)]; % to put ones died during migration
% wij(isnan(wij)) = 0;

% % migDeath (Eqn 2 from Heino & Hanski 2001 Am Nat) -- Sj specific to k
% % destination
% Sjk = area_dist;
% migSurv = Sjk.^2./(0.01+Sjk.^2); %migSurv = migSurv./repmat(sum(migSurv,1),npatch,1);
% phi_jk = (area_dist./repmat(sum(Sjk,1),npatch,1));
% wij = repmat(area_modified_eps',npatch,1).*phi_jk.*migSurv; % basically replace (1-migDeath) with migSurv
% wij(logical(eye(size(wij)))) = (1-area_modified_eps);
% wij = [wij;(1-sum(wij,1))]; % to put ones died during migration
% wij(isnan(wij)) = 0;

wij = (1-migDeath)*repmat(area_modified_eps',npatch,1).*area_dist./repmat(sum(area_dist,1),npatch,1);
wij(logical(eye(size(wij)))) = (1-area_modified_eps); 
wij = [wij;(area_modified_eps)'.*migDeath.*ones(1,npatch)]; % to put ones died during migration
wij(isnan(wij)) = 0;