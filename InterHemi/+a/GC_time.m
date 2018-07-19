function [gc,p,s] = GC_time(X,params)

u.v2struct(params);

% VAR Model Estimation
ptic('\n*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(X,morder,regmode);
ptoc;

% Autocovariance
ptic('*** var_to_autocov... ');
[G,info] = var_to_autocov(A,SIG,momax);
ptoc;

% Calculate time-domain pairwise-conditional causalities - this just requires
% the autocovariance sequence.

ptic('*** autocov_to_pwcgc... ');
F = autocov_to_pwcgc(G);
ptoc;

% Check for failed GC calculation. If so, return NaN
if(isbad(F,false))
    gc = nan(size(A,1),size(A,2));
    p = nan(size(A,1),size(A,2));
    s = zeros(size(A,1),size(A,2));
    return;
end
% assert(~isbad(F,false),'GC calculation failed');

% Significance test using theoretical null distribution, adjusting for multiple
% hypotheses.

pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,statt); % take careful note of arguments!
sig  = significance(pval,alpha_sig,mhtc);

% Plot time-domain causal graph, p-values and significance.

% figure(fig); clf;
% subplot(3,3,(win-1)*3+1);
% plot_pw(F);
% title('Pairwise-conditional GC'); colorbar;
% subplot(3,3,(win-1)*3+2);
% plot_pw(pval);
% title('p-values'); colorbar;
% subplot(3,3,(win-1)*3+3);
% plot_pw(sig);
% title(['Significant at p = ' num2str(alpha_sig)])

% For good measure we calculate Seth's causal density (cd) measure - the mean
% pairwise-conditional causality. We don't have a theoretical sampling
% distribution for this.

cd = mean(F(~isnan(F)));

fprintf('\ncausal density = %f\n',cd);

gc= F; p = pval; s= sig;