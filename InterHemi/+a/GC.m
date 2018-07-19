function [gc,f] = GC(X,params)

u.v2struct(params);

% VAR Model Estimation
ptic('\n*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(X,morder,regmode);
ptoc;

% Autocovariance
ptic('*** var_to_autocov... ');
[G,~] = var_to_autocov(A,SIG,acmaxlags);
ptoc;

% Time domain Granger Causality 
% ptic('*** autocov_to_pwcgc... ');
% F = autocov_to_pwcgc(G);
% ptoc;
% figure;
% plot_pw(F);
% keyboard

% Significance test using theoretical null distribution, adjusting for multiple
% hypotheses.

% NEXT LINE KEEPS BREAKING FOR SOME REASON. TSTAT?
% pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
% sig  = significance(pval,alpha,mhtc);

% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)

% Calculate spectral pairwise-conditional causalities at given frequency
% resolution - again, this only requires the autocovariance sequence.

ptic('\n*** autocov_to_spwcgc... ');
f = autocov_to_spwcgc(G,fres);
ptoc;

% Check for failed spectral GC calculation

assert(~isbad(f,false),'spectral GC calculation failed');

% Plot spectral causal graph.
fres = size(f,3)-1;
lam = sfreqs(fres,fs)';
idx = lam >= frange(1) & lam <= frange(2);
lam = lam(idx);
f = f(:,:,idx);

gc(1,:) = squeeze(f(2,1,:));
gc(2,:) = squeeze(f(1,2,:));

f = lam;
% figure(3); clf;
% plot_spw(f,fs,[0,100]);

end