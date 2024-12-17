function [thr,bic,stats] = GaMRed_hist(x,y,K,draw,SW)
%GaMRed  Estimating noise components using Gaussian mixture model
% Input:
% x,y - binned data
% K - number of Gaussian components
% K_noise - number of noise components from left
% draw - if draw results
% Output:
% thr - noise threshold
% bic - Bayesian Information Criterion for estimated model
% stats - additional statistics
%
% Author: Michal Marczyk
% Michal.Marczyk@polsl.pl

if nargin < 3
    error('Insufficient number of arguments.')
end
if nargin < 4
    draw = 0;
end

[x,ind] = sort(x);
y = y(ind);
N = sum(y);       %nb of measurements
bic = Inf;

% remove white or black regions
y(1) = 0;
y(end) = 0;

% remove no signal at the beginning or the end
ind = find(y>0);
ind = ind(1):ind(end);
x = x(ind);
y = y(ind);

%initial conditions
if K==1
    alpha_init = 1;
    mi_init = mean(x);
    sigma_init = std(x);
else
    [alpha_init,mi_init,sigma_init] = gmm_init_dp_hist(x,y,K);
end

if draw
    disp('Starting values')
    disp(num2str([alpha_init; mi_init; sigma_init]))
end

while bic == Inf || bic == 0
    
    % EM algorithm
    [alpha,mi,sigma,logL] = EM_iter_hist(x,y,alpha_init,mi_init,sigma_init,SW);
    
    %calculating BIC
    bic = -2*logL + (3*K - 1)*log(N);
    if bic == Inf || bic == 0
        disp('EM crash. Repeat calculations.')
    end
end

[mi,ind] = sort(mi);
alpha = alpha(ind);
sigma = sigma(ind);

if draw
    disp('Final values')
    disp(num2str([alpha; mi; sigma]))
end

% find threshold betweeen components
if K == 1
    thr = min(x) - 1e-10;
elseif K == 2
    thr = find_thr(x,alpha,mi,sigma,[0;1],draw);
else
    temp = [alpha;mi;sigma]';
    idx = kmeans(temp,2,'emptyaction','singleton','replicates',50);
    thr = find_thr(x,alpha,mi,sigma,idx-1,draw);
end

if ~exist('thr','var')
    thr = NaN;
end

stats.thr = thr;
stats.alpha = alpha;
stats.mu = mi;
stats.K = K;
stats.sigma = sigma;
stats.logL = logL;

end %end function

function [pp_est,mu_est,sig_est,logL] = EM_iter_hist(x,y,alpha,mu,sig,SW)

x = x(:); y = y(:)'; alpha = alpha(:)'; mu = mu(:)'; sig = sig(:)';

N = length(y);  % no. of x values
n = sum(y);     %no. of measurements
sig2 = max(sig.^2, SW^2);
change = Inf;
count = 1;
SW = SW^2;     %minimum variance
eps_change = 1e-6;
KS = length(alpha);
while change > eps_change && count < 10000
    old_alpha = alpha;
    old_sig2 = sig2;
    
    f = zeros(KS,N); sig = sqrt(sig2);
    for a=1:KS
        f(a,:) = norm_pdf(x,mu(a),sig(a));
    end
    px = alpha * f;
    px(isnan(px) | px==0) = 5e-324;
    for a=1:KS
        pk = ((alpha(a)*f(a,:)).*y)./px;
        denom = sum(pk);
        mu(a) = (pk*x)/denom;
        sig2num = sum(pk*((x-mu(a)).^2));
        sig2(a) = max(SW,sig2num/denom);
        alpha(a) = denom/n;
    end
    
    change = sum(abs(alpha-old_alpha)) + sum(((abs(sig2-old_sig2))./sig2))/(length(alpha));
    count = count+1;
end

% RETURN RESULTS
logL = sum(log(px).*y);
[mu_est,ind] = sort(mu);
sig_est = sqrt(sig2(ind));
pp_est = alpha(ind);

end %end function

function y = norm_pdf(x,mu,sigma)
y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);
end

function thr = find_thr(data,alpha,mi,sigma,idx,draw)
% alpha,mi,sigma - components parameters
% idx index for informative/non-informative components

idx = logical(idx);

%generate data with better precision
K = length(mi);
f_temp = zeros(1e7,K);
x_temp = linspace(min(data),max(data),1e7)';
for k=1:K; f_temp(:,k) = alpha(k)*norm_pdf(x_temp,mi(k),sigma(k)); end

%find GMM for informative and non-informative components
f1 = sum(f_temp(:,~idx),2);
f2 = sum(f_temp(:,idx),2);

%calculate difference of f1 and f2 and find its global minimum
f_diff = abs(f1-f2);
[~,ind1] = max(f1);
[~,ind2] = max(f2);
[~,ind] = sort(f_diff);

ind(ind < ind1 | ind > ind2) = [];
if isempty(ind)
    [~,ind] = sort(f_diff);
    a = 1;
    thr_ind = ind(a);
    while thr_ind < ind1 || thr_ind > ind2
        thr_ind = ind(a+1);
    end
else
    thr_ind = ind(1);
end
thr = x_temp(thr_ind);

if draw
    figure; subplot(2,1,1);
    plot(x_temp,f1,'g',x_temp,f2,'r')
    lines = findobj(gca,'Type','Line');
    set(lines,'LineWidth',2)
    set(get(gca,'Ylabel'),'FontSize',14)
    xlabel('Variable'); ylabel('Model')
    
    subplot(2,1,2);
    plot(x_temp,f_diff,'r')
    lines = findobj(gca,'Type','Line');
    set(lines,'LineWidth',2)
    set(get(gca,'Ylabel'),'FontSize',14)
    xlabel('Variable'); ylabel('Models difference')
    title(['Threshold:',num2str(thr)]);
end
end  %end function

function [alpha,mu,sigma] = gmm_init_dp_hist(x,y,K)
% GMM_INIT_DP(data,K)
% Compute initial conditions for GMM by using dynamic programming for
% approximate signal (by operation of binning).
% Input:
% x,y - sample to partition
% K - number of partitions [1x1]
% Output:
% alpha - weights
% mu - means
% sigma - standard deviations

%parameters
par1 = 0.1; %for robustness (fine for data in range 0-20)
par2 = 10;   %minimum no. of points in signal fragment

% initialize
s_corr = ((x(2) - x(1))^2)/12;  %sheppards correction for binned data
K = K - 1;
N = length(x);
p_opt_idx = zeros(1,N);
p_aux = zeros(1,N);
opt_pals = zeros(K,N);
for a=1:N
    invec = x(a:N);
    yinvec = y(a:N);
    if sum(yinvec) <= par2
        p_opt_idx(a) = inf;
    else
        wwec=yinvec/(sum(yinvec));
        var_bin = sum(((invec-sum(invec.*wwec)).^2).*wwec);
        if var_bin > s_corr
            p_opt_idx(a)=(par1+sqrt(var_bin-s_corr))/(max(invec)-min(invec));
        else
            p_opt_idx(a) = inf;
        end
    end
end

% aux_mx
aux_mx = zeros(N,N);
for a=1:N-1
    for b=a+1:N
        invec = x(a:b-1);
        yinvec = y(a:b-1);
        if sum(yinvec) <= par2
            aux_mx(a,b) = inf;
        else
            wwec = yinvec/(sum(yinvec));
            var_bin = sum(((invec-sum(invec.*wwec)).^2).*wwec);
            if var_bin > s_corr
                aux_mx(a,b) = (par1+sqrt(var_bin-s_corr))/(max(invec)-min(invec));
            else
                aux_mx(a,b) = inf;
            end
            
        end
    end
end

%iterate
for kster = 1:K
    % kster
    for a=1:N-kster
        for b=a+1:N-kster+1
            p_aux(b) =  aux_mx(a,b) + p_opt_idx(b);
        end
        [mm,ix] = min(p_aux(a+1:N-kster+1));
        p_opt_idx(a) = mm;
        opt_pals(kster,a) = a + ix(1);
    end
end

%restore optimal decisions
opt_part = zeros(1,K);
opt_part(1) = opt_pals(K,1);
for kster = K-1:-1:1
    opt_part(K-kster+1) = opt_pals(kster,opt_part(K-kster));
end

%find initial conditions
opt_part=[1 opt_part N+1];
alpha = zeros(1,K); mu = alpha; sigma = alpha;
for a=1:K+1
    invec = x(opt_part(a):opt_part(a+1)-1);
    yinvec = y(opt_part(a):opt_part(a+1)-1);
    wwec = yinvec/(sum(yinvec));
    alpha(a) = sum(yinvec)/sum(y);
    mu(a) = sum(invec.*wwec);
    sigma(a)= sqrt(sum(((invec-sum(invec.*wwec)).^2).*wwec)-s_corr);
end
end  %end function