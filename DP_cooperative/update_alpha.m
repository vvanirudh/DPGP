function alpha = update_alpha (k, n)

% inverse gamma prior: a^(-3/2) exp (-1/(2a))
% likelihood: a^k *Gamma(a) / Gamma(n+a) 
% posterior: a^(k-3/2) * exp(-1/2a) * Gamma (a) / Gamma(n+a)
a = 0.01:0.01:3.0;
% p_a = a.^(k-9/2-1) .* exp(-1./ (8*a)) .* gamma(a) ./gamma(a+n)
% for numerical accuracy 
% note gamma(a+n) = gamma(a) a (a+1) ... (a+n-1)
% implies gamma(a+n) / gamma(n) = gamma(a) * prod ((a+i-1)/i) for i=1:n
n_tmp = 1:n;
[A, N_tmp] = meshgrid(a,n_tmp);
A_add_N = (N_tmp + A -1) ./ N_tmp;

p_a = a.^(k-9/2-1) .* exp(-1./ (8*a)) ./ prod(A_add_N);
p_a = p_a / sum(p_a);
%plot (a,p_a,'r');

sample = randsample(length(a), 1, true, p_a);
alpha = a(sample);
%figure;
%hist(samples);
% alpha = 0.5;
end