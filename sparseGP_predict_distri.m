% current distribution: N(pos, pos_var)
% model: sparseGP
% return: prediction after dt
function [mu, var]= sparseGP_predict_distri(sparseGP, pos, pos_var)


W = [sparseGP.hyperparam(1)^2, 0; 0, sparseGP.hyperparam(2)^2];
W_inv = [sparseGP.hyperparam(1)^(-2), 0; 0, sparseGP.hyperparam(2)^(-2)];

I = eye(2);

numBV = sparseGP.numBV;
pos_basis = [sparseGP.BV_x; sparseGP.BV_y]';


hyperparam = sparseGP.hyperparam;
kx = Gaussian_kernel(pos(1), pos(2), sparseGP.BV_x, sparseGP.BV_y, hyperparam,0);
kxx = Gaussian_kernel(pos(1), pos(2), pos(1), pos(2), hyperparam,0);
beta = sparseGP.alpha;


Q = ones(numBV);
q = ones(numBV,1);
K = Gaussian_kernel(sparseGP.BV_x, sparseGP.BV_y, sparseGP.BV_x, sparseGP.BV_y, hyperparam,1);
%K
%sparseGP.BV_speed
%beta = K\sparseGP.BV_speed';



for i = 1:numBV
    q(i) = exp(-0.5 * (pos - pos_basis(i,:))  * (inv(pos_var + W) * (pos - pos_basis(i,:))'));
    
    for j = i:numBV
        xb = (pos_basis(i,:) + pos_basis(j,:)) / 2.0;
        x_diff = (pos_basis(i,:) - pos_basis(j,:));
    
        Q(i,j) = exp(-0.5 * (xb - pos) * (inv(0.5*W+pos_var) * (xb - pos)')) ...
        * exp(-0.5 * x_diff * 0.5*W_inv * x_diff');
        Q(j,i) = Q(i,j);
    end
end
q = det(W_inv * pos_var + I) ^ -0.5 * q;
Q = det(2 * W_inv * pos_var + I) ^ -0.5 * Q;



mu = q'* beta;
%var = kxx + trace((beta * beta' - sparseGP.Q) * Q) - mu^2;
var = kxx + trace((beta * beta' + sparseGP.C) * Q) - mu^2;

end