function dist = Neural_Network_kernel(tau1, tau2, hyperparam, ...
                                      ifAddNoiseInput)

    sigma0 = hyperparam(5);    
    sigma = hyperparam(6);
    sigma_noise = hyperparam(4);
    
    dist = zeros(length(tau1), length(tau2));
    
    if (nargin==4)
        ifAddNoise = ifAddNoiseInput;
    else
        ifAddNoise = 1;
    end    
    
    Sigma = diag([sigma0^2; sigma^2], 0);
    
    %loop version
    for i=1:length(tau1)
        for j=1:length(tau2)
            tau1_aug = [1; tau1(i)]; tau2_aug = [1; tau2(j)];
            
            dist(i,j) = 2/pi * ((2*tau1_aug'*Sigma*tau2_aug) / ...
                                sqrt((1+2*tau1_aug'*Sigma* ...
                                      tau1_aug)*(1+2*tau2_aug'* ...
                                                 Sigma*tau2_aug)));
            
            if (tau1(i)==tau2(j) && ifAddNoise==1)
                dist(i,j) = dist(i,j) + sigma_noise^2;
            end
            
        end
    end

end
