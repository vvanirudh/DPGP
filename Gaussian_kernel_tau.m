function dist = Gaussian_kernel_tau(in1, in2, hyperparam, ...
                                    ifAddNoiseInput)

    lx = hyperparam(1);
    ly = hyperparam(2);
    ltau = hyperparam(3);
    
    sigma_input = hyperparam(4);
    sigma_noise = hyperparam(5);
    
    if (nargin == 4)
        ifAddNoise = ifAddNoiseInput;
    else
        ifAddNoise = 1;
    end
    
    x1 = in1(:,1); y1 = in1(:, 2); tau1 = in1(:,3);
    x2 = in2(:,1); y2 = in2(:, 2); tau2 = in2(:,3);
    
    for i=1:length(x1)
        for j=1:length(x2)
            
            dist(i,j) = sigma_input^2 * exp( - (x1(i) - x2(j))^2 / ...
                                             (2*lx^2) - (y1(i) - ...
                                                         y2(j))^2 / ...
                                             (2*ly^2) - (tau1(i) - ...
                                                         tau2(j))^2 ...
                                             / (2*ltau^2));
            
            if (x1(i) == x2(j) && y1(i) == y2(j) && tau1(i) == ...
                tau2(j))
                dist(i,j) = dist(i,j) + sigma_noise^2;
            end
            
        end
    end

end
