function sparseGP = build_sparseGP_tau(x_obs, y_obs, tau_obs, speed, ...
                                   hyperparam_0, ifShift)

    
    if (nargin == 5)
        ifShift = 0;
    end
    
    if (ifShift >=1)
        x_obs = [x_obs; x_obs+0.5; x_obs-0.5];
        y_obs = [y_obs; y_obs+0.5; y_obs-0.5];
        tau_obs = [tau_obs; tau_obs+0.5; tau_obs-0.5];
    end
    
    % permutate the order of input trajectory
    reorder = randperm(length(x_obs));
    x_obs = x_obs(reorder);
    y_obs = y_obs(reorder);
    tau_obs = tau_obs(reorder);
    speed = speed(reorder);
    
    sparseGP.x_data = x_obs;
    sparseGP.y_data = y_obs;
    sparseGP.tau_data = tau_obs;
    sparseGP.speed = speed;
    sparseGP.hyperparam = hyperparam_0;
    
    threshold = 0.1;
    
    n = length(x_obs);
    sparseGP.budget = n;
    sparseGP.numBV = 0;
    
    % default value
    sparseGP.budget = 5;
    
    %% initialization
    gkt = Gaussian_kernel_tau(x_obs(1), y_obs(1), tau_obs(1), x_obs(1), ...
                              y_obs(1), tau_obs(1), ...
                              sparseGP.hyperparam);
    gkt_0 = Gaussian_kernel_tau(x_obs(1), y_obs(1), tau_obs(1), x_obs(1), ...
                              y_obs(1), tau_obs(1), ...
                              sparseGP.hyperparam, 0);
    
    
    alpha = speed(1) / gkt;
    
    C = -1 / gkt;
    Q = 1 / gkt_0;

    %% include into BV
    sparseGP.BV_x = x_obs(1);
    sparseGP.BV_y = y_obs(1);
    sparseGP.BV_tau = tau_obs(1);
    sparseGP.BV_speed = speed(1);
    numBV = 1;
    sigma_noise = sparseGP.hyperparam(5);
    
    %% main loop
    for i=2:n
        %in_bv = [sparseGP.BV_x', sparseGP.BV_y', sparseGP.BV_tau'];
        %in_i = [x_obs(i), y_obs(i), tau_obs(i)];                
        
        kx = Gaussian_kernel_tau(sparseGP.BV_x, sparseGP.BV_y, ...
                                 sparseGP.BV_tau, x_obs(i), y_obs(i), ...
                                 tau_obs(i), sparseGP.hyperparam, 0);
        
        kxx = Gaussian_kernel_tau(x_obs(i), y_obs(i), tau_obs(i), ...
                                  x_obs(i), y_obs(i), tau_obs(i), ...
                                  sparseGP.hyperparam, 0);
                
        Kt_xx = kxx + kx' * C * kx;
        %keyboard()
        q = (speed(i) - alpha'*kx) / (Kt_xx + sigma_noise^2);
        r = -1 / (Kt_xx + sigma_noise^2);
        
        e_hat = Q * kx;
        gamma = kxx - kx' * Q * kx;
        Q = [Q, zeros(numBV,1); zeros(1,numBV),0] + 1 / gamma * ([e_hat; -1]) ...
            * ([e_hat; -1])';
        
        % adding a new BV
        if (numBV < sparseGP.budget || gamma > threshold)
            % update alpha and C      
            s = [C*kx; 1];
            alpha = [alpha; 0] + q * s;
            C = [C, zeros(numBV,1); zeros(1,numBV),0] + r * s*s';
            
            numBV = numBV + 1;
            sparseGP.BV_x(numBV) = x_obs(i);
            sparseGP.BV_y(numBV) = y_obs(i);
            sparseGP.BV_tau(numBV) = tau_obs(i);
            sparseGP.BV_speed(numBV) = speed(i);
            %where is the update for Q???
        % approximate update
        else
            s = C*kx + e_hat;
            alpha = alpha + q * s;
            C = C + r * s*s';
            Q = Q(1:end-1,1:end-1) - (Q(1:end-1,end) * Q(1:end-1,end)') / (Q(end,end));

        end
        
        % remove a BV
        if (numBV > sparseGP.budget)
           % compute score
           epsilon = zeros(numBV,1);
           for j=1:numBV
               epsilon(j) = abs(alpha(j) / Q(j,j));
           end

           [min_value, min_index] = min(epsilon);
           threshold = min_value;
           
           % reduce matrix size
           % swap row and column
           C = [C(1:min_index-1,:); C(min_index+1:end,:); C(min_index,:)];
           C = [C(:,1:min_index-1), C(:,min_index+1:end), C(:,min_index)];
           Q = [Q(1:min_index-1,:); Q(min_index+1:end,:); Q(min_index,:)];
           Q = [Q(:,1:min_index-1), Q(:,min_index+1:end), Q(:,min_index)];
           alpha = [alpha(1:min_index-1); alpha(min_index+1:end); alpha(min_index)];
           % update (order matters!)
           alpha = alpha(1:end-1) - alpha(end) * (Q(1:end-1,end)) / (Q(end,end));
           C = C(1:end-1,1:end-1) + C(end,end) * Q(1:end-1,end) * Q(1:end-1,end)' / Q(end,end)^2 - ...
               1 / Q(end,end) * (Q(1:end-1,end) * C(1:end-1,end)' + C(1:end-1,end) * Q(1:end-1,end)');
           Q = Q(1:end-1,1:end-1) - (Q(1:end-1,end) * Q(1:end-1,end)') / (Q(end,end));
           
           numBV = numBV - 1;
           sparseGP.BV_x = [sparseGP.BV_x(1:min_index-1), sparseGP.BV_x(min_index+1:end)];
           sparseGP.BV_y = [sparseGP.BV_y(1:min_index-1), sparseGP.BV_y(min_index+1:end)];
           sparseGP.BV_tau = [sparseGP.BV_tau(1:min_index-1), sparseGP.BV_tau(min_index+1:end)];
           sparseGP.BV_speed = [sparseGP.BV_speed(1:min_index-1), sparseGP.BV_speed(min_index+1:end)];
           
        end
        
        speed_eff = -C \ alpha;
        
    end
    
    
    sparseGP.speed_eff = speed_eff;
    sparseGP.alpha = alpha;
    sparseGP.numBV = numBV;
    sparseGP.C = C;
    sparseGP.s = s;
    sparseGP.kxx = kxx;
    sparseGP.kx = kx;
    sparseGP.q = q;
    sparseGP.Q = Q;
end
