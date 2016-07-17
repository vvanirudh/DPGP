function dist = Gaussian_kernel(x1,y1,x2,y2, hyperparam, ifAddNoiseInput)

    lx = hyperparam(1);
    ly = hyperparam(2);
    sigma_input = hyperparam(3);
    sigma_noise = hyperparam(4);

    dist = zeros(length(x1), length(x2));
    
    if (nargin == 6)
        ifAddNoise = ifAddNoiseInput;
    else 
        ifAddNoise = 1;
    end

    % matrix version
%     [X1,X2] = meshgrid(x1, x2);
%     [Y1,Y2] = meshgrid(y1, y2);
%     X1 = X1';
%     X2 = X2';
%     Y1 = Y1';
%     Y2 = Y2';
%     
% 
%     dist = sigma_input^2 * exp( -(X1-X2).^2 / (2*lx^2) ...
%                                 -(Y1-Y2).^2 / (2*ly^2));
% 
%     if (ifAddNoise == 1)
%         min_n = min(length(x1), length(x2));   
%         same_entry = (x1(1:min_n) == x2(1:min_n)) .* (y1(1:min_n) == y2(1:min_n));
%         dist(1:min_n, 1:min_n) = dist(1:min_n, 1:min_n) + sigma_noise^2 * diag(same_entry); 
%     end
      
    % loop version
    for i = 1:length(x1)
        for j = 1: length(x2)
                dist(i,j) = sigma_input^2 * exp(-(x1(i)-x2(j))^2 / (2*lx^2) ...
                                                -(y1(i)-y2(j))^2 / (2*ly^2));
                if (x1(i) == x2(j) && y1(i) == y2(j) && ifAddNoise == 1)
                    dist(i,j) = dist(i,j) + sigma_noise^2;
                end
        end
    end

end