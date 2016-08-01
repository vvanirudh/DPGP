function dist = Gaussian_kernel_two(p11, p21, p12, p22, hyperparam, ...
                                    ifAddNoiseInput)

lx1 = hyperparam(1);
lx2 = hyperparam(2);

ly1 = hyperparam(3);
ly2 = hyperparam(4);

sigma_input = hyperparam(5);
sigma_noise = hyperparam(6);

dist = zeros(size(p11, 1), size(p12, 1));

if (nargin == 6)
    ifAddNoise = ifAddNoiseInput;
else
    ifAddNoise = 1;
end

for i=1:size(p11, 1)
    for j = 1:size(p12, 1)
        
        dist(i,j) = sigma_input^2 * exp(-(p11(i,1) - p12(i,1))^2 / ...
                                        (2*lx1^2) - (p21(i,1) - ...
                                                     p22(i,1))^2 / ...
                                        (2*lx2^2) - (p11(i,2) - ...
                                                     p12(i,2))^2 / ...
                                        (2*ly1^2) - (p21(i,2) - ...
                                                     p22(i,2))^2 / ...
                                        (2*ly2^2));
        if (p11(i,1) == p12(i,1) && p11(i,2) == p12(i,2) && p21(i,1) ...
            == p22(i,1) && p21(i,2)==p22(i,2))
            dist(i,j) = dist(i,j) + sigma_noise^2;
        end
        
    end
end

end