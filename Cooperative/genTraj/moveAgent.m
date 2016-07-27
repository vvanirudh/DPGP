%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to move agent using quadratic
% interpolation
% by Anirudh Vemula, Jul 25, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x2_new = moveAgent(x1, x2, ind, threshold)

% Function always assumes that agent 2 has to be moved w.r.t agent 1
x2_new = x2;

% Length of the trajectory
l2 = length(x2);

if x1(ind) > x2(ind)    
    x2_new(ind) = x2(ind) - threshold;        
else
    x2_new(ind) = x2(ind) + threshold;
end

quad = polyfit([ind-3 ind ind+3], [x2_new(ind-3) x2_new(ind) ...
                    x2_new(ind+3)], 2);


x2_new(ind-1) = polyval(quad, ind-1);
x2_new(ind+1) = polyval(quad, ind+1);
x2_new(ind-2) = polyval(quad, ind-2);
x2_new(ind+2) = polyval(quad, ind+2);

end