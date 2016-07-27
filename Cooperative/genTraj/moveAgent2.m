%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to try a different kind of collision
% avoidance
% by Anirudh Vemula, Jul 26, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x2_new = moveAgent2(x1, x2, ind, threshold)

% Function always assumes that agent 2 has to be moved w.r.t agent 1    
x2_new = x2;

% Length of the trajectory
l2 = length(x2);

if x1(ind) > x2(ind)
    x2_new(ind) = x2(ind) - threshold;
else
    x2_new(ind) = x2(ind) + threshold;
end

quad = polyfit([ind-3 ind l2-1 l2], [x2_new(ind-3) x2_new(ind) ...
                    x2_new(l2-1) x2_new(l2)], 3);

x2_new(ind-3:end) = polyval(quad, ind-3:l2);

end
