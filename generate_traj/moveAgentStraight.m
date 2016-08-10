%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to try a different kind of collision
% avoidance
% by Anirudh Vemula, Aug 6, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x2_new = moveAgentStraight(x1, x2, ind, threshold)

% Function always assumes that agent 2 has to be moved w.r.t agent 1    
x2_new = x2;

% Length of the trajectory
l2 = length(x2);

if x1(ind) > x2(ind)
    x2_new(ind:end) = x2(ind:end) - threshold;
    
else
    x2_new(ind:end) = x2(ind:end) + threshold;
end

quad = polyfit([ind-5 ind-3 ind], [x2_new(ind-5) x2_new(ind-3) ...
                    x2_new(ind)], 2);

x2_new(ind-5:ind) = polyval(quad, ind-5:ind);

end
