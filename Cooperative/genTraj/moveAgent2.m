%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to try a different kind of collision
% avoidance
% by Anirudh Vemula, Jul 26, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x2_new = moveAgent2(x1, x2, ind, threshold)

x2_new = x2;

if x1(ind) > x2(ind)
    x2_new(ind) = x2(ind) - threshold;
else
    x2_new(ind) = x2(ind) + threshold;
end

quad = polyfit([ind-2 ind-1 ind], [x2_new(ind-2) x2_new(ind-1) ...
                    x2_new(ind)], 2);
end
