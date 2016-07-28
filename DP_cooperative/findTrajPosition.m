%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to find the position in the trajectory
% corresponding to a given timestep
% by Anirudh Vemula, Jul 21, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, y] = findTrajPosition(traj, t)

% Construct the cumulative time array
cum_t = cumsum(traj.dt);

% Find the index with cum_t > current_time
ind = find(cum_t > t);

% Return the (x,y) position
if ~isempty(ind)
    x = traj.x(ind(1));
    y = traj.y(ind(1));
else
    x = traj.x(end);
    y = traj.y(end);
end

end