function rmsError = predError(traj, t_start, t_end, numBest)

if nargin == 3
    numBest = 1;
end

global sparseGPs
    traj_len = length(traj.x);
    t = traj.dt;
    dt = traj.dt;
    t(1) = 0;
    time_fin = sum(dt(1:end-1));
    fin_pre_ind = 0;
    for i = 2:traj_len
       t(i) = t(i-1) + traj.dt(i-1);
       if t(i-1) + t_end < time_fin && t(i) + t_end > time_fin
           fin_pre_ind = i-1;
       end
    end
    
    pos_var_in = [0.5, 0; 0, 0.5];
    rmsError = ones(2,fin_pre_ind);
    tmp_dt_start_ind = 1;
    tmp_dt_end_ind = 1;
    
    
    speed_vec_x = traj.dx_dt;
    speed_vec_y = traj.dy_dt;
    speed_vec_x(3:end) = (speed_vec_x(1:end-2) + speed_vec_x(2:end-1) + speed_vec_x(3:end))/3;
    speed_vec_y(3:end) = (speed_vec_y(1:end-2) + speed_vec_y(2:end-1) + speed_vec_y(3:end))/3;
    speed_vec = sqrt(traj.dx_dt.^2 + traj.dy_dt.^2);
    speed_vec(3:end) = (speed_vec(1:end-2) + speed_vec(2:end-1) + speed_vec(3:end))/3;
    
    
    for i = 1:fin_pre_ind
       tmp_traj = {};
       tmp_traj.x = traj.x(max(1,i-5):i);
       tmp_traj.y = traj.y(max(1,i-5):i);
       tmp_traj.dx_dt = traj.dx_dt(max(1,i-5):i);
       tmp_traj.dy_dt = traj.dy_dt(max(1,i-5):i);
       tmp_traj.dt = traj.dt(max(1,i-5):i);
%        tmp_traj.x = traj.x(1:i);
%        tmp_traj.y = traj.y(1:i);
%        tmp_traj.dx_dt = traj.dx_dt(1:i);
%        tmp_traj.dy_dt = traj.dy_dt(1:i);
%        tmp_traj.dt = traj.dt(1:i);
       if isfield(traj, 'DP_alpha')
          tmp_traj.DP_alpha = traj.DP_alpha;
       end
       [ind_order, ifNew] = findBestPattern(tmp_traj);
       fprintf(sprintf('best index order: '));
       for k = 1:numBest
           if k == numBest && ifNew <= numBest
            fprintf(sprintf('new'));  
           else
            fprintf(sprintf('%d \t', ind_order(k)));
           end
       end
       fprintf('\n');
       
       %length(sparseGP_x.x_data)
       %sparseGP_x.speed

       % find index for starting RMS comparison
       for j = tmp_dt_start_ind:traj_len
           if ((t(j) - t(i)) > t_start)
               tmp_dt_start_ind = j;
               break;
           end
       end

       % find index for final RMS comparison
       for k = tmp_dt_end_ind:traj_len
           if ((t(k) - t(i) < t_end) &&  (t(k+1) - t(i) > t_end))
               tmp_dt_end_ind = max(k,tmp_dt_start_ind);
               break;
           end
       end

       % prediction


       %speed = sqrt(traj.dx_dt(i)^2 + traj.dy_dt(i)^2);
       speed = speed_vec(i);
       
       tmp_dt = dt(i:tmp_dt_end_ind);
       num_steps = length(tmp_dt);
       pos = [traj.x(i), traj.y(i)];
       min_error = 9999 * ones(1,numBest);
       for k = 1:numBest
           ind_tmp = ind_order(k);
           sparseGP_x = sparseGPs(ind_tmp).sparseGP_x;
           sparseGP_y = sparseGPs(ind_tmp).sparseGP_y;
           if k == numBest && ifNew <= numBest
               speed_x = speed_vec_x(i);
               speed_y = speed_vec_y(i);
               pred_traj = linearTrajPred(pos, pos_var_in, tmp_dt, num_steps, speed_x, speed_y);
           else
               pred_traj = sparseGP_trajPredict(sparseGP_x, sparseGP_y, pos, pos_var_in, tmp_dt, num_steps,speed);
           end
           tmpError =  sum((traj.x(tmp_dt_start_ind: tmp_dt_end_ind) - pred_traj(tmp_dt_start_ind-i+1:end-1,1)).^2) + ...
           sum((traj.y(tmp_dt_start_ind: tmp_dt_end_ind) - pred_traj(tmp_dt_start_ind-i+1:end-1,2)).^2);
           min_error(k) = sqrt(tmpError / (tmp_dt_end_ind - tmp_dt_start_ind + 1));
       end
       [min_value , min_ind] = sort(min_error, 'ascend');
       rmsError(1,i) = min_value(1);
       rmsError(2,i) = t(i);
       sparseGP_x = sparseGPs(ind_order(min_ind(1))).sparseGP_x;
       sparseGP_y = sparseGPs(ind_order(min_ind(1))).sparseGP_y;
       if ifNew <= numBest && min_ind(1) == numBest
           pred_traj = linearTrajPred(pos, pos_var_in, tmp_dt, num_steps, speed_x, speed_y);
           fprintf('best index (true among candidates) : new \n'); 
       else
           pred_traj = sparseGP_trajPredict(sparseGP_x, sparseGP_y, pos, pos_var_in, tmp_dt, num_steps,speed);
           fprintf(sprintf('actual best index among candidates : %d\n', ind_order(min_ind(1)))); 
       end
    end
    
      predPlot(pred_traj, sparseGP_x, sparseGP_y);
      hold on
      plot(traj.x, traj.y,'k-^');
      plot(traj.x(1), traj.y(1),'m*');
      rmsError

    
    % interpolate to 0.5s increment
    t_inc = 0.5;
    [~,lengthRms] = size(rmsError);
    if lengthRms == 1
        return 
    end

    t_int_ind = floor(rmsError(2,end-1)/t_inc);
     %t_int_ind = 11;
    t_new = 0:t_inc:t_inc*t_int_ind;
    rmsErrorNew = zeros(2,t_int_ind+1);
    rmsErrorNew(2,:) = t_new;
    rmsErrorNew(1,1) = rmsError(1,1);
    cur_ind = 1;
    for i = 2:t_int_ind+1
       for j =  cur_ind:fin_pre_ind-1
          if t_new(i) > t(j) &&  t_new(i) <= t(j+1)
              cur_ind = j;
              rmsErrorNew(1,i) = rmsError(1,j) + (rmsError(1,j+1)-rmsError(1,j)) * ...
                  (t_new(i) - t(j)) / (t(j+1) - t(j));
          end
       end    
    end
    rmsError = rmsErrorNew;
end































