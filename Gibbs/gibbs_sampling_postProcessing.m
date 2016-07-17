% cluster
% burn_in
% splicing

function [avgSample, mode, config, config_count] = gibbs_sampling_postProcessing(cluster, burn_in, splicing)
% % burn in
cluster_p = cluster(burn_in+1:end,:);

% splicing
cluster_p = cluster_p(1:splicing:end,:);

% reordering
[m n] = size(cluster_p);
for i = 1:m
   cluster_p(i,:) = reorder_vec(cluster_p(i,:));
end

% find the mode
config = [];
config_count = [];
cluster_pp = cluster_p;

mm = m;

while mm > 0

   config = [config; cluster_pp(1,:)];

   diff = max(abs(cluster_pp - repmat(config(end,:),mm,1))');

   ind = find(diff ~= 0);

   cluster_pp = cluster_pp(ind,:);
   config_count = [config_count; mm-length(ind)];
   [mm, n] = size(cluster_pp);
end
[value, value_ind] = max(config_count);
mode = config(value_ind,:);
avgSample = averageSample(cluster_p);
end

function vec_p = reorder_vec(vec)
    n = length(vec);
    vec_p = 9999 * ones(1,n);
    num = 1;
    for i = 1:n
       if vec_p(i) == 9999;
          ind = find(vec == vec(i));
          vec_p(ind) = num;
          num = num + 1;
          
       end
    end   
end