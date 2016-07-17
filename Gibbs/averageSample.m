function sample = averageSample(cluster)

% each row is a column
[num_sample, num_traj] = size(cluster);

similarity_matrix = zeros(num_sample);

for i = 1:num_sample
    for j = 1:i
        similarity_matrix(i,j) = similarityDist(cluster(i,:)', cluster(j,:)');
        similarity_matrix(i,j) = similarity_matrix(j,i);
    end
end

% best configuration
similarity_score = sum(similarity_matrix, 1);

[value, ind] = max(similarity_score);

% throw away the bottom 30% samples
reduced_num_samples = floor(0.5 * num_sample);
[order, order_ind] = sort(similarity_score);
reduced_sim_matrix = similarity_matrix(order_ind(1:reduced_num_samples),...
    order_ind(1:reduced_num_samples));
reduced_sim_score = sum(reduced_sim_matrix,1);
[value, ind] = max(reduced_sim_score);

sample = cluster(order_ind(ind),:);
end