function simLevel = similarityDist(v1, v2)
    
    num_clus_1 = max(v1);
    num_clus_2 = max(v2);
    n = length(v1);
    
    label1 = zeros(n,num_clus_1);
    label2 = zeros(n,num_clus_2);
    
    for i = 1:num_clus_1
       label1(:,i) = (v1==i); 
    end
    
    for i = 1:num_clus_2
       label2(:,i) = (v2==i); 
    end
    
    tmp = label1' * label2;
    simLevel = 0;
    for i = 1:num_clus_1
        for j = 1:num_clus_2
            simLevel = simLevel + tmp(i,j) * sum(label1(:,i) == label2(:,j));     
        end
    end
end