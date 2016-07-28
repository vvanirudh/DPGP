function  rmsErrorTmp = predErrorTrajs(trajs, t_start, t_end, numBest)
if nargin == 3
    numBest = 1
end
n = numel(trajs.data)

count = zeros(500,1);
rmsErrorTmp = NaN(n+1,500);
for i = 1:n
    i
    tmpErr = predError(trajs.data(i), t_start, t_end, numBest)
    [~, tmplength] = size(tmpErr);
    count(1:tmplength) = count(1:tmplength) + 1;
    rmsErrorTmp(i,1:tmplength) = tmpErr(1,:);
end
rmsErrorTmp(end,:) = count;

filename = sprintf('predictedError_%d', numBest);
save(filename,'rmsErrorTmp');

% for i = 1:length(count)
%     if count(i)<5
%        max_length = i-1;
%        break; 
%     end
% end
% 
% rmsErrors = zeros(3, max_length);
% rmsErrors(3,:) = 0:0.5:(max_length-1)*0.5;
% 
% for i = 1:max_length
%    errorVec = sort(rmsErrorTmp(:,i));
%    %errorVec = rmsErrorTmp(:,i);
%    rmsErrors(1,i) = mean(errorVec(1:floor(count(i))));
%    rmsErrors(2,i) = std(errorVec(1:floor(count(i))));
% end
% 
% save('predictedError','rmsErrors');

end