
numPrediction = 1;
filename_DPAR = sprintf('../DP_expt/predictedError_%d.mat', numPrediction);
filename_DPGP = sprintf('../DP/predictedError_%d.mat', numPrediction);
filename_label = sprintf('../DP_expt/predictedError_label_%d.mat', numPrediction);
exptRmsErrorRaw_tmp = load(filename_DPAR );
DPGPRmsErrorRaw_tmp = load(filename_DPGP);
exptRmsErrorRaw = exptRmsErrorRaw_tmp.rmsErrorTmp(1:end-1,:);
DPGPRmsErrorRaw = DPGPRmsErrorRaw_tmp.rmsErrorTmp(1:end-1,:);

exptRmsErrorRaw_label_tmp = load(filename_label);
exptRmsErrorRaw_label = exptRmsErrorRaw_label_tmp.rmsErrorTmp(1:end-1,:);

count = exptRmsErrorRaw_tmp.rmsErrorTmp(end,:);
% remove outliers (where both clustering performed poorly)
[num_traj, traj_length] = size(exptRmsErrorRaw);
% for i = 1:num_traj
%     for j = 1:traj_length
%        if isnan(exptRmsErrorRaw(i,j))
%            break
%        end
%        
%        if exptRmsErrorRaw(i,j) > 3 && DPGPRmsErrorRaw(i,j) > 3
%           exptRmsErrorRaw(i,j) = NaN;
%           DPGPRmsErrorRaw(i,j) = NaN;
%           count(j) = count(j) - 1;
%        end
%     end 
% end


for i = 1:length(count)
    if count(i)<5
       max_length = i-1;
       break; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exptRmsError = zeros(6, max_length);
exptRmsError(6,:) = 0:0.5:(max_length-1)*0.5;
for i = 1:max_length
   %errorVec = sort(rmsErrorTmp(:,i));
   errorVec = sort(exptRmsErrorRaw(:,i));
   exptRmsError(1,i) = mean(errorVec(1:floor(count(i))));
   exptRmsError(2,i) = std(errorVec(1:floor(count(i))));
   exptRmsError(3,i) = prctile(errorVec(1:floor(count(i))),50);
   exptRmsError(4,i) = prctile(errorVec(1:floor(count(i))),25);
   exptRmsError(5,i) = prctile(errorVec(1:floor(count(i))),75);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exptRmsError_label = zeros(6, max_length);
exptRmsError_label(6,:) = 0:0.5:(max_length-1)*0.5;
for i = 1:max_length
   %errorVec = sort(rmsErrorTmp(:,i));
   errorVec = sort(exptRmsErrorRaw_label(:,i));
   exptRmsError_label(1,i) = mean(errorVec(1:floor(count(i))));
   exptRmsError_label(2,i) = std(errorVec(1:floor(count(i))));
   exptRmsError_label(3,i) = prctile(errorVec(1:floor(count(i))),50);
   exptRmsError_label(4,i) = prctile(errorVec(1:floor(count(i))),25);
   exptRmsError_label(5,i) = prctile(errorVec(1:floor(count(i))),75);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DPGPRmsError = zeros(6, max_length);
DPGPRmsError(6,:) = 0:0.5:(max_length-1)*0.5;
for i = 1:max_length
   %errorVec = sort(rmsErrorTmp(:,i));
   errorVec = sort(DPGPRmsErrorRaw(:,i));
   DPGPRmsError(1,i) = mean(errorVec(1:floor(count(i))));
   DPGPRmsError(2,i) = std(errorVec(1:floor(count(i))));
   DPGPRmsError(3,i) = prctile(errorVec(1:floor(count(i))),50);
   DPGPRmsError(4,i) = prctile(errorVec(1:floor(count(i))),25);
   DPGPRmsError(5,i) = prctile(errorVec(1:floor(count(i))),75);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
diffRmsErrorRaw =  (exptRmsErrorRaw - DPGPRmsErrorRaw);
diffRmsErrorRaw = diffRmsErrorRaw(:, 1:max_length);
diffRmsError = zeros(6, max_length);
diffRmsError(6,:) = 0:0.5:(max_length-1)*0.5;
for i = 1:max_length
   %errorVec = sort(rmsErrorTmp(:,i));
   errorVec = sort(diffRmsErrorRaw(:,i));
   diffRmsError(1,i) = mean(errorVec(1:floor(count(i))));
   diffRmsError(2,i) = std(errorVec(1:floor(count(i))));
   diffRmsError(3,i) = prctile(errorVec(1:floor(count(i))),50);
   diffRmsError(4,i) = prctile(errorVec(1:floor(count(i))),25);
   diffRmsError(5,i) = prctile(errorVec(1:floor(count(i))),75);
end

diffRmsErrorRaw_outliers = diffRmsErrorRaw;
diffRmsError_outliers = zeros(3, max_length);
diffRmsError_outliers(3,:) = 0:0.5:(max_length-1)*0.5;
for i = 1:max_length
   %errorVec = sort(rmsErrorTmp(:,i));
   count_sub = 0;
   for j = 1:num_traj
      if (diffRmsErrorRaw_outliers(j,i) > -0.1 && ...
         diffRmsErrorRaw_outliers(j,i) < 0.1) 
         count_sub = count_sub +1;
         diffRmsErrorRaw_outliers(j,i) = NaN;
      end
   errorVec = sort(diffRmsErrorRaw_outliers(:,i));
   diffRmsError_outliers(1,i) = mean(errorVec(1:floor(count(i)-count_sub)));
   diffRmsError_outliers(2,i) = std(errorVec(1:floor(count(i)-count_sub))); 
   end
end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% figure
% line1 = shadedErrorBar(exptRmsError(3,:),...
%     exptRmsError(1,:), exptRmsError(2,:),'b',1);
% hold on
% line2 = shadedErrorBar(DPGPRmsError(3,:), ...
%     DPGPRmsError(1,:), DPGPRmsError(2,:),'r',1);
% %line3 = shadedErrorBar(singleRmsError.rmsErrors(3,:), ...
% %    singleRmsError.rmsErrors(1,:), singleRmsError.rmsErrors(2,:),'g',1);
% 
% title('RMS prediction error for the next 3-5s');
% xlabel('time since started tracking a pedestrian (s)');
% ylabel('RMS error (m)');
% %legend([line1.mainLine, line2.mainLine, line3.mainLine], 'new method', 'DPGP', 'all singleton clusters');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
x = exptRmsError(6,:);
line1 = plot(x, exptRmsError(3,:),'b');
line2 = plot(x, DPGPRmsError(3,:),'r');
shade1 = fill([x, fliplr(x)], [exptRmsError(4,:), fliplr(exptRmsError(5,:))], 'b',...
    'facealpha', 0.2, 'edgecolor', 'w');
shade2 = fill([x, fliplr(x)], [DPGPRmsError(4,:), fliplr(DPGPRmsError(5,:))], 'r',...
    'facealpha', 0.2, 'edgecolor', 'w');
legend([line1, line2, shade1], 'DPAR', 'DPGP', '25-75 prctl');
xlabel('time since started tracking a pedestrian (s)');
ylabel('RMS error (m)');
title('RMS prediction error for the next 3-5s');
hold off

figure
boxplot(diffRmsErrorRaw);
set(gca,'XTickLabel',{' '})
set(gca,'XTick',[1:5:max_length])
set(gca,'XTickLabel',(1:5:max_length)/2-0.5)
%set(gca,'XTickLabel',{' '})
title('RMS prediction error diff (DPAR - DPGP) for the next 3-5s');
xlabel('time since started tracking a pedestrian (s)');
ylabel('RMS error (m)');



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% for hand labeled clustering assignment
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
diffRmsErrorRaw_label =  (exptRmsErrorRaw - exptRmsErrorRaw_label);
diffRmsErrorRaw_label = diffRmsErrorRaw_label(:, 1:max_length);
diffRmsError_label = zeros(6, max_length);
diffRmsError_label(6,:) = 0:0.5:(max_length-1)*0.5;
for i = 1:max_length
   %errorVec = sort(rmsErrorTmp(:,i));
   errorVec = sort(diffRmsErrorRaw(:,i));
   diffRmsError_label(1,i) = mean(errorVec(1:floor(count(i))));
   diffRmsError_label(2,i) = std(errorVec(1:floor(count(i))));
   diffRmsError_label(3,i) = prctile(errorVec(1:floor(count(i))),50);
   diffRmsError_label(4,i) = prctile(errorVec(1:floor(count(i))),25);
   diffRmsError_label(5,i) = prctile(errorVec(1:floor(count(i))),75);
end


figure
hold on
x = exptRmsError(6,:);
line1 = plot(x, exptRmsError(3,:),'b');
line2 = plot(x, exptRmsError_label(3,:),'g');

shade1 = fill([x, fliplr(x)], [exptRmsError(4,:), fliplr(exptRmsError(5,:))], 'b',...
    'facealpha', 0.2, 'edgecolor', 'w');
shade2 = fill([x, fliplr(x)], [exptRmsError_label(4,:), fliplr(exptRmsError_label(5,:))], 'g',...
    'facealpha', 0.2, 'edgecolor', 'w');
legend([line1, line2], 'DPAR', 'hand labeled');
xlabel('time since started tracking a pedestrian (s)');
ylabel('RMS error (m)');
title('RMS prediction error for the next 3-5s');
hold off

figure
boxplot(diffRmsErrorRaw_label);
set(gca,'XTickLabel',{' '})
set(gca,'XTick',[1:5:max_length])
set(gca,'XTickLabel',(1:5:max_length)/2-0.5)
%set(gca,'XTickLabel',{' '})
title('RMS prediction error diff (DPAR - true label) for the next 3-5s');
xlabel('time since started tracking a pedestrian (s)');
ylabel('RMS error (m)');

% figure
% hold on
% x = diffRmsError_outliers(3,:);
% line1 = plot(x, diffRmsError_outliers(1,:),'b');
% shade1 = fill([x, fliplr(x)], [diffRmsError_outliers(1,:) + diffRmsError_outliers(2,:), ...
%     fliplr(diffRmsError_outliers(1,:) - diffRmsError_outliers(2,:))], 'b',...
%     'facealpha', 0.2, 'edgecolor', 'w');
% xlabel('time since started tracking a pedestrian (s)');
% ylabel('RMS error (m)');
% title('RMS prediction error diff (new method - DPGP) for the next 3-5s');
% hold off



