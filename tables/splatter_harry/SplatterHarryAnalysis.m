cd('tables/splatter_harry');
%%
LMS_LMinusMContrastPre = csvread('/Users/spitschan/Documents/MATLAB/Projects/Spitschan201x_MaxMel/tables/splatter_harry/preValidationStats_LMSStimulation_LMinusMContrast.txt',1);
LMS_LMinusMContrastPost = csvread('/Users/spitschan/Documents/MATLAB/Projects/Spitschan201x_MaxMel/tables/splatter_harry/postValidationStats_LMSStimulation_LMinusMContrast.txt',1);

[~, s] = sort(LMS_LMinusMContrastPre(:, 1));
LMS_LMinusMContrastPre = LMS_LMinusMContrastPre(s, 3:7);

[~, s] = sort(LMS_LMinusMContrastPost(:, 2));
LMS_LMinusMContrastPost = LMS_LMinusMContrastPost(s, 3:7);

LMS_LMinusMContrast = [LMS_LMinusMContrastPre LMS_LMinusMContrastPost]
mean(LMS_LMinusMContrast, 2);