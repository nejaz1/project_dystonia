%% Bold activity in contralateral hemisphere
% |*S1 - left hand*|
clear all;
d = '/Volumes/External/data/FingerPattern_dystonia/RegionOfInterest';
D = load(fullfile(d,'mean_betas.mat'));
x = getrow(D,D.stimType==0 & D.regType==1 & D.hand~=D.regSide & D.hand==0);
anovaMixed(x.meanAct,x.sn,'within',x.digit,{'digit'},'between',x.group,{'group'});

%%
% |*S1 - right hand*|
x = getrow(D,D.stimType==0 & D.regType==1 & D.hand~=D.regSide & D.hand==1);
anovaMixed(x.meanAct,x.sn,'within',x.digit,{'digit'},'between',x.group,{'group'});

%%
% |*M1 - left hand*|
x = getrow(D,D.stimType==0 & D.regType==2 & D.hand~=D.regSide & D.hand==0);
anovaMixed(x.meanAct,x.sn,'within',x.digit,{'digit'},'between',x.group,{'group'});

%%
% |*M1 - right hand*|
x = getrow(D,D.stimType==0 & D.regType==2 & D.hand~=D.regSide & D.hand==1);
anovaMixed(x.meanAct,x.sn,'within',x.digit,{'digit'},'between',x.group,{'group'});

%% Dystonic vs non-dystonic (aervaeCoG Distances)
% |*S1 M1 - sensory condition*|
D = load('spatial_distances.mat');
D.mDist = mean(D.dist,2); 

fprintf('S1 average distance:\n'); 
ttestDirect(D.mDist,D.group,2,'independent','subset',D.region==1 & D.metric==6 & D.condition==2); 
fprintf('M1 average distance:\n'); 
ttestDirect(D.mDist,D.group,2,'independent','subset',D.region==2 & D.metric==6 & D.condition==2); 

%% Dystonic vs non-dystonic (Distances patter)
% |*S1 M1 - sensory condition*|
D = load('spatial_distances.mat');
S1= getrow(D,D.region==1 & D.condition==2 & D.metric==6 & D.group==1); 
S2= getrow(D,D.region==1 & D.condition==2 & D.metric==6 & D.group==2);
anova_rm({S1.dist,S2.dist});

% Note that stimtype == 1 is the correct choice for sensory. 
fprintf('M1 pattern:\n'); 
S1= getrow(D,D.region==2 & D.condition==2 & D.metric==6 & D.group==1); 
S2= getrow(D,D.region==2 & D.condition==2 & D.metric==6 & D.group==2);
anova_rm({S1.dist,S2.dist});

%% Dystonic vs non-dystonic (RSA distances)
% average distances 
D = load('reg_distance_raw.mat');
D.mDist = mean(D.dist,2); 

% Note that stimtype == 1 is the correct choice for sensory. 
fprintf('S1 average distance:\n'); 
ttestDirect(D.mDist,D.group,2,'independent','subset',D.region==1 & D.stimtype==1); 
fprintf('M1 average distance:\n'); 
ttestDirect(D.mDist,D.group,2,'independent','subset',D.region==2 & D.stimtype==1); 

%% %% Dystonic vs non-dystonic (RSA distances)
% average distances 
D = load('reg_distance_raw.mat');
% Normalization would be this, however, it's more straightforward looking
% just at the interaction. 
% D.dist = bsxfun(@rdivide,D.dist,sqrt(mean(D.dist.^2,2)))
fprintf('S1 pattern:\n'); 
S1= getrow(D,D.region==1 & D.stimtype==1 & D.group==1); 
S2= getrow(D,D.region==1 & D.stimtype==1 & D.group==2); 
anova_rm({S1.dist,S2.dist});

% Note that stimtype == 1 is the correct choice for sensory. 
fprintf('M1 pattern:\n'); 
S1= getrow(D,D.region==2 & D.stimtype==1 & D.group==1); 
S2= getrow(D,D.region==2 & D.stimtype==1 & D.group==2); 
anova_rm({S1.dist,S2.dist});

