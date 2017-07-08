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

%% Dystonic vs non-dystonic (CoG Distances)
% |*S1 - sensory condition*|
clear all;
d = '/Volumes/External/data/FingerPattern_dystonia/RegionOfInterest';
D = load(fullfile(d,'spatial_distances.mat'));
S = getrow(D,D.region==1 & D.metric==3 & D.condition==2);
MANOVA1(S.group,S.dist);

%%
% |*M1 - sensory condition*|
S = getrow(D,D.region==2 & D.metric==3 & D.condition==2);
MANOVA1(S.group,S.dist);
