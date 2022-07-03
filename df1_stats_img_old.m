%% Statistics on Imaging measures reported in Sadnika et al. (2022)
%% Intact finger representation within primary sensorimotor cortex of Musician?s Dystonia
% Find the correct directory for distance files: 
clear all;
% d = '/Volumes/diedrichsen_data$/data/FingerPattern/FingerPattern_dystonia/RegionOfInterest';
cd('./summary_data'); 


%% Mean RSA distances for control and dystonias groups
st = 1  % Stimulus type 1: passive 0: Active 
D = load('reg_distance_raw.mat');
D.mDist = sqrt(mean(D.dist,2)); 

% Loop over the regions 
region = {'S1','M1'}; 
for r=1:2
    T=getrow(D,D.region==r & D.stimtype==st & D.hand ==1); 
    fprintf('%s average distance:\n',region{r});
    for i = 1:2
        MEAN(i) = mean(T.mDist(T.group==i));
        N(i) = sum(T.group==i); 
        SEM(i)=std(T.mDist(T.group==i))/sqrt(N(i));
    end; 
    fprintf('Controls: %2.3f SEM: %2.3f\n',MEAN(1),SEM(1)); 
    fprintf('Patients: %2.3f SEM: %2.3f\n',MEAN(2),SEM(2));
    ttestDirect(T.mDist,T.group,2,'independent');   
    % Bayes Factor against a 29% reduction from controls: 
    % Standard error of difference 
    SE = sqrt(var(T.mDist(T.group==1))/N(1)+var(T.mDist(T.group==2))/N(2));
    % Expected effect: 
    expEff = (MEAN(1)-MEAN(1)*(1-0.29))/SE; 
    trueEff = (MEAN(1)-MEAN(2))/SE; 
    % Bayes factor
    BayesFact = tpdf(-trueEff,15)/tpdf(expEff-trueEff,15); 
    fprintf('bayes factor against a 29 percent reduction: %f \n\n',BayesFact);        
    
end;     
%% Dystonic vs non-dystonic patterns (RSA distances)
% average distances 
st = 1  % Stimulus type 1: passive 0: Active 
D = load('reg_distance_raw.mat');
% Normalization would be this, however, it's more straightforward looking
% just at the interaction. 
% D.dist = bsxfun(@rdivide,D.dist,sqrt(mean(D.dist.^2,2)))
fprintf('S1 pattern:\n'); 
S1= getrow(D,D.region==1 & D.stimtype==st & D.hand==1 & D.group==1); 
S2= getrow(D,D.region==1 & D.stimtype==st & D.hand==1 & D.group==2); 
anova_rm({S1.dist,S2.dist});

% Note that stimtype == 1 is the correct choice for sensory. 
fprintf('M1 pattern:\n'); 
S1= getrow(D,D.region==2 & D.stimtype==st & D.group==1); 
S2= getrow(D,D.region==2 & D.stimtype==st & D.group==2); 
anova_rm({S1.dist,S2.dist});

%% Affected vs. non-affected finger distances. 
D = load('reg_distance_raw.mat');
S= dload('subject_list.txt');
T=getrow(D,D.region==2 & D.stimtype==1 & D.group==2);
S = getrow(S,S.group==2 * S.scan==1); 
% Calculate whether a distance should be affected or not (1=unaffected)
S.DistAff=[(1-S.R1).*(1-S.R2),...
           (1-S.R1).*(1-S.R3),...
           (1-S.R1).*(1-S.R4),...
           (1-S.R1).*(1-S.R5),...
           (1-S.R2).*(1-S.R3),...
           (1-S.R2).*(1-S.R4),...
           (1-S.R2).*(1-S.R5),...
           (1-S.R3).*(1-S.R4),...
           (1-S.R3).*(1-S.R5),...
           (1-S.R4).*(1-S.R5)]; 
nDist = T.dist-mean(T.dist,1);     % Subtract out mean difference
T.affect = sum(nDist.*(S.DistAff==0),2)./sum((S.DistAff==0),2);
T.unaffect = sum(nDist.*(S.DistAff==1),2)./ sum((S.DistAff==0),2);
ttest(T.unaffect,T.affect,2,'paired')


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



