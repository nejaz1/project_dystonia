function ANA=Individjhu2b_subj_dystonia(subjname,fig,block,trial,varargin)

% Analysis of a subject/session from the DYSTONIA Individuation test 

%   Programming: Ejaz 2015
%   Bugs: Sadnicka 2015
%   USAGE: 
%       Individuation2b_subj('JHP_0001_W0', 1);      Shows all trials of JHP_0001_W0.
%       Individuation2b_subj('JHP_0001_W0', 1,3);    Shows all trials of block 3
%       Individuation2b_subj('JHP_0001_W0', 1,3,1);  Shows trial 1 of block 3
%       Individuation2b_subj('JHP_0001_W0', 0);      generates the _ana.mat file (analyzed output) 
% 
 
if nargin<2
    fig=0;
end;

datafilename = ['IN2b_' subjname '.mat'];
outfilename  = ['IN2b_' subjname '_ana.mat'];
mvcfile      = ['IN2mvc_' subjname '.mat'];
mvc_MOVfile  = ['IN2mvc_' subjname '_00.mat'];
anaDir       = '/Volumes/External/data/FingerPattern_dystonia/Individuation_EMG/analysis';

ANA = [];
D   = load(datafilename);

% Add the mvc into the data structure 
MVC     = load(mvcfile); 
MVC_MOV = load(mvc_MOVfile);
for t=1:length(MVC.TN)
    dg = (MVC.hand(t)-1).*5+MVC.digit(t);
    mvcd(t) = dg;
    mvcf(t) = prctile(MVC_MOV.MOV{t}(:,dg+3),95);
end
mvc  = MVC.mvcf';
mvcf = pivottable(mvcd',[],mvcf','mean')';

D.mvcf = repmat(mvcf,[size(D.BN,1) 1]);
D.mvc  = repmat(mvc,[size(D.BN,1) 1]);

% processing individuation trials
if (nargin<5) 
    s = 1;
    trials = [s:length(D.BN)];
else
    if (nargin<6)
        s = find(D.BN==block & D.TN==1); 
        trials = [s:length(D.BN)];
    else
        s = find((D.BN==block) & (D.TN==trial));
        trials = [s];
    end;
end;

oldblock = -1;

for i=trials
   if (oldblock~=D.BN(i)) 
        oldblock = D.BN(i);
        load(['IN2b_' subjname '_' num2str(D.BN(i),'%2.2d') '.mat']); %data label
   end;
   fprintf('Block %d, Trial %d\n',D.BN(i),D.TN(i));
   
   C   = Individjhu2b_trial_dystonia(MOV{D.TN(i)},getrow(D,i),fig);
   ANA = addstruct(ANA,C,'row','force'); 
end;
D = ANA; 
save(fullfile(anaDir,outfilename),'-struct','D');  %added struct to make data label changeable via load (used as a fcn)
