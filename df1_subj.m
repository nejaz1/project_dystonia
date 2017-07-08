function ANA=df1_subj(subjname,fig,block,trial,varargin)
% function ANA=Individjhu2b_subj(subjname,fig,block,trial,varargin)
% Analysis of a subject/session from the SMARTS Individuation test 
%   USAGE: 
%       Individuation2b_subj('JHP0001_W0', 1);      Shows all trials of JHP0001_W0.
%       Individuation2b_subj('JHP0001_W0', 1,3);    Shows all trials of block 3
%       Individuation2b_subj('JHP0001_W0', 1,3,1);  Shows trial 1 of block 3
%       Individuation2b_subj('JHP0001_W0', 0);      generates the _ana.mat file (analyzed output) 
% 
% When you run it, CD to the analyze directory! 
if nargin<2
    fig=0;
end;

datafilename=['../data/' subjname '/IN2b_' subjname '.mat'];
outfilename=['IN2b_' subjname '_ana.mat'];
mvcfile=['../data/' subjname '/IN2mvc_' subjname '.mat'];

ANA=[];
D=load(datafilename);

% Add the mvc into the data structure 
MVC=load(mvcfile); 
D.mvc=repmat(MVC.mvcf',size(D.BN,1),1); 


if (nargin<3) 
    s=1;
    trials=[s:length(D.BN)];
else
    if (nargin<4)
        s=find(D.BN==block & D.TN==1); 
        trials=[s:length(D.BN)];
    else
        s=find((D.BN==block) & (D.TN==trial));
        trials=[s];
    end;
end;
oldblock=-1;

for i=trials
   if (oldblock~=D.BN(i)) 
        oldblock=D.BN(i);
        load(['../data/' subjname '/IN2b_' subjname '_' num2str(D.BN(i),'%2.2d') '.mat']); %data label
   end;
   if fig
       fprintf('Block %d, Trial %d\n',D.BN(i),D.TN(i));
   end
   C=df1_trial(MOV{D.TN(i)},getrow(D,i),fig);
   ANA=addstruct(ANA,C,'row','force'); 
end;

D=ANA; 
save(outfilename,'-struct','D');  %added struct to make data label changeable via load (used as a fcn)
