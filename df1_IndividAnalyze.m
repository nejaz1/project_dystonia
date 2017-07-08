function varargout=df1_imana_anita(what,varargin)


baseDir= fullfile('C:/Projects/FingerPattern_dystonia');

% baseDir=        fullfile('/Users','tob', 'Projects','FingerPattern_dystonia');
% baseDir=        fullfile('/media','DATA', 'Projects','FingerPattern_dystonia');
% baseDir=        fullfile('/Volumes/MacintoshHD2/fingerPattern_dystonia');
% baseDir=        '/Users/naveed/Documents/data/FingerPattern_dystonia';
% baseDir=        fullfile('~/Projects/fingerPattern_dystonia');
% baseDir=        fullfile('/Volumes/MotorControl/project/FingerPattern_dystonia);
baseDir = '/Users/joern/Projects/fingerPattern_dystonia';
% baseDir = '/Volumes/MacintoshHD2/fingerPattern_dystonia';
behaviourDir    = fullfile(baseDir, 'Behavioural_data');
behaviourEMGDir = fullfile(baseDir, 'Individuation_EMG');
groupData=      fullfile(baseDir, 'group_data');
groupDir=       fullfile(baseDir, 'group_analysis');
anatomicalDir=  fullfile(baseDir, 'anatomicals');
freesurferDir=  fullfile(baseDir, 'surfaceFreesurfer');
caretDir=       fullfile(baseDir, 'surfaceCaret');
regDir=         fullfile(baseDir, 'RegionOfInterest');
glmName= {'glm_firstlevel_1','glm_firstlevel_2'};
glmDir= fullfile(baseDir,glmName{1});


subj_MTname={   'MT02554', 'MT02613', 'MT02614', 'MT02689', 'MT02815', ...
    'MT02816', 'MT02852', 'MT02851', 'MT02855', 'MT02853', ...
    'MT02850', 'MT02859', 'MT02860'};%'MT02481','MT02481',

subj_name={     'd01', 's01', 'd02', 's02', 'd04', ...
    's04', 'd06', 's03', 'd07', 'd08', ...
    'd09', 'd10', 'd11'};

subj_group=[    2 1 2 1 2 ...
    1 2 1 2 2 ...
    2 2 2]';  % 2 for dystonic 1 for control

subj_instrum=[  1 1 1 1 1 ...
    2 2 1 2 2 ...
    1 2 2]; % 1 for piano 2 for guitar

subj_clin = [   % columns 1:10 represent LH 1(thumb):5(little finger) then RH 6 (thumb):10(littlefinger).  0 no symptoms.  1 dystonic
    0 0 0 0 0 1 0 0 0 0; %d01
    0 0 0 0 0 0 0 0 0 0; %s01
    0 0 0 0 0 1 0 0 0 0; %d02
    0 0 0 0 0 0 0 0 0 0; %s02
    0 1 0 0 0 0 0 1 1 0; %d04
    0 0 0 0 0 0 0 0 0 0; %s04
    0 0 0 0 0 1 1 0 0 0; %d06
    0 0 0 0 0 0 0 0 0 0; %s03
    0 0 0 0 0 0 1 1 1 0; %d07
    0 0 0 0 0 0 1 1 1 0; %d08
    0 0 0 0 0 0 0 0 1 1; %d09
    0 0 0 0 0 0 0 1 1 0; %d10
    0 0 0 0 0 1 0 0 0 0] ;%d11

subj_NumVol= [  150 144 144 144 144 ...
    144 144 144 144 144 ...
    144 144 144]; % number of volumes in a functional run =======05.09.2012 %[0 0 ...] for s01-1 s01-2

dummyScans= 3; % number of scans that are deleted from the functional data =======05.09.2012

subj_NumVol= subj_NumVol-dummyScans; % delete dummy scans from the total number of volumes in the functional data =======05.09.2012

fingerPairs={'12','13','14','15','23','24','25','34','35','45'};
radius= 12;
numVox= 60;
run={'01','02','03','04','05','06','07', '08'};
hem={'lh','rh'};
hemName={'LeftHem','RightHem'};
condition= {'', ''};

maskname={'mask.img', 'mask_suit.nii'};

atlasA={'i','x'};
atlasname={'fsaverage','fsaverage_sym'};

numregions_surf=8;
numregions_suit=0;
numregions=numregions_surf+numregions_suit;
regname={'S1','M1','PMd','PMv','SMA','V12','AIP','OPJ','Cant','Cinf'};

switch(what)
    case 'subj' 
        cd(fullfile(behaviourEMGDir,'analyze')); 
        for i=1:length(subj_name) 
            df1_subj(subj_name{i},0);
        end; 
        
    case 'force_scanner'
        sn=varargin{1};
        cd(fullfile(behaviourDir,subj_name{sn}))
        D=dload(['DF1_', subj_name{sn}, '.dat']);
        barplot(D.hand, D.maxForce, 'split', [D.stimType D.digit], 'subset', [D.announce==0& D.digit~=0], 'facecolor',[0.5 0.5 0.5]);
        set(gcf,'PaperPosition',[1 1 6 2]);
        wysiwyg;
        saveas(gcf, 'maxForce', 'jpg')
        
    case 'make_alldat'     % case for Anna behavioural analysis
        S=[];
        cd(fullfile(behaviourEMGDir,'analysis')); 
        for i=1:length(subj_name)
            Si = load(fullfile(sprintf('IN2b_%s_ana.mat',subj_name{i})));
            Si.SN = zeros(length(Si.BN),1)+i;
            Si.group =  zeros(length(Si.BN),1)+subj_group(i);
            S=addstruct(S,Si);
        end;
        save(fullfile(behaviourEMGDir,'analysis','IN2b_alldat.mat'),'-struct','S');
        
        
    case 'make_individgraph_bysubject'
        
        % for this case all subjects are included even (3 currently do not
        % have accompanying imaging data)
        CAT.linecolor = {'r', 'b', 'g', 'k', 'm'};
        CAT.markercolor = {'r', 'b', 'g', 'k', 'm'};
        CAT.markerfill = {'r', 'b', 'g', 'k', 'm'};
        
        % 11 dystonia subjects
        D = []
        cd(fullfile(baseDir,'Individuation_EMG','analysis'));
        D = load('IN2b_alldat.mat');
        for i=unique(D.SN)'
            d = getrow(D, D.SN==i)
            
            h = figure
            subplot(1,2,1),
            xyplot(d.actForce, d.meanDevP, d.targetForce, 'split', d.digit, 'subset', d.hand==1, 'CAT', CAT, 'leg', 'auto');
            title('LH'), axis([0 40 0 2.5]), xlabel('actForce'), ylabel('meanDevP')
            subplot(1,2,2)
            xyplot(d.actForce, d.meanDevP, d.targetForce, 'split', d.digit, 'subset', d.hand==2, 'CAT', CAT, 'leg', 'auto');
            title('RH'), axis([0 40 0 2.5]), xlabel('actForce'), ylabel('meanDevP')
            % saveas (h, sprintf('d%d',i), 'bmp')
        end
        
        % for 5 control musicians
       
        
%     case 'make_summDat'
%         
%         % 11 dystonia subjects
%         D = [];
%         cd C:\Projects\FingerPattern_dystonia\Individuation_EMG\analysis
%         D = load('alldat_dys.mat');
%         for i=1:11;
%             sprintf('d%d',i) = getrow(D, D.SN==i);
%             for j = 1:5;
%                 tmp = getrow(sprintf('d%d',i), (sprintf('d%d.hand',i))==1);
%                 sprintf('L%d_%d',j,i) = getrow(tmp, (d.digit == j));
%             end
%             
%             for j= 1:5;
%                 tmp = getrow(sprintf('d%d',i), (sprintf('d%d.hand',i))==2);
%                 sprintf('R%d_%d',j,i) = getrow(tmp, (d.digit == j));
%             end
%         end
%         
%         % Creates new matrix with clinical data for this case in which columns 1:10 represent
%         % LH 1(thumb):5(little finger) then RH 6 (thumb):10(littlefinger).
%         % 0 no symptoms.  1 dystonic
%         
%         dys_clin = [   %
%             0 0 0 0 0 1 0 0 0 0; %d01
%             0 0 0 0 0 1 0 0 0 0; %d02
%             0 0 0 0 0 0 1 0 0 0; %d03
%             0 1 0 0 0 0 0 1 1 0; %d04
%             0 0 0 0 0 0 0 1 1 1; %d05
%             0 0 0 0 0 1 1 0 0 0; %d06
%             0 0 0 0 0 0 1 1 1 0; %d07
%             0 0 0 0 0 0 1 1 1 0; %d08
%             0 0 0 0 0 0 0 0 1 1; %d09
%             0 0 0 0 0 0 0 1 1 0; %d10
%             0 0 0 0 0 1 0 0 0 0] %d11
%         
%         % sprintf('summDat.d%d', i) = dys_clin(i,1:10)
%         
%         
    otherwise
        error('no such case!')
end

function  data = df1_acc(LDAinput,d,h,s,r)
% function  out = sl1_calAcc(LDAinput, varargin)
% LDAinput is a P x N input matrix (P voxels, N trials)
LDAinput= LDAinput(~isnan(sum(LDAinput, 2)), :);

for hand=[0 1] % hand 0=> left 1=> right
    for stimType= [0 1] % stimType  0=> of 1=> passive mov
        if ~(hand==0 && stimType==1)
            indx=find(h==hand & s==stimType);
            X=LDAinput(:,indx);
            [U,S,V]=svd(X,0);
            for i=1:max(r)
                testI{i}=find(r(indx)==i);
                trainI{i}=find(r(indx)~=i);
            end;
            %----do the classification (p*n matrix as input)
            data(hand+1+stimType,1)=crossval_takemultipleout(@classify_lda_KclassesQuicker,S*V',d(indx),trainI,testI);
            %hand+1+stimType %Check
        end
    end
end

function COG=df1_COG(data,coord)
i=find(~isnan(data));
data=data(i,:);coord=coord(i,:);
low=0; % prctile(data,10);
d=data-low;
d(d<0)=0;
COG(:,1)=sum(coord(:,1).*d)./sum(d);
COG(:,2)=sum(coord(:,2).*d)./sum(d);
COG(:,3)=sum(coord(:,3).*d)./sum(d);

function [C u res] = df1_decomp_struct(y,c,h,b) %var=
%DATA needs to be arranged as first 8 sensory blocks followed by 8 motor blocks
Kb=10;K=5;B=8;
a=0.95;

[P,N]=size(y);
tt=c+(h)*5;
y=y(~isnan(y(:,1)),:);
[y,Sw]=mva_prewhiten(y,tt,'diagonal',1);

Zc=[kron(eye(2),ones(5,1)), [eye(10)]];
Z=repmat(Zc,8,1);
Z=[Z kron(eye(16),ones(5,1))];

F=repmat(eye(K),B,1);
E=kron(eye(B),ones(K,1));
Z=[[ones(K*B,1) zeros(K*B,1) F zeros(K*B,K) E zeros(K*B,B)];...
    [zeros(K*B,1) ones(K*B,1) zeros(K*B,K) F zeros(K*B,B) E]];
Ac={};
Ac{1}=blockdiag([1 0;0 0],zeros(Kb,Kb),zeros(B*2,B*2));                                         % Var_a_condition sensory
Ac{2}=blockdiag([0 0;0 1],zeros(Kb,Kb),zeros(B*2,B*2));                                         % Var_a_condition motor
Ac{3}=blockdiag(zeros(2,2),eye(K),zeros(K,K),zeros(B*2,B*2));                                   % var_b_finger sensory
Ac{4}=blockdiag(zeros(2,2),zeros(K,K),eye(K),zeros(B*2,B*2));                                   % var_b_finger motor
Ac{5}=Ac{1}+Ac{2}+a*blockdiag((1-eye(2)),zeros(Kb,Kb),zeros(B*2,B*2));                          % Cov_a condition
Ac{6}=blockdiag(zeros(2,2),eye(Kb)+a*(diag(ones(K,1),K)+diag(ones(K,1),-K)),zeros(B*2,B*2));    % Cov_b_finger
Ac{7}=blockdiag(zeros(2,2),zeros(Kb,Kb),eye(B),zeros(B));                                       % Var_Block: sensory
Ac{8}=blockdiag(zeros(2,2),zeros(Kb,Kb),zeros(B),eye(B));                                       % Var_Block: motor
Ac{9}=Ac{7}+Ac{8}+a*blockdiag(zeros(2,2),zeros(Kb,Kb),diag(ones(B,1),B)+diag(ones(B,1),-B));    % Cov_s

[G,h,u,l,n,jumpI]=mvpattern_covcomp(y',Z,'Ac',Ac,'meanS',1,'num_iter',1000,'Ac',Ac,'TolL',1e-4,'accel_method','Aitken');
C.Dvar_cs= G(1,1);
C.Dvar_cm= G(2,2);
C.Dcov_c= G(1,2);
C.Dvar_fs= G(3,3);
C.Dvar_fm= G(8,8);
C.Dcov_f= G(3,8);
C.Dvar_bs= G(13,13);
C.Dvar_bm= G(21,21);
C.Dcov_b=G(13,21);
C.Dvar_e= h(end, end);

P=size(y',2);
b=pinv(Z*P)*sum(y',2);
r=y'-repmat(Z*b,1,P);
res=r-Z*u;

function err=exp_error(sigma,x,r); % Error under the exponential model(FWHM)
rp=exp(-x.^2/(2*sigma.^2));
err=sum((r-rp).^2);

function rp=exp_pred(sigma,x); % Error under the exponential model(FWHM)
rp=exp(-x.^2/(2*sigma.^2));

