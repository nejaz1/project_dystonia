function varargout=df1_imana(what,varargin)
% baseDir=        fullfile('/Users','tob', 'Projects','FingerPattern_dystonia');
%  baseDir=        fullfile('/media','DATA', 'Projects','FingerPattern_dystonia');
% baseDir=        fullfile('/Volumes/MacintoshHD2/fingerPattern_dystonia');
% baseDir=        '/Users/naveed/Documents/data/FingerPattern_dystonia';
% baseDir=        fullfile('~/Projects/fingerPattern_dystonia');
% % baseDir=        fullfile('/Volumes/MotorControl/project/FingerPattern_dystonia);
% baseDir = '/Users/joern/Projects/fingerPattern_dystonia';
% baseDir = '/Volumes/Naveed/data/FingerPattern_dystonia';
% baseDir = 'J:\data\FingerPattern_dystonia';
%baseDir = '/Volumes/ANNA/data/FingerPattern_dystonia';
%baseDir = '/Volumes/ANNA/data/FingerPattern_dystonia';
% baseDir = '/Users/naveed/Documents/data/FingerPattern_dystonia';
% baseDir = '/Volumes/External/data/FingerPattern_dystonia';
% baseDir = '/Volumes/Naveed/data/FingerPattern_dystonia';
%baseDir = '/Volumes/MotorControl/data/FingerPattern_dystonia';
% baseDir = '/srv/diedrichsen/FingerPattern_dystonia';
baseDir = '/Volumes/diedrichsen_data$/data/FingerPattern/FingerPattern_dystonia';

behaviourDir    = fullfile(baseDir, 'Behavioural_data');
emgDir          = fullfile(baseDir, 'Individuation_EMG/data');
groupData=      fullfile(baseDir, 'group_data');
groupDir=       fullfile(baseDir, 'group_analysis');
anatomicalDir=  fullfile(baseDir, 'anatomicals');
freesurferDir=  fullfile(baseDir, 'surfaceFreesurfer');
caretDir=       fullfile(baseDir, 'surfaceCaret');
regDir=         fullfile(baseDir, 'RegionOfInterest');
glmName= {'glm_firstlevel_1','glm_firstlevel_2'};
glmDir= fullfile(baseDir,glmName{1});
analysisDir     = [baseDir '/analysis'];
figureDir       = [baseDir '/Individuation_EMG/analysis/figures'];
statsDir        = [analysisDir '/stats'];
codeDir         = '~/Matlab/Projects/project_dystonia/';

colours     = {[0,0,0],[0.6 0.6 1],[0 0 1],[0.6 1 0.6],[0 1 0]};
sty_grp     = colours([3,5]);
sty_control = colours(3);
sty_patient = colours([5]);


% All these variables should go into the patient_list.txt document
subj_MTname={   'MT02554', 'MT02613', 'MT02614', 'MT02689', 'MT02815', ...
    'MT02816', 'MT02852', 'MT02851', 'MT02855', 'MT02853', ...
    'MT02850', 'MT02859', 'MT02860','MT03687','MT03688','MT03686'};%'MT02481','MT02481',
subj_name={ 'd01', 's01', 'd02', 's02', 'd04', ...
    's04', 'd06', 's03', 'd07', 'd08', ...
    'd09', 'd10', 'd11', 's06', 's07',...
    's08','s05'};

subj_group=[2 1 2 1 2 ...
    1 2 1 2 2 ...
    2 2 2 1 1 ...
    1 1 ]';  % 2 for dystonic 1 for pianists
subj_NumVol= [  150 144 144 144 144 ...
    144 144 144 144 144 ...
    144 144 144 146 146 ...
    146 141];

emg_subj={'d01','d02','d03','d04','d05',...
    'd06','d07','d08','d09','d10',...
    'd11',...
    's01','s02','s03','s04','s05',...
    's06','s08','s09'};
emg_subjID={'021112','031212','021112','021112','071112',...
    '071112','151112','161112','031212','191112',...
    '191112',...
    '141212','101212','011112','011112','151112',...
    '240714','240714','240714'};
emg_NumRuns=[10 10 10 10 11 10 10 10 10 10,...
    10 10 10 10 10 10 10 11 10 10 10];

% Other variables
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
numregions_suit=2;
numregions=numregions_surf+numregions_suit;
regname={'S1','M1','PMd','PMv','SMA','V12','AIP','OPJ','Cant','Cinf'};
regSide=[zeros(1,10) ones(1,10)];
regType=[1:10 1:10];

use3D=1;

ens_label       = {'1/2','1/3','1/4','1/5',...
    '2/1','2/3','2/4','2/5',...
    '3/1','3/2','3/4','3/5',...
    '4/1','4/2','4/3','4/5',...
    '5/1','5/2','5/3','5/4'};


% set default plotting style
style.file(fullfile(codeDir,'df_style.m'));
style.use('default');

% open science container file
% ost.load(fullfile(baseDir,'dystonia.ost'));

switch(what)
    case 'preprocess'
        %df1_imana('preprocess', 2)
        sn=varargin{1};
        for s=sn
            df1_imana('make_nii', s)
            %                         df1_imana('setAC',s)
            %                         df1_imana('segmentation',s);
            df1_imana('slice_timing',s);
            %                         df1_imana('makefieldmap',s);
            %                         df1_imana('make_realign_unwarp',s);
            %                         df1_imana('move_images',s);
            %                         df1_imana('meanimage_bias_correction',s);
            %             figure(s); df1_imana('plot_movementparameters',s);
            %                         df1_imana('coreg', s);
            %             disp('check the coreg1');keyboard;
            %                         df1_imana('coreg_2', s);
            %             disp('check the coreg2');keyboard;
            %                         df1_imana('coreg_3', s);
            %             disp('check the coreg3');keyboard;
            %                         df1_imana('make_samealign',s);
            %                         df1_imana('check_samealign',s);
            %                         disp('check same alignment');keyboard;
            %                         df1_imana('make_maskImage',s)
            %                         df1_imana('make_glm',s)
            %                         df1_imana('estimate_glm',s, 1)
            %                         df1_imana('contrast',s,1)
            %                         df1_imana('SUIT_isolate',s)
            %                         df1_imana('SUIT_normalize',s)
            %                         df1_imana('make_mask_suit',1)
            %                         df1_imana('MVA_search',s, 60)
            %             df1_imana('MVA_do',s, 60)
            %             df1_imana('MVA_zValue',s,1)
            %             df1_imana('MNI_normalization_write',s, 1)
            %                           df1_imana('surf_freesurfer', s)
            %                           df1_imana('surf_xhemireg',s)
            %             df1_imana('surf_map_ico',s,1)
            %             df1_imana('surf_make_caret',s,1)
            %             df1_imana('surf_MVA_search',s, 1, 1); %LEFT
            %             df1_imana('surf_MVA_search',s, 2, 1); %RIGHT
            %             df1_imana('surf_MVA_do',s);
        end;
        % df1_imana('MNI_mask_anatomical', sn)
        % df1_imana('MNI_smooth',sn)
        % df1_imana('MNI_make_rfx',sn)
        % df1_imana('MNI_rfx_spm', 1:5)
    case 'make_nii' %________________STEP 1____________________________
        % df1_imana('make_nii', 1)
        % delete the first 3 images
        sn=varargin{1};
        %----go to data folder
        cd(fullfile(baseDir, 'imaging_data_raw',subj_name{sn}));
        %----loop over runs
        for i=1:8
            if strcmp(subj_name{sn}, 's01-2')
                spmj_tar2nii ([subj_MTname{sn} '.r', num2str(i) '.tar'], [subj_name{sn} '_run' run{i} '.nii']);
            else
                spmj_tar2nii ([subj_MTname{sn} '.r', num2str(i) '.tar'], [subj_name{sn} '_run' run{i} '.nii'], 'startTR', dummyScans+1,'use3D',1)
            end
        end
        spmj_tar2nii([subj_MTname{sn} '.f1.tar'], [subj_name{sn} '_magnitude.nii']);
        spmj_tar2nii([subj_MTname{sn} '.f2.tar'], [subj_name{sn} '_phase.nii']);
        spmj_tar2nii([subj_MTname{sn} '.epi.tar'], [subj_name{sn} '_epi.nii']);
        %----fieldmap
        %---- make directory
        mkdir(fullfile(baseDir, 'fieldmaps', subj_name{sn}));
        source = fullfile(baseDir, 'imaging_data_raw', subj_name{sn}, [subj_name{sn} '_magnitude.nii']);
        dest = fullfile(baseDir, 'fieldmaps', subj_name{sn}, [subj_name{sn} '_magnitude.nii']);
        movefile(source,dest);
        source = fullfile(baseDir, 'imaging_data_raw', subj_name{sn}, [subj_name{sn} '_phase.nii']);
        dest = fullfile(baseDir, 'fieldmaps', subj_name{sn}, [subj_name{sn} '_phase.nii']);
        movefile(source,dest);
        %----anatomicals
        spmj_tar2nii([subj_MTname{sn} '.ana.tar'], [subj_name{sn} '_anatomical.nii']);
        source = fullfile(baseDir, 'imaging_data_raw',subj_name{sn}, [subj_name{sn} '_anatomical.nii']);
        dest = fullfile(baseDir, 'anatomicals',subj_name{sn}, [subj_name{sn} '_anatomical.nii']);
        %movefile(source,dest);
        %----set rotations to zero
        V=spm_vol(source);
        R= spm_read_vols(V);
        V.mat= round(V.mat);
        %---- make directory
        mkdir(fullfile(baseDir, 'anatomicals', subj_name{sn}));
        %---LPI resclice
        spmj_reslice_LPI(source,'name', dest);
        %---set tranlation to zero
        V= spm_vol(dest);
        R= spm_read_vols(V);
        V.mat(1:3,4)= [0 0 0];
        spm_write_vol(V,R);
        %         %----mean epi    ==========04.09.2012
        images= {[subj_name{sn} '_epi.nii,2'] [subj_name{sn} '_epi.nii,3']};
        spmj_imcalc_mtx(images,[subj_name{sn} '_ep.nii'],'mean(X)');
        movefile([subj_name{sn} '_ep.nii'], [subj_name{sn} '_epi.nii']);
    case 'setAC' %___________________STEP 1.0__________________________
        %df1_imana('setAC',1)
        sn=varargin{1};
        ana_name= fullfile(baseDir, 'anatomicals',subj_name{sn}, [subj_name{sn} '_anatomical.nii']);
        %----coords of the AC
        %[ 0 0 0 ;     %s01-1 look it up 0    0    0 ;      %s01-2 look it up
        AC_coord=[-91 -132 -142;      %d01
            -92 -126 -152;      %s01
            -90 -130 -158;      %d02
            -90 -135 -161;      %s02
            -87 -131 -148;      %d04
            -87 -134 -155;      %s04
            -91 -137 -146;      %d06
            -88 -128 -153;      %s03
            -87 -131 -155;      %d07
            -91 -132 -154;      %d08
            -88 -123 -147;      %d09
            -91 -134 -150;      %d10
            -95 -140 -164;      %d11
            -87 -130 -168;      %s06
            -93 -128 -168;      %s07
            -88 -132 -163;      %s08
            -86 -136 -162];     %s05
        %----set the AC
        V= spm_vol(ana_name);
        R= spm_read_vols(V);
        V.mat(1:3,4)= AC_coord(sn,:);
        V.mat
        spm_write_vol(V,R);
    case 'segmentation' %____________STEP 1.1__________________________
        % df1_imana('segmentation',1)
        sn=varargin{1};
        for s=sn
            spmj_segmentation(fullfile(anatomicalDir, subj_name{s}, [subj_name{s}, '_anatomical.nii']))
        end;
    case 'slice_timing' %____________STEP 2____________________________
        %df1_imana('slice_timing',1)
        % [2:2:32 1:2:32] interleave
        prefix='';
        sn=varargin{1};
        for r= 1:8
            for i=1:(subj_NumVol(sn))
                if use3D
                    N{i} = [fullfile(baseDir, 'imaging_data_raw',subj_name{sn}, [prefix subj_name{sn},'_run',run{r},'_',num2str(i),'.nii,'])];
                else
                    N{i} = [fullfile(baseDir, 'imaging_data_raw',subj_name{sn}, [prefix subj_name{sn},'_run',run{r},'.nii,',num2str(i)])];
                end;
            end;
            J.scans{r} = N;
        end
        J.nslices = 32;
        J.tr = 2.72;
        J.ta = 2.635;
        J.so = [2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31];
        J.refslice = 1;
        J.prefix = 'a';
        matlabbatch{1}.spm.temporal.st= J;
        spm_jobman('run',matlabbatch)
    case 'makefieldmap' %____________STEP 3____________________________
        %df1_imana('makefieldmap', 1)
        prefix='a';
        %prefix='';
        sn=varargin{1};
        spmj_makefieldmap(baseDir, subj_name{sn}, run,'prefix',prefix, 'subfolderRawdata', '','use3D',use3D);
    case 'make_realign_unwarp' %_____STEP 4____________________________
        % df1_imana('make_realign_unwarp',2)
        %prefix='ab'; %with bias corred epi's
        prefix='a'; %with slice time corrected images
        %prefix=''; % on raw images
        sn=varargin{1}
        spmj_realign_unwarp(baseDir, subj_name{sn}, run, 1, subj_NumVol(sn),'prefix',prefix, 'subfolderRawdata','') %76
        % maybe produce mean.nii
        % spm_imcalc_ui({},'mean.nii','mean(X)',{1,[],[],[]})
    case 'move_images' %_____________STEP 5____________________________
        % df1_imana('move_images',2)
        prefix='ua'; %'ua'
        sn=varargin{1};
        mkdir(fullfile(baseDir, 'imaging_data',subj_name{sn}))
        for j=1:numel(run);
            source = fullfile(baseDir, 'imaging_data_raw',subj_name{sn}, [prefix,subj_name{sn},'_run',run{j},'.nii']);
            dest = fullfile(baseDir, 'imaging_data',subj_name{sn}, [prefix,subj_name{sn},'_run',run{j},'.nii']);
            movefile(source,dest);
            source = fullfile(baseDir, 'imaging_data_raw', subj_name{sn},['rp_a' subj_name{sn},'_run',run{j},'.txt']); %rp_ab
            dest = fullfile(baseDir, 'imaging_data',subj_name{sn}, ['rp_a' subj_name{sn},'_run',run{j},'.txt']); %rp_ab
            movefile(source,dest);
        end;
    case 'meanimage_bias_correction' %STEP 5.1_________________________
        % df1_imana('meanimage_bias_correction',1)
        sn=varargin{1};
        prefix='ua';
        %prefix='u';
        %----full epi mean image
        P{1}={fullfile(baseDir, 'imaging_data_raw',subj_name{sn}, [subj_name{sn},'_epi.nii'])};
        spm_biascorrect(P{1});
        source = fullfile(baseDir, 'imaging_data_raw',subj_name{sn},['b',subj_name{sn},'_epi.nii']);
        dest = fullfile(baseDir, 'imaging_data',subj_name{sn}, [subj_name{sn}, '_epi.nii']);
        copyfile(source,dest);
        %----functional mean image
        P{1}={fullfile(baseDir, 'imaging_data_raw',subj_name{sn},['mean' prefix subj_name{sn},'_run',run{1},'.nii'])};
        spm_biascorrect(P{1});
        source = fullfile(baseDir, 'imaging_data_raw',subj_name{sn},['bmean' prefix subj_name{sn},'_run',run{1},'.nii']);
        dest = fullfile(baseDir, 'imaging_data',subj_name{sn}, ['meanepi_' subj_name{sn} '.nii']);
        copyfile(source,dest);
    case 'plot_movementparameters' %___________________________________
        % df1_imana('plot_movementparameters',2)
        sn=varargin{1};
        X=[];
        for j=1:numel(run);
            x= dlmread (fullfile(baseDir, 'imaging_data',subj_name{sn}, ['rp_a' subj_name{sn},'_run',run{j},'.txt']));
            X= [X; x];
        end
        set(0,'DefaultAxesColorOrder',[0 0 0], 'DefaultAxesLineStyleOrder','-|.|:')
        subplot(2,1,1); plot(X(:,1:3) )
        legend('x', 'y', 'z', 'location' , 'EastOutside')
        
        subplot(2,1,2); plot(X(:,4:6)*180/pi)
        legend('pitch', 'roll', 'yaw', 'location' , 'EastOutside')
    case 'coreg' %___________________STEP 6____________________________
        % df1_imana('coreg',1)
        % coregtool;
        sn=varargin{1};
        %----full_epi => anatomicals
        J.ref = {fullfile(baseDir,'anatomicals',subj_name{sn},[subj_name{sn},'_anatomical.nii'])};
        J.source = {fullfile(baseDir, 'imaging_data',subj_name{sn}, [subj_name{sn},'_epi.nii'])};
        J.other = {''};
        J.eoptions.cost_fun = 'nmi';
        J.eoptions.sep = [4 2];
        J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        J.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estimate=J;
        spm_jobman('run',matlabbatch);
    case 'coreg_2' %_________________STEP 6.1__________________________
        % df1_imana('coreg_2',1)
        % coregtool;
        sn=varargin{1};
        %----func => full_epi
        J.ref = {fullfile(baseDir, 'imaging_data',subj_name{sn}, [subj_name{sn},'_epi.nii'])};
        J.source = {fullfile(baseDir, 'imaging_data',subj_name{sn}, ['meanepi_',subj_name{sn},'.nii'])};
        J.other = {''};
        J.eoptions.cost_fun = 'nmi';
        J.eoptions.sep = [4 2];
        J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        J.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estimate=J;
        spm_jobman('run',matlabbatch);
    case 'coreg_3' %_________________STEP 6.2__________________________
        % df1_imana('coreg_3',1)
        % coregtool;
        sn=varargin{1};
        %----func => full_epi
        J.ref = {fullfile(baseDir, 'anatomicals',subj_name{sn}, [subj_name{sn},'_anatomical.nii'])};
        J.source = {fullfile(baseDir, 'imaging_data',subj_name{sn}, ['meanepi_',subj_name{sn},'.nii'])};
        J.other = {''};
        J.eoptions.cost_fun = 'nmi';
        J.eoptions.sep = [4 2];
        J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        J.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estimate=J;
        spm_jobman('run',matlabbatch);
    case 'coreg_4' % using bias corrected image
        % df1_imana('coreg_3',1)
        % coregtool;
        sn=varargin{1};
        %----func => full_epi
        J.ref = {fullfile(baseDir, 'anatomicals',subj_name{sn}, [subj_name{sn},'_anatomical.nii'])};
        J.source = {fullfile(baseDir, 'imaging_data',subj_name{sn}, ['rbmeanepi_',subj_name{sn},'.nii'])};
        J.other = {''};
        J.eoptions.cost_fun = 'nmi';
        J.eoptions.sep = [4 2];
        J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        J.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estimate=J;
        spm_jobman('run',matlabbatch);
    case 'make_samealign' %__________STEP 7____________________________
        % df1_imana('make_samealign',1)
        %prefix='ra';
        prefix='ua';
        %prefix='u';
        sn=varargin{1};
        Q={};
        cd(fullfile(baseDir, 'imaging_data',subj_name{sn}));
        
        % subj s05 only has 7 runs
        if (sn==17)
            run=run(1:7);
        end;
        for r= 1:numel(run)
            for i=1:(subj_NumVol(sn))
                Q{end+1} = [fullfile(baseDir, 'imaging_data',subj_name{sn}, [prefix, subj_name{sn},'_run',run{r},'.nii,',num2str(i)])];
            end
        end
        if exist(fullfile(baseDir, 'imaging_data', subj_name{sn}, ['rbmeanepi_' subj_name{sn} '.nii']),'file')
            P{1} = fullfile(baseDir, 'imaging_data', subj_name{sn}, ['rbmeanepi_' subj_name{sn} '.nii']);
        else
            P{1} = fullfile(baseDir, 'imaging_data', subj_name{sn}, ['meanepi_' subj_name{sn} '.nii']);
        end;
        spmj_makesamealign_nifti(char(P),char(Q))
        %delete(char(P))
    case 'check_samealign' %_________STEP 8____________________________
        % df1_imana('check_samealign',1)
        prefix='ua';
        %prefix='u';
        sn=varargin{1};
        Q={};
        cd(fullfile(baseDir, 'imaging_data',subj_name{sn}));
        % subj s05 only has 7 runs
        if (sn==17)
            run=run(1:7);
        end;
        for r= 1:numel(run)
            for i=1:(subj_NumVol(sn))
                Q{end+1} = [fullfile(baseDir, 'imaging_data',subj_name{sn}, [prefix, subj_name{sn},'_run',run{r},'.nii,',num2str(i)])];
            end
        end
        if exist(fullfile(baseDir, 'imaging_data', subj_name{sn}, ['rbmeanepi_' subj_name{sn} '.nii']),'file')
            P{1} = fullfile(baseDir, 'imaging_data', subj_name{sn}, ['rbmeanepi_' subj_name{sn} '.nii']);
        else
            P{1} = fullfile(baseDir, 'imaging_data', subj_name{sn}, ['meanepi_' subj_name{sn} '.nii']);
        end;
        spmj_checksamealign(char(P),char(Q))
    case 'make_maskImage' %__________STEP 8.1__________________________
        % df1_imana('make_maskImage',1)
        % Do this on the parcelation maps
        sn=varargin{1};
        for s=sn
            cd(fullfile(baseDir,'imaging_data', subj_name{s}));
            %----load mean image
            if exist(fullfile(baseDir, 'imaging_data',subj_name{s}, ['rbmeanepi_' subj_name{s} '.nii']),'file')
                nam{1} = fullfile(baseDir, 'imaging_data',subj_name{s}, ['rbmeanepi_' subj_name{s} '.nii']);
            else
                nam{1}=  fullfile(baseDir, 'imaging_data',subj_name{s}, ['meanepi_' subj_name{s} '.nii']);
            end;
            %---- mask with c1,c2, & c3 from segmentation
            nam{2}= fullfile(anatomicalDir, subj_name{s},  ['c1' subj_name{s}, '_anatomical.nii'])
            nam{3}= fullfile(anatomicalDir, subj_name{s},  ['c2' subj_name{s}, '_anatomical.nii'])
            nam{4}= fullfile(anatomicalDir, subj_name{s},  ['c3' subj_name{s}, '_anatomical.nii'])
            spm_imcalc_ui(nam, 'rmask_noskull.nii', 'i1>1 & (i2+i3+i4)>0.2')
            %----reslice volume
            % V= spm_vol(nam{1});
            % spmj_reslice_vol('mask_noskull.nii', V.dim, V.mat, 'rmask_noskull.nii')
            
        end
    case 'make_glm_1' %________________STEP 9____________________________
        %df1_imana('make_glm',1)
        % set threshold for the masking
        
        subj=varargin{1};
        timePerSlice= 2.72/32; % time per slice in sec
        delay=  2.72;                   % 1 TR is seconds
        dur=    3.2*2.72;               % In seconds
        
        for sn=subj
            %----
            prefix= 'ua'; %'a' for new sequence
            numVolumes= subj_NumVol(sn); %76 for new seq
            %             num_dummys=3;
            D=dload( fullfile(behaviourDir,subj_name{sn}, ['DF1_',subj_name{sn},'.dat']));
            %D=getrow(D,D.day==7 & D.lastTrial==0); % Take scan data
            %sequences= unique(D.seqType);
            T=[];
            %correct for the number of dummy scans (not in nii file)
            %TR= D.startTR-num_dummys;
            Slice= D.startSlice-dummyScans*32;
            
            J.dir = {fullfile(baseDir,'glm_firstlevel_1', subj_name{sn})};
            if (~exist(J.dir{1},'dir'))
                mkdir(J.dir{1});
            end;
            
            J.timing.units = 'secs';% 'scans';%
            J.timing.RT = 2.72;
            J.timing.fmri_t = 16;
            J.timing.fmri_t0 = 1;
            
            % subj s05 only has 7 runs
            if ( subj==17)
                D=getrow(D,D.BN<8);
                run=run(1:7);
            end;
            
            for r=1:numel(run) %run
                for i=1:numVolumes
                    N{i} = fullfile(baseDir, 'imaging_data',subj_name{sn},[prefix,subj_name{sn},'_run',run{r},'.nii,',num2str(i)]); % ua uab 'ra'
                end;
                J.sess(r).scans= N;
                J.sess(r).cond=[];
                for h=[0 1] % hand 0=> left 1=> right
                    for s= [0 1] % stimType  0=> of 1=> passive mov
                        for d=1:5 % digit
                            if ~(h==0 && s==1)
                                idx_digit=find(D.BN==r  & D.hand==h & D.stimType==s & D.digit==d & D.announce==1);
                                J.sess(r).cond(end+1).name = ['B%dH%dC%dD%d',r,h,s,d];
                                J.sess(r).cond(end).onset = Slice(idx_digit)*timePerSlice-delay;
                                %----print information on screen
                                fprintf('run: %d \t hand: %d \t stimType: %d \t digit: %d\n',r,h,s,d);
                                fprintf('onset time in seconds: %2.2fs \t\t minutes: %2.2f \t\t slice of tgt: %2.2f \n',[([Slice(idx_digit)- delay]*timePerSlice) ([Slice(idx_digit)- delay]*timePerSlice)/60 (Slice(idx_digit)+ dummyScans*32)]');
                                %----
                                J.sess(r).cond(end).duration =  dur;
                                J.sess(r).cond(end).tmod = 0;
                                J.sess(r).cond(end).pmod = struct('name', {}, 'param', {}, 'poly', {});
                                S.sn=sn;
                                S.run=r;
                                S.numEvents=length(idx_digit);
                                S.hand=h;
                                S.regType=1;  % Digit regressor
                                S.digit=d; % Sequences laled after its type
                                S.stimType=s;         % Sequences labled 1-4
                                T=addstruct(T,S);
                            end
                        end
                    end
                end
                J.sess(r).multi = {''};
                J.sess(r).regress = struct('name', {}, 'val', {});
                J.sess(r).multi_reg = {''};
                J.sess(r).hpf = 128;
            end;
            J.fact = struct('name', {}, 'levels', {});
            J.bases.hrf.derivs = [0 0];
            J.volt = 1;
            J.global = 'None';
            J.mask = {fullfile(baseDir, 'imaging_data',subj_name{sn},'rmask_noSkull.nii')};
            J.cvi =  'wls';
            matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec=J;
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
            
            spm_jobman('run',matlabbatch);
        end;
        % varargout={J};
    case 'make_glm_2' %__________One regressor per run, but with extra regressors for error trials
        %df1_imana('make_glm',1)
        % set threshold for the masking
        
        subj=varargin{1};
        timePerSlice= 2.72/32; % time per slice in sec
        delay=  2.72;                   % 1 TR is seconds
        dur=    3.2*2.72;               % In seconds
        
        for sn=subj
            %----
            prefix= 'ua'; %'a' for new sequence
            numVolumes= subj_NumVol(sn); %76 for new seq
            %             num_dummys=3;
            D=dload( fullfile(behaviourDir,subj_name{sn}, ['DF1_',subj_name{sn},'.dat']));
            %D=getrow(D,D.day==7 & D.lastTrial==0); % Take scan data
            %sequences= unique(D.seqType);
            T=[];
            %correct for the number of dummy scans (not in nii file)
            %TR= D.startTR-num_dummys;
            Slice= D.startSlice-dummyScans*32;
            
            J.dir = {fullfile(baseDir,'glm_firstlevel_2', subj_name{sn})};
            if (~exist(J.dir{1},'dir'))
                mkdir(J.dir{1});
            end;
            
            J.timing.units = 'secs';% 'scans';%
            J.timing.RT = 2.72;
            J.timing.fmri_t = 16;
            J.timing.fmri_t0 = 1;
            
            for r=1:numel(run) %run
                for i=1:numVolumes
                    N{i} = fullfile(baseDir, 'imaging_data',subj_name{sn},[prefix,subj_name{sn},'_run',run{r},'.nii,',num2str(i)]); % ua uab 'ra'
                end;
                J.sess(r).scans= N;
                J.sess(r).cond=[];
                indx=find(D.BN==r  & D.announce==1);
                
                for i=1:length(indx) % hand 0=> left 1=> right
                    j=indx(i);
                    J.sess(r).cond(end+1).name = ['B%dH%dC%dD%d',r,D.hand(j),D.stimType(j),D.digit(j)];
                    J.sess(r).cond(end).onset = Slice(j)*timePerSlice-delay;
                    %----print information on screen
                    J.sess(r).cond(end).duration =  dur;
                    J.sess(r).cond(end).tmod = 0;
                    J.sess(r).cond(end).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    S.sn=sn;
                    S.run=r;
                    S.numEvents=1;
                    S.hand=D.hand(j);
                    S.regType=1;  % Digit regressor
                    S.digit=D.digit(j); % Sequences laled after its type
                    S.stimType=D.stimType(j);         % Sequences labled 1-4
                    T=addstruct(T,S);
                end
                J.sess(r).multi = {''};
                J.sess(r).regress = struct('name', {}, 'val', {});
                J.sess(r).multi_reg = {''};
                J.sess(r).hpf = inf;
            end;
            J.fact = struct('name', {}, 'levels', {});
            J.bases.hrf.derivs = [0 0];
            J.volt = 1;
            J.global = 'None';
            J.mask = {fullfile(baseDir, 'imaging_data',subj_name{sn},'rmask_noSkull.nii')};
            J.cvi =  'wls';
            matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec=J;
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
            
            spm_jobman('run',matlabbatch);
        end;
        varargout={J};
    case 'estimate_glm' %____________STEP 10___________________________
        % df1_imana('estimate_glm',1, 1)
        sn=varargin{1};
        glmType=varargin{2};
        for s=sn
            matlabbatch{1}.spm.tools.rwls.fmri_rwls_est.spmmat = { fullfile(baseDir,glmName{glmType},subj_name{s} ,'SPM.mat')};
            matlabbatch{1}.spm.tools.rwls.fmri_rwls_est.method.Classical = 1;
            spm_jobman('run',matlabbatch);
        end
        
        % for checking
        % load('SPM')
        % spm_rwls_resstats(SPM)
    case 'hrf_getdata' %_______________________________________________
        % to check the model quality of the glm
        %df1_imana('hrf_getdata',1)
        sn=varargin{1};
        pre=4; post=12;
        T=[];
        for s=sn
            cd(fullfile(baseDir,'glm_firstlevel_1', subj_name{s}));
            load('SPM.mat');
            SPM=spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s}));
            
            load(fullfile(regDir,sprintf('%s_regions.mat',subj_name{s})));
            R={R{2},R{10}};
            [y_raw, y_adj, y_hat, y_res,B] = region_getts(SPM,R);
            D=spmj_get_ons_struct(SPM);
            for r=1:size(y_raw,2)
                for i=1:size(D.block,1);
                    D.y_adj(i,:)=cut(y_adj(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                    D.y_hat(i,:)=cut(y_hat(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                    D.y_res(i,:)=cut(y_res(:,r),pre,round(D.ons(i)),post,'padding','nan')';
                end;
                vec=ones(size(D.event,1),1);
                D.region=vec*r;
                D.SN=vec*s;
                T=addstruct(T,D);
            end;
            fprintf('%d\n',s);
        end;
        cd(regDir);
        save ROI_M1_IRF1.mat T
    case 'hrf_plot' %__________________________________________________
        % df1_imana('hrf_plot',1)
        sn=varargin{1};
        load(fullfile(regDir,'ROI_M1_IRF1.mat'));
        T=getrow(T,T.SN==sn);
        T.type=floor((T.event-1)/5)+1;
        for r=1:2
            for t=1:3
                subplot(3,2,t*2+r-2);
                traceplot([-4:12],T.y_adj,'errorfcn','stderr','subset',T.region==r & T.type==t);
                hold on;
                traceplot([-4:12],T.y_hat,'errorfcn','stderr','subset',T.region==r & T.type==t,'linestyle',':');
                hold off;
                
            end;
        end;
        
    case 'contrast_main' %________________STEP 11
        %df1_imana('contrast',1, 1)
        sn=varargin{1};
        glmType=varargin{2};
        for subj=sn
            cd(fullfile(baseDir,glmName{glmType}, subj_name{subj}));
            load SPM;
            T=load('SPM_info.mat');
            SPM=rmfield(SPM,'xCon');
            
            con=zeros(1,size(SPM.xX.X,2));
            con(T.regType==1)=1;
            con=con/sum(con);
            SPM.xCon(1)=spm_FcUtil('Set','t_task', 'T', 'c',con',SPM.xX.xKXs);
            
            %____F contrast between seq
            for h=[0 1] % hand 0=> left 1=> right
                for s= [0 1] % stimType  0=> of 1=> passive mov
                    if ~(h==0 && s==1)
                        con=zeros(4,size(SPM.xX.X,2));
                        con(1,T.digit==1 & T.hand==h & T.stimType==s)=1;con(1,T.digit==2 & T.hand==h & T.stimType==s)=-1;
                        con(2,T.digit==2 & T.hand==h & T.stimType==s)=1;con(2,T.digit==3 & T.hand==h & T.stimType==s)=-1;
                        con(3,T.digit==3 & T.hand==h & T.stimType==s)=1;con(3,T.digit==4 & T.hand==h & T.stimType==s)=-1;
                        con(4,T.digit==4 & T.hand==h & T.stimType==s)=1;con(4,T.digit==5 & T.hand==h & T.stimType==s)=-1;
                        %%[[repmat([1:15]',8,1); [1:8]'] con'] %check
                        if (any(sum(con,2)~=0))
                            keyboard;
                        end;
                        conName= sprintf('F_hand_%i_stimType_%i', h,s);
                        SPM.xCon(h+s+2)=spm_FcUtil('Set',conName, 'F', 'c',con',SPM.xX.xKXs);
                        %h+s+2 %check
                    end
                end
                
            end;
            
            %             %_____t contrast task against rest
            for h=[0 1] % hand 0=> left 1=> right
                for s= [0 1] % stimType  0=> of 1=> passive mov
                    if ~(h==0 && s==1)
                        con=zeros(1,size(SPM.xX.X,2));
                        con(T.hand==h & T.stimType==s)=1;
                        con=con/sum(con);
                        conName= sprintf('T_hem_%i_stimType_%i', h,s);
                        SPM.xCon(h+s+5)=spm_FcUtil('Set',conName, 'T', 'c',con',SPM.xX.xKXs);
                    end
                end
            end;
            
            %             % Optional: t-contrast for each finger / hand against rest
            for h=[0 1] % hand 0=> left 1=> right
                for s= [0 1] % stimType  0=> of 1=> passive mov
                    if ~(h==0 && s==1)
                        for d=1:5
                            con=zeros(1,size(SPM.xX.X,2));
                            con(1,T.digit==d & T.hand==h & T.stimType==s)=1;
                            con=con/sum(con);
                            conName= sprintf('T_hand_%i_stimType_%i_digit_%i', h,s,d);
                            SPM.xCon(end+1)=spm_FcUtil('Set',conName, 'T', 'c',con',SPM.xX.xKXs);
                        end
                    end
                end
            end
            
            %____do the constrasts
            SPM=spm_contrasts(SPM,[1:length(SPM.xCon)]);
            save SPM SPM;
        end
    case 'copy_main_contrast'
        glmType = 1;
        for sn=unique(subj_name)
            disp(sn);
            cd(fullfile(baseDir,glmName{glmType},sn{1}));
            mkdir('contrast_main');
            
            d1=dir('spm*');
            for i=1:length(d1)
                copyfile(d1(i).name,fullfile(cd,'contrast_main',d1(i).name));
            end;
            
            d2=dir('ess_*');
            for i=1:length(d2)
                copyfile(d2(i).name,fullfile(cd,'contrast_main',d2(i).name));
            end;
            
            d3=dir('con_*');
            for i=1:length(d3)
                copyfile(d3(i).name,fullfile(cd,'contrast_main',d3(i).name));
            end;
        end;
    case 'contrast_splithalf' %________________STEP 11
        sn      = varargin{1};
        glmType = varargin{2};
        
        for subj=sn
            cd(fullfile(baseDir,glmName{glmType}, subj_name{subj}));
            disp(subj_name{subj});
            load SPM;
            T=load('SPM_info.mat');
            SPM=rmfield(SPM,'xCon');
            
            % ODD RUN: t-contrast for each finger / hand against rest
            idx = 1;
            for h=[0 1] % hand 0=> left 1=> right
                for s= [0 1] % stimType  0=> of 1=> passive mov
                    if ~(h==0 && s==1)
                        for d=1:5
                            con=zeros(1,size(SPM.xX.X,2));
                            con(1,T.digit==d & T.hand==h & T.stimType==s)=1;
                            
                            odd = con(1:length(T.run))' .* mod(T.run,2);
                            con(1:length(T.run)) = odd;
                            
                            con=con/sum(con);
                            conName= sprintf('ODD_T_hand_%i_stimType_%i_digit_%i', h,s,d);
                            SPM.xCon(idx)=spm_FcUtil('Set',conName, 'T', 'c',con',SPM.xX.xKXs);
                            idx = idx + 1;
                        end
                    end
                end
            end;
            
            % EVEN RUN: t-contrast for each finger / hand against rest
            for h=[0 1] % hand 0=> left 1=> right
                for s= [0 1] % stimType  0=> of 1=> passive mov
                    if ~(h==0 && s==1)
                        for d=1:5
                            con=zeros(1,size(SPM.xX.X,2));
                            con(1,T.digit==d & T.hand==h & T.stimType==s)=1;
                            
                            even    = con(1:length(T.run))' .* (1-mod(T.run,2));
                            con(1:length(T.run))    = even;
                            
                            con=con/sum(con);
                            conName= sprintf('EVEN_T_hand_%i_stimType_%i_digit_%i', h,s,d);
                            SPM.xCon(idx)=spm_FcUtil('Set',conName, 'T', 'c',con',SPM.xX.xKXs);
                            idx = idx + 1;
                        end
                    end
                end
            end;
            
            %____do the constrasts
            SPM=spm_contrasts(SPM,[1:length(SPM.xCon)]);
            
            % move the split-half contrasts to the relevant folder
            mkdir('contrast_splithalf');
            
            d1=dir('spm*');
            for i=1:length(d1)
                movefile(d1(i).name,fullfile(cd,'contrast_splithalf',d1(i).name));
            end;
            
            d2=dir('ess_*');
            for i=1:length(d2)
                movefile(d2(i).name,fullfile(cd,'contrast_splithalf',d2(i).name));
            end;
            
            d3=dir('con_*');
            for i=1:length(d3)
                movefile(d3(i).name,fullfile(cd,'contrast_splithalf',d3(i).name));
            end;
        end;
    case 'percent_signal'  %Calculate percent signal changes
        sn=varargin{1};
        glm=varargin{2};
        outname={'psc_L_motor.nii','psc_R_motor.nii','psc_R_sens.nii'};
        for s=sn
            cd(fullfile(baseDir,glmName{glm}, subj_name{s}));
            load SPM;
            T=load('SPM_info.mat');
            X=(SPM.xX.X(:,SPM.xX.iC));
            h=median(max(X)); % Height of response;
            P={};
            for p=SPM.xX.iB
                P{end+1}=sprintf('beta_%4.4d.img',p);       % get the intercepts and use them to calculate the basline (mean images)
            end;
            for con=1:3
                P{end+1}=sprintf('con_%4.4d.img',con+4);
                if (length(P)==9)
                    formula=sprintf('100.*%f.*i9./((i1+i2+i3+i4+i5+i6+i7+i8)/8)',h);
                else
                    formula=sprintf('100.*%f.*i8./((i1+i2+i3+i4+i5+i6+i7)/8)',h);
                end;
                spm_imcalc_ui(P,outname{con},formula,...
                    {0,[],spm_type(16),[]});        % Calculate percent signal change
            end;
            fprintf('%d: %3.3f\n',s,h);
        end;
    case 'EMG_make_emg'
        s = 1:length(emg_subj);
        for sn = s
            % (1) converting mvc from smr to mat files
            mvcFile = fullfile(emgDir,emg_subj{sn},[emg_subj{sn} '_' emg_subjID{sn} '_MVC.smr']);
            fid = fopen(mvcFile);
            SONImport(fid);
            fclose(fid);
            
            source = fullfile(emgDir,emg_subj{sn},[emg_subj{sn} '_' emg_subjID{sn} '_MVC.mat']);
            dest   = fullfile(emgDir,emg_subj{sn},[emg_subj{sn} '_mvc.mat']);
            movefile(source,dest);
            disp('MVC file converted');
            
            % (2) convert individual runs from smr to mat files
            for i = 1:emg_NumRuns(sn)
                dataFile = fullfile(emgDir,emg_subj{sn},[emg_subj{sn} '_' emg_subjID{sn} '_run' num2str(i) '.smr']);
                fid = fopen(dataFile);
                SONImport(fid);
                fclose(fid);
                
                source = fullfile(emgDir,emg_subj{sn},[emg_subj{sn} '_' emg_subjID{sn} '_run' num2str(i) '.mat']);
                dest   = fullfile(emgDir,emg_subj{sn},[emg_subj{sn} '_run' num2str(i) '.mat']);
                movefile(source,dest);
                disp(['Run ' num2str(i) ' converted']);
            end;
        end;
    case 'SUIT_isolate'             % Isolate the Cerebellum
        %df1_imana('SUIT_isolate',1)
        sn=varargin{1};
        %----make suit DIR and copy the anatomical over
        if (~exist(fullfile(anatomicalDir,subj_name{sn},'suit'),'dir'))
            mkdir(fullfile(anatomicalDir,subj_name{sn},'suit'));
            copyfile(fullfile(anatomicalDir,subj_name{sn},[subj_name{sn} '_anatomical.nii']), fullfile(anatomicalDir,subj_name{sn},'suit',[subj_name{sn} '_anatomical.nii']));
        end
        %----got to suit DIR
        cd(fullfile(anatomicalDir,subj_name{sn},'suit'))
        %----start SPM
        spm fmri
        %----do SUIT isolate
        suit_isolate([subj_name{sn} '_anatomical.nii'])
        %----Now adjust ...anatomical_pcereb_corr.nii by hand
    case 'SUIT_normalize'           % Normalise the the template
        %df1_imana('SUIT_normalize',1)
        sn=varargin{1};
        %----got to suit DIR
        cd(fullfile(anatomicalDir,subj_name{sn},'suit'))
        %----run suit normalize
        suit_normalize(['c_',subj_name{sn},'_anatomical.nii'], 'mask',['c_',subj_name{sn},'_anatomical_pcereb_corr.nii'])
    case 'SUIT_make_mask'           %Make the mask for the cerebellum
        sn=varargin{1};
        mask_all= fullfile(baseDir, glmName{1},subj_name{sn}, 'mask.img');
        mask_cerebllum= fullfile(anatomicalDir,subj_name{sn},'suit',['c_',subj_name{sn},'_anatomical_pcereb_corr.nii'])
        mask_suit= fullfile(baseDir, glmName{1},subj_name{sn}, 'mask_suit.nii');
        %----calc suit mask
        spm_imcalc_ui({mask_all,mask_cerebllum}, mask_suit, 'i1>0 & i2>0');
    case 'SUIT_reslice'
        files={'dist_L_motor.nii','dist_R_motor.nii','dist_R_sens.nii'};
        sn=varargin{1};
        mask   = fullfile(anatomicalDir,subj_name{sn},'suit',['c_',subj_name{sn},'_anatomical_pcereb_corr.nii'])
        params = fullfile(anatomicalDir,subj_name{sn},'suit',['mc_',subj_name{sn},'_anatomical_snc.mat'])
        for i=1:length(files)
            P=fullfile(glmDir,subj_name{sn},files{i});
            O=fullfile(baseDir,'SUIT_group_data',[subj_name{sn} '_' files{i}]);
            suit_reslice(P,params,'outfilename',O,'mask',mask);
        end;
    case 'MVA_search_surf'             % Volume-based search light for surface
        sn=varargin{1};
        atlas=2;
        refDir= fullfile(baseDir,glmName{1});
        
        for s=sn
            for h=1:2
                caret_subjDIR = fullfile(caretDir,[atlasA{atlas},subj_name{s}],hemName{h});
                coord_pial= caret_load(fullfile(caret_subjDIR, [hem{h} '.PIAL.coord']));
                coord_white= caret_load(fullfile(caret_subjDIR, [hem{h} '.WHITE.coord']));
                topo= caret_load(fullfile(caret_subjDIR, [hem{h} '.CLOSED.topo']));
                surf(h).c1=coord_white.data';
                surf(h).c2=coord_pial.data';
                surf(h).f=topo.data';
            end;
            clear coord_pial coord_white coord_caret;
            volDef= spm_vol(fullfile(refDir, subj_name{s},'mask.img'));
            volDef.mask=spm_read_vols(volDef);
            [LI,voxmin,voxmax,voxel,node,surfindx,depth,radvox]= lmva_voxelselection_surf(surf, [6 80],volDef);
            LI=LI';
            save(fullfile(refDir,subj_name{s}, 'vol_roi_80vox.mat'), 'LI','voxmin','voxmax','voxel','radvox');
            save(fullfile(refDir,subj_name{s}, 'vol_surf.mat'), 'voxel','node','surfindx','depth');
            
            %             V=volDef;
            %             X=zeros(V.dim);
            %             depth=depth+1;
            %             depth(isnan(depth))=1;
            %             X(voxel)=depth;
            %             V=volDef;
            %             V.dt=[4 0];
            %             V.pinfo=[3/100 0 0]';
            %             V.fname=fullfile(refDir,subj_name{s},'vol_roi_depth.nii');
            %             spm_write_vol(V,X);
        end;
    case 'MVA_search_cereb'            % Add cerebellar searchlight
        sn=varargin{1};
        atlas=2;
        refDir= fullfile(baseDir,glmName{1});
        
        for s=sn
            % Get cerebellar mask
            cereb_mask=fullfile(anatomicalDir,subj_name{s},'suit',['c_' subj_name{s} '_anatomical_pcereb_corr.nii']);
            Vcereb_mask=spm_vol(cereb_mask);
            
            % Get gray matter segmentation
            gray_matter=fullfile(anatomicalDir,subj_name{s},'suit',[subj_name{s} '_anatomical_seg1.nii']);
            Vgray_matter=spm_vol(gray_matter);
            
            % Function mask image
            volDef= spm_vol(fullfile(refDir, subj_name{s},'mask.img'));
            volDef.mask=spm_read_vols(volDef);
            indx=find(volDef.mask>0);
            voxels=surfing_inds2subs(volDef.dim,indx);
            voxels=[voxels ones(size(voxels,1),1)];
            
            voxelsAna=(inv(Vcereb_mask.mat)*volDef.mat*voxels')';
            A=spm_sample_vol(Vcereb_mask,voxelsAna(:,1),voxelsAna(:,2),voxelsAna(:,3),1);
            
            voxelsSeg=(inv(Vgray_matter.mat)*volDef.mat*voxels')';
            B=spm_sample_vol(Vgray_matter,voxelsSeg(:,1),voxelsSeg(:,2),voxelsSeg(:,3),1);
            
            volDef.mask(indx(A<0.5 | B<0.2))=0;
            
            % Load the cortical searchlights
            searchlightname=fullfile(refDir,subj_name{s},'vol_roi_80vox.mat');
            S=load(searchlightname);
            
            % Kill the voxels that belong to cortical searchlights
            volDef.mask(S.voxel)=0;
            
            % Compute the search lights for the cerebellum and add them to
            % the file
            [LI,voxmin,voxmax,voxel,radvox]= lmva_voxelselection_volume([20 80],volDef);
            S.LI=[S.LI;LI];
            S.voxel=[S.voxel;voxel];
            S.voxmin=[S.voxmin;voxmin];
            S.voxmax=[S.voxmax;voxmax];
            S.radvox=[S.radvox;radvox];
            save(fullfile(refDir,subj_name{s}, 'vol_roi_80vox.mat'), '-struct','S');
        end;
    case 'MVA_do' %_________________% Do the multivariate classification as searchlight
        % df1_imana('MVA_do',1:4, 60,2)
        sn=varargin{1};
        
        workDir= fullfile(baseDir, glmName{1});
        ldafunction= @df1_calcAcc;
        
        for s=sn
            cd(fullfile(workDir, subj_name{s}));
            T= load(fullfile(baseDir,glmName{1}, subj_name{s},'SPM_info.mat'));
            load('SPM');
            %define beta images
            beta_images=  {SPM.Vbeta(SPM.xX.iC).fname}';
            
            lda_names= {fullfile(workDir, subj_name{s}, 'acc_L_motor.nii'),...
                fullfile(workDir, subj_name{s}, 'acc_R_motor.nii'),...
                fullfile(workDir, subj_name{s}, 'acc_R_sens.nii')};
            %----do LDA
            lmva_spm('vol_roi_80vox.mat',beta_images,lda_names,ldafunction,'params',{T.run',T.hand',T.stimType',T.digit'})%
        end
    case 'MVA_do_dist'             % Does the distance analysis as the search light
        sn=varargin{1};
        glm=1;
        ldafunction= @df1_calcDist;
        
        for s=sn
            glmDirSubj=fullfile(baseDir,glmName{glm},subj_name{s});
            load(fullfile(glmDirSubj,'SPM.mat'));
            D= load(fullfile(glmDirSubj,'SPM_info.mat'));
            SPM=spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s}));
            
            outfiles = {fullfile(glmDirSubj, 'dist_L_motor.nii'),...
                fullfile(glmDirSubj, 'dist_R_motor.nii'),...
                fullfile(glmDirSubj, 'dist_R_sens.nii')};
            
            surfaces=load(fullfile(glmDirSubj,['vol_roi_80vox','.mat']));
            lmva_spm(surfaces,SPM.xY.P,outfiles,ldafunction,'params',{SPM,D});
        end;
    case 'MVA_zValue' %______________DF____________________________________
        %df1_imana('MVA_zValue',3:5,1)
        sn=varargin{1};
        glmType=varargin{2}
        if glmType==1
            workDir= fullfile(baseDir, 'glm_firstlevel_1');
            numRegressors= 32;
        elseif glmType==4
            workDir= fullfile(baseDir, 'glm_firstlevel_3');
            numRegressors= 32*3;
        end
        mu= 1/4;
        
        sigma= sqrt(1/4 * 3/4 * 1/numRegressors);
        images= {['lda_L_motor',num2str(numVox),'.nii'], ['lda_R_motor',num2str(numVox),'.nii'], ['lda_R_sens',num2str(numVox),'.nii']}
        outimages= {['zlda_L_motor',num2str(numVox),'.nii'], ['zlda_R_motor',num2str(numVox),'.nii'], ['zlda_R_sens',num2str(numVox),'.nii']}
        
        for s=sn;
            for j=1:numel(images)
                cd(fullfile(workDir, subj_name{s}));
                spmj_imcalc_mtx(images{j}, outimages{j},...
                    sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask!
            end
        end
        %------------------------------------------
        %----group statistic-----------------------
        %------------------------------------------
    case 'percent_signal_change'
        % df1_imana('percent_signal_change', 1:11, 'untrained')
        % do it for the glm_firstlevel only
        sn=varargin{1}; seqCondition=varargin{2};
        scale=1.2600; % Height of the HRF in the design Matrix plit regressors and look plot(SPM.xX.X(1:100,4))
        for s=sn
            cd(fullfile(baseDir,'glm_firstlevel_1', subj_name{s}, seqCondition));
            load SPM;
            P={};
            for p=SPM.xX.iB
                P{end+1}=sprintf('beta_%4.4d.img',p);       % get the intercepts and use them to calculate the basline (mean images)
            end;
            %P_intercept=sprintf('perc_00_%s.img',subj_name{s}); % Percent 0 is baseline
            P_intercept=sprintf('mean_intercept.nii'); % Percent 0 is baseline
            spm_imcalc_ui(char(P),P_intercept,'mean(X)',{1,[],spm_type(16),[]});
            for con=6
                conname=sprintf('con_%4.4d.img',con);
                outname=sprintf('psc_con_%4.4d.nii',con); % ,subj_name{s}
                formula=sprintf('100.*%f.*i2./i1',scale);
                spm_imcalc_ui(char({P_intercept,conname}),outname,formula,...
                    {0,[],spm_type(16),[]});        % Calculate percent signal change
            end;
        end
    case 'MNI_normalization_write' %___________________________________
        %df1_imana('MNI_normalization_write',3:5, 1)
        sn=varargin{1}; glmType=varargin{2}
        if glmType==1
            workDir= fullfile(baseDir, 'glm_firstlevel_1');
            outPrefix=[];
            images={'con_0001.img','con_0005.img','con_0006.img','con_0007.img','mask.img', 'lda_L_motor60.nii', 'lda_R_motor60.nii', 'lda_R_sens60.nii', ...
                'zlda_L_motor60.nii', 'zlda_R_motor60.nii', 'zlda_R_sens60.nii'} %
            %         elseif glmType==2
            %             workDir= fullfile(baseDir, 'glm_firstlevel_3');
            %             outPrefix='_3reg';
            %             images={''} %
        end
        for s=sn
            defor= fullfile(anatomicalDir, subj_name{s}, [subj_name{s}, '_anatomical_seg_sn.mat']);
            for j=1:numel(images)
                [d,name,ext]=spm_fileparts(images{j});
                sn_images{j}= fullfile(workDir,subj_name{s},images{j});
                
                out_images{j}= fullfile(groupData,[name, '_' subj_name{s}, outPrefix,  '.nii']);
            end
            spmj_normalization_write(defor, sn_images,'outimages',out_images);
            
            % Do the anatomical
            sn_images={}; out_images={};
            sn_images{1}=fullfile(baseDir, 'anatomicals', subj_name{s},[subj_name{s},'_anatomical.nii']);
            out_images{1}=fullfile(groupData,[subj_name{s}, '_anatomical.nii']);
            spmj_normalization_write(defor, sn_images,'outimages',out_images);
        end;
    case 'MNI_mask_anatomical' %__________________________________________________
        %df1_imana('MNI_mask_anatomical', 3:5)
        sn=varargin{1};
        cd(groupData);
        i=1;
        for s=sn
            images{i}= fullfile(groupData,['mask_',subj_name{s},  '.nii']);
            ana_images{i}= fullfile(groupData,[subj_name{s}, '_anatomical.nii']);
            i=i+1;
        end
        spmj_imcalc_mtx(images,'mask_avrg.nii','nanmean(X)');
        spmj_imcalc_mtx(ana_images,'avrg_anatomical.nii','nanmean(X)');
        spmj_imcalc_mtx('mask_avrg.nii','mask_thres.nii','X>0.4');
    case 'MNI_smooth' %________________________________________________
        %df1_imana('MNI_smooth',3:5)
        sn=varargin{1}
        images= {'zlda_L_motor60', 'zlda_R_motor60', 'zlda_R_sens60','con_0001','con_0005','con_0006','con_0007'}
        nam=[];
        for s=sn
            for j=1:numel(images)
                nam{end+1} = fullfile(groupData, [images{j},'_',subj_name{s},'.nii,1']);
            end
        end
        matlabbatch{1}.spm.spatial.smooth.data =  cellstr(char(nam));
        matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6 ];
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run',matlabbatch);
    case 'MNI_make_rfx'
        %df1_imana('MNI_make_rfx',[1:3])
        sn=varargin{1};
        images= {'szlda_L_160', 'szlda_R_160','scon_0004','scon_0005','scon_0006'}
        nam=[];
        for j=1:numel(images)
            mkdir(fullfile(groupDir,[images{j}])); %,'_4reg'
            matlabbatch{1}.spm.stats.factorial_design.dir = {fullfile(groupDir,[images{j}])} %,'_4reg'
            scans=[];
            for s=sn
                scans{end+1} = fullfile(groupData, [images{j},'_',subj_name{s}, '.nii,1']); %_4reg
            end
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans= scans;
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {fullfile(groupData,'mask_thres.nii,1')};
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            spm_jobman('run',matlabbatch);
        end
    case 'MNI_rfx_spm'
        %df1_imana('MNI_rfx_spm', 1:5)
        con= {'szlda_L_160', 'szlda_R_160','scon_0004','scon_0005','scon_0006'}
        contr=varargin{1};
        for k=contr
            dirname=fullfile(groupDir,[con{k}]);
            cd(dirname);
            load SPM;
            SPM=spm_spm(SPM);
            spmj_rfx_contrast(SPM);
            cd ..
        end;
        %------------------------------------------
        %----surface analysis----------------------
        %------------------------------------------
    case 'surf_makeall'
        sn=varargin{1};
        df1_imana('surf_freesurfer', sn);
        df1_imana('surf_xhemireg',sn);
        df1_imana('surf_map_ico',sn,2);
        df1_imana('surf_make_caret',sn,2);
    case 'surf_freesurfer' %_________________1_________________________
        %df1_imana('surf_freesurfer', 5)
        sn=varargin{1};
        freesurfer_reconall(freesurferDir,subj_name{sn},fullfile(anatomicalDir,subj_name{sn},[subj_name{sn} '_anatomical.nii']));
    case 'surf_xhemireg'
        % df1_imana('surf_xhemireg',1)
        sn=varargin{1};
        for i=sn
            freesurfer_registerXhem({subj_name{i}},freesurferDir,'hemisphere',[1:2]);
        end;
    case 'surf_map_ico' %____________DF________________________________________
        % df1_imana('surf_map_ico',3:5,2)
        sn=varargin{1};
        atlas=varargin{2};
        if (atlas==2)
            for i=sn
                freesurfer_mapicosahedron_xhem(subj_name{i},freesurferDir,'smoothing',1,'hemisphere',[1:2]);
            end;
        else
            for i=sn
                freesurfer_mapicosahedron(subj_name{i},freesurferDir,'smoothing',1);
            end;
        end;
        %=========================================================================
        % Transfer data into Caret-format
        %=========================================================================
    case 'surf_make_caret' %_________DF________________________________________
        % df1_imana('surf_make_caret',3:5,2)
        sn=varargin{1};
        atlas=varargin{2};
        for i=sn
            caret_importfreesurfer([atlasA{atlas} subj_name{i}],freesurferDir,caretDir);
        end;
        %=========================================================================
        % Map functional volumes on surface over caret
        %=========================================================================
        
    case 'surf_map_con' %_________DF______________New mappint algorithm_____
        %df1_imana('surf_map_con', 3:5)
        % map volume images to metric file and save them in individual surface folder
        sn=varargin{1};
        hemisphere=[1:2];
        atlas=2;
        
        vararginoptions({varargin{2:end}},{'atlas','hemisphere'});
        glm=1;
        fileList={'con_0005.img','con_0006.img', 'con_0007.img', 'spmT_0005.img', 'spmT_0006.img', 'spmT_0007.img'};    %right hand sensory 20:24
        for s=sn
            for h=hemisphere
                caretSDir = fullfile(caretDir,[atlasA{atlas},subj_name{s}],hemName{h});
                specname=fullfile(caretSDir,[atlasA{atlas},subj_name{s} '.' hem{h}   '.spec']);
                white=fullfile(caretSDir,[hem{h} '.WHITE.coord']);
                pial=fullfile(caretSDir,[hem{h} '.PIAL.coord']);
                
                C1=caret_load(white);
                C2=caret_load(pial);
                
                for f=1:length(fileList)
                    images{f}=fullfile(baseDir, glmName{glm},subj_name{s},fileList{f});
                end;
                M=caret_vol2surf_own(C1.data,C2.data,images,'ignore_zeros',1);
                caret_save(fullfile(caretSDir,[subj_name{s} '_func.metric']),M);
            end;
        end;
    case 'surf_map_dist'   % Map the distances between fingers across the cortex
        %df1_imana('surf_map_con', 3:5)
        % map volume images to metric file and save them in individual surface folder
        sn=varargin{1};
        hemisphere=[1:2];
        atlas=2;
        
        vararginoptions({varargin{2:end}},{'atlas','hemisphere'});
        glm=1;
        fileList={'dist_L_motor.nii','dist_R_motor.nii', 'dist_R_sens.nii'};    %right hand sensory 20:24
        for s=sn
            for h=hemisphere
                fprintf('%d %d\n',s,h);
                caretSDir = fullfile(caretDir,[atlasA{atlas},subj_name{s}],hemName{h});
                specname=fullfile(caretSDir,[atlasA{atlas},subj_name{s} '.' hem{h}   '.spec']);
                white=fullfile(caretSDir,[hem{h} '.WHITE.coord']);
                pial=fullfile(caretSDir,[hem{h} '.PIAL.coord']);
                
                C1=caret_load(white);
                C2=caret_load(pial);
                
                for f=1:length(fileList)
                    images{f}=fullfile(baseDir, glmName{glm},subj_name{s},fileList{f});
                end;
                M=caret_vol2surf_own(C1.data,C2.data,images,'ignore_zeros',1);
                caret_save(fullfile(caretSDir,[subj_name{s} '_dist.metric']),M);
            end;
        end;
        
    case 'surf_map_finger' %______DF______________New mappint algorithm_____
        %df1_imana('surf_map_finger', 3:5)
        % map volume images to metric file and save them in individual surface folder
        sn=varargin{1};
        hemisphere=[1:2];
        atlas=2;
        
        vararginoptions({varargin{2:end}},{'atlas','hemisphere'});
        glm=1;
        fileList={...
            'spmT_0008.img','spmT_0009.img','spmT_0010.img','spmT_0011.img','spmT_0012.img', ... %left hand motor    1:5
            'spmT_0013.img','spmT_0014.img','spmT_0015.img','spmT_0016.img','spmT_0017.img', ... %right hand motor   6:10
            'spmT_0018.img','spmT_0019.img','spmT_0020.img','spmT_0021.img','spmT_0022.img'};    %right hand sensory 11:15
        for s=sn
            for h=hemisphere
                caretSDir = fullfile(caretDir,[atlasA{atlas},subj_name{s}],hemName{h});
                specname=fullfile(caretSDir,[atlasA{atlas},subj_name{s} '.' hem{h}   '.spec']);
                white=fullfile(caretSDir,[hem{h} '.WHITE.coord']);
                pial=fullfile(caretSDir,[hem{h} '.PIAL.coord']);
                topo=fullfile(caretSDir,[hem{h} '.CLOSED.topo']);
                
                C1=caret_load(white);
                C2=caret_load(pial);
                
                for f=1:length(fileList)
                    images{f}=fullfile(baseDir, glmName{glm},subj_name{s},fileList{f});
                end;
                % M=caret_vol2surf_own(C1.data,C2.data,images,'ignore_zeros',1);
                M=caret_vol2surf_own(C1.data,C2.data,images,'topo',topo,'exclude_thres',0.75,'ignore_zeros',1);
                caret_save(fullfile(caretSDir,[subj_name{s} '_finger.metric']),M);
            end;
        end;
    case 'surf_map_finger_splithalf'
        % df1_imana('surf_map_finger', 3:5)
        % map volume images to metric file and save them in individual surface folder
        sn=varargin{1};
        hemisphere=[1:2];
        atlas=2;
        
        vararginoptions({varargin{2:end}},{'atlas','hemisphere'});
        glm=1;
        
        % map contrasts from odd experimental runs
        fprintf('ODD\n---\n');
        fileList={...
            'spmT_0001.nii','spmT_0002.nii','spmT_0003.nii','spmT_0004.nii','spmT_0005.nii', ... %left hand motor    1:5
            'spmT_0006.nii','spmT_0007.nii','spmT_0008.nii','spmT_0009.nii','spmT_0010.nii', ... %right hand motor   6:10
            'spmT_0011.nii','spmT_0012.nii','spmT_0013.nii','spmT_0014.nii','spmT_0015.nii'};    %right hand sensory 11:15
        for s=sn
            disp(s);
            for h=hemisphere
                caretSDir = fullfile(caretDir,[atlasA{atlas},subj_name{s}],hemName{h});
                specname=fullfile(caretSDir,[atlasA{atlas},subj_name{s} '.' hem{h}   '.spec']);
                white=fullfile(caretSDir,[hem{h} '.WHITE.coord']);
                pial=fullfile(caretSDir,[hem{h} '.PIAL.coord']);
                topo=fullfile(caretSDir,[hem{h} '.CLOSED.topo']);
                
                C1=caret_load(white);
                C2=caret_load(pial);
                
                for f=1:length(fileList)
                    images{f}=fullfile(baseDir, glmName{glm},subj_name{s},'contrast_splithalf',fileList{f});
                end;
                % M=caret_vol2surf_own(C1.data,C2.data,images,'ignore_zeros',1);
                M=caret_vol2surf_own(C1.data,C2.data,images,'topo',topo,'exclude_thres',0.75,'ignore_zeros',1);
                caret_save(fullfile(caretSDir,[subj_name{s} '_finger_odd.metric']),M);
            end;
        end;
        
        % map contrasts from even experimental runs
        fprintf('EVEN\n----\n');
        fileList={...
            'spmT_0016.nii','spmT_0017.nii','spmT_0018.nii','spmT_0019.nii','spmT_0020.nii', ... %left hand motor    16:20
            'spmT_0021.nii','spmT_0022.nii','spmT_0023.nii','spmT_0024.nii','spmT_0025.nii', ... %right hand motor   21:25
            'spmT_0026.nii','spmT_0027.nii','spmT_0028.nii','spmT_0029.nii','spmT_0030.nii'};    %right hand sensory 26:30
        for s=sn
            disp(s);
            for h=hemisphere
                caretSDir = fullfile(caretDir,[atlasA{atlas},subj_name{s}],hemName{h});
                specname=fullfile(caretSDir,[atlasA{atlas},subj_name{s} '.' hem{h}   '.spec']);
                white=fullfile(caretSDir,[hem{h} '.WHITE.coord']);
                pial=fullfile(caretSDir,[hem{h} '.PIAL.coord']);
                topo=fullfile(caretSDir,[hem{h} '.CLOSED.topo']);
                
                C1=caret_load(white);
                C2=caret_load(pial);
                
                for f=1:length(fileList)
                    images{f}=fullfile(baseDir, glmName{glm},subj_name{s},'contrast_splithalf',fileList{f});
                end;
                % M=caret_vol2surf_own(C1.data,C2.data,images,'ignore_zeros',1);
                M=caret_vol2surf_own(C1.data,C2.data,images,'topo',topo,'exclude_thres',0.75,'ignore_zeros',1);
                caret_save(fullfile(caretSDir,[subj_name{s} '_finger_even.metric']),M);
            end;
        end;
    case 'surf_avrgsurfshape' %___DF______________Generates group metric____
        % Make group metric for surface shape into symmetric template (left template)
        %df1_imana('surf_avrgsurfshape')
        INname={'surface_shape','surface_shape','surface_shape'};
        OUTname={'curv.surface_shape','sulc.surface_shape','area.surface_shape'};
        inputcol=[1 2 3];
        atlas=2;
        sn=[3:length(subj_name)];
        
        vararginoptions(varargin,{'atlas','sn'});
        for h=1:2
            surfaceGroupDir=[caretDir filesep atlasname{atlas}  filesep hemName{h} ];
            cd(surfaceGroupDir);
            for j=1:length(INname); %----loop over each input metric file and make a group metric file
                k=1;
                for i=sn; %----define names of subj metric files
                    infilenames{j}{k}=[caretDir filesep atlasA{atlas} subj_name{i} filesep hemName{h} filesep hem{h} '.' INname{j}];
                    k=k+1;
                end;
                %----name for group metric file in average surface folder
                outfilenames{j}=[surfaceGroupDir filesep hem{h} '.' OUTname{j}];
                %----make the group metric
                fprintf('hem: %i  image: %i \n', h,j);
                caret_metricpermute(infilenames{j},'outfilenames',outfilenames(j),'inputcol',inputcol(j));
            end;
        end;
    case 'surf_fingerpatterns' %__DF______________Plot finger patterns______
        %df1_imana('surf_fingerpatterns' , 6, 1)
        %close all;
        sn=varargin{1};
        h=varargin{2};
        groupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
        cd(groupDir);
        border=fullfile(caretDir,'fsaverage_sym',hemName{h},['CS.border']);
        switch(h)
            case 1
                coord='lh.FLAT.coord';
                topo='lh.CUT.topo';
                data='lh.surface_shape';
                xlims=[-20 10];
                ylims=[-15 30];
            case 2
                coord='rh.FLAT.coord';
                topo='rh.CUT.topo';
                data='rh.surface_shape';
                xlims=[-10 20];
                ylims=[-15 30];
                
        end;
        
        B=caret_load(border);
        
        data=fullfile(caretDir,['x' subj_name{sn}],hemName{h},[subj_name{sn} '_finger.metric']);
        sshape=fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.surface_shape']);
        %         subplot(2,3,1);
        M=caret_plotflatmap('col',2,'data',sshape,'border',B.Border,'topo',topo,'coord',coord,'xlims',xlims,'ylims',ylims);
        %         colormap('gray')
        %         figure
        %         colormap('default')
        for i=1:15 %----left motor----right motor----right sensory
            subplot(3,5,i);
            [M,d]=caret_plotflatmap('M',M,'col',i,'data',data,'cscale',[-6 12],...
                'border',B.Border,'topo',topo,'coord',coord);
            maxT(i)=max(d(:));
        end;
        mm=max(maxT);
        for i=1:15
            subplot(3,5,i);
            caxis([-mm/2 mm]);
        end;
        set(gcf,'PaperPosition',[1 1 12 7]);
        wysiwyg;
        saveas(gcf, [subj_name{sn},'_',hemName{h}], 'jpg')
    case 'surf_spatial_cog'   %_____ToDo______________________________________
        %df1_imana('surf_spatial_cog', 3, 1)
        sn=varargin{1};
        reg=varargin{2};%----which region are you looking at
        
        switch reg
            case 1 %M1
                %lims{1}=[-70 -30 30 70 ];
                lims{1}=[-1 0 0  1 ];
                %lims{1}=[-100 0 0 100 ];
                lims{2}=[-1 0 0  -1 ];
                %lims{2}=[0 40 45 85];
                %lims{2}=[-100 0 0 100 ];
                paintROI=1; %----Which ROI should be used
                flat{1}= '.FLAT.coord';
                flat{2}= '.CUT.topo';
                metric_file{1}= '.func.metric';
                %metric_file{1}= '.zacc_trained_stats.metric';
                %metric_file{2}= '.zacc_untrained_stats.metric';
                group_th=2.5; %---- threshold for the plotted group average activation
                subj_th=0.8 %0.4; % ---- additional threshold for the individual ROI region in which the COG is calculated
        end
        AllT=[];
        drawpoints=1;
        drawborder=0;
        
        set(gcf,'PaperPosition',[2 2 8 4]);
        wysiwyg;
        for h=1:2
            T=[];
            groupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} filesep];
            cd(groupDir);
            %----load the surface information
            coord=([groupDir hem{h} flat{1}]);
            topo=([groupDir hem{h} flat{2}]);
            FLAT=caret_load(coord);
            %----load the shape file for visualisation
            S=caret_load([groupDir hem{h} '.surface_shape']);
            %----go to subj folder
            cd([caretDir filesep 'x' subj_name{sn} filesep hemName{h}])
            %----load the functional data
            C=caret_load([subj_name{sn},'_func.metric']);
            %----get data
            DATA= C.data(:,10);
            %---set small values to zeros => they will not be presented in the RGB map
            %DATA(DATA<group_th)=0;
            %----
            DSCALE=[0.97 5.7;0 0;0.97 1.5];%DSCALE=[0.97 5.7;0 0;0.97 5.7];DSCALE=[0.97 5.7;0 0;0.97 1.5];DSCALE=[2.45 5.7;0 0;0.97 1.5]; % 50% and 40%
            
            %             %----load the ROI file in which SMA and dPM is defined
            %             P=caret_load([hem{h} '.ROI_cog.paint']);
            %----in case you like to draw boarders as well
            drawborder=0;
            if (drawborder)
                B=caret_load([hem{h} '.borderdisp_tobi2.border']);%'.dispborder.border']);
            else
                B.Border=[];
            end;
            
            subplot(1,2,h);
            %----plot the backround image with the functional data in the backround
            caret_plotflatmap_rgb('coord',coord,'topo',topo,'underlay',S.data(:,2),...
                'data',DATA,...
                'dscale',DSCALE,...
                'uscale',[-3 1],...
                'xlims',lims{h}(1:2),'ylims',lims{h}(3:4),'alpha',0.8,...
                'border',B.Border);
            %----Determine individual COGs
            for s=subj
                D.SN=s;
                D.hem=h;
                %----IMPORTANT!!!!
                %----index to ROI  and
                %----if subj_th is not NaN also threshold on both tr and untr data to select the nodes
                %----that are going into the COG calculation
                % version
                %
                CT.data((CT.data(:,s)< subj_th), s)=NaN;
                CU.data((CU.data(:,s)< subj_th), s)=NaN;
                indx=(P.data==paintROI & ~isnan(CT.data(:,s)) & ~isnan(CU.data(:,s)) );
                %----alterantive
                % CT.data((CT.data(:,s)< subj_th & (DATA(:,1)~=0| DATA(:,3)~=0)), s)=NaN;
                % CU.data((CU.data(:,s)< subj_th & (DATA(:,1)~=0| DATA(:,3)~=0)), s)=NaN;
                % NEW do it on the group mask
                % indx=(P.data==paintROI & ~isnan(CT.data(:,s))  & ~isnan(CU.data(:,s)) & (DATA(:,1)~=0 | DATA(:,3)~=0));
                %sum(indx)
                %----calculate the COGs in the subset of nodes!
                D.COG_U=sl1_COG(CU.data(indx,s),FLAT.data(indx,:));
                D.COG_T=sl1_COG(CT.data(indx,s),FLAT.data(indx,:));
                T=addstruct(T,D);
            end;
            %----add the cog's to the image
            if (drawpoints)
                hold on;
                c1=plot(T.COG_U(:,1),T.COG_U(:,2),'ko');
                c2=plot(T.COG_T(:,1),T.COG_T(:,2),'ko'); %+*.xsd^v><ph
                set(c1,'MarkerFaceColor',[0.3 0.3 1]);
                set(c2,'MarkerFaceColor',[1 0.3 0.3]);
                hold off;
            end;
            %----add average COG
            hold on;
            c1=plot(mean(T.COG_U(:,1)),mean(T.COG_U(:,2)),'ks');
            c2=plot(mean(T.COG_T(:,1)),mean(T.COG_T(:,2)),'ks'); %+*.xsd^v><ph
            set(c1,'MarkerFaceColor',[0.3 0.3 1]);
            set(c2,'MarkerFaceColor',[1 0.3 0.3]);
            hold off;
            
            axis equal;
            fprintf('Anterior\n');
            X=[T.COG_U-T.COG_T];
            T2Hot1(X(:,1:2),0.1)
            
            AllT=addstruct(AllT,T);
            [T.COG_U T.COG_T T.COG_U(:,1)-T.COG_T(:,1)]
        end;
        
        varargout={AllT};
    case 'surf_make_ROIpaint'   %___DF______________Generates surface ROI paint file
        %df1_imana('surf_make_ROIpaint')
        for h=1:2 % loop over hemispheres
            groupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
            cd(groupDir);
            M=caret_load([hem{h} '.propatlas.metric']);
            P=caret_load(['lateral.paint']);  % Column 2 is indicator for lateral surface
            S1=sum(M.data(:,[1 2 3 4]),2);
            M1=sum(M.data(:,[7 8]),2);
            PM=M.data(:,9);
            [Prop,ROI]=max([S1 M1 PM],[],2);
            ROI(Prop<0.3)=0;
            ROI(ROI==1 & P.data(:,2)==0)=0;
            ROI(ROI==2 & P.data(:,2)==0)=0;
            ROI(ROI==3 & P.data(:,2)==0)=4;
            ROI(ROI==3 & P.data(:,2)==0)=0;
            ROI(P.data(:,1)==1)=5;
            
            % define Hand area - based on cutoff on flat area and propatlas
            C=caret_load([hem{h} '.FLAT.coord']);
            ROI(:,2)=double(M1>0.2);
            ROI(C.data(:,2)<-15,2)=0;
            ROI(C.data(:,2)>25,2)=0;
            % Save Paint file
            Paint=caret_struct('paint','data',ROI,'paintnames',{'S1','M1','PM','SMA','SPL'},'column_name',{'ROI','handROI'});
            caret_save(['ROIs_BA.paint'],Paint);
        end
    case 'surf_makeGroup'   %_______
        atlas=2;
        
        INname={'func','func','func',...
            'dist','dist','dist'};
        OUTname={'func_L_motor','func_R_motor',...
            'func_R_sens',...
            'dist_L_motor','dist_R_motor','dist_R_sens'};
        
        inputcol= [1 2 3 1 2 3];
        replaceNaN=[1 1 1 1 1 1];
        vararginoptions(varargin,{'atlas'});
        
        for h=1:2
            surfaceGroupDir=[caretDir filesep atlasname{atlas} filesep hemName{h} ];
            cd(surfaceGroupDir);
            for j=1:length(INname); %----loop over each input metric file and make a group metric file
                for i=1:length(subj_name); %----define names of subj metric files
                    infilenames{j}{i}=[caretDir filesep atlasA{atlas} subj_name{i} filesep hemName{h} filesep subj_name{i} '_' INname{j} '.metric'];
                end;
                %----name for group metric file in average surface folder
                outfilenames{j}=[surfaceGroupDir filesep hem{h} '.' OUTname{j} '.metric'];
                %----make the group metric
                fprintf('hem: %i  image: %i \n', h,j);
                caret_metricpermute(infilenames{j},'outfilenames',outfilenames(j),'inputcol',inputcol(j),'replaceNaNs',replaceNaN(j));
            end;
        end;
    case 'surf_zacc'  % Make the zacc files (on group level) and remove the NaN
        %df1_imana('surf_zacc')
        atlas =2;
        INname={'acc_left','acc_right'};
        OUTname={'zacc_left','zacc_right'};
        vararginoptions(varargin,{'atlas'});
        
        for h=1:2
            surfaceGroupDir=[caretDir filesep atlasname{atlas} filesep hemName{h} ];
            cd(surfaceGroupDir);
            for i=1:length(INname)
                iname=[hem{h} '.' INname{i} '.metric'];
                oname=[hem{h} '.' OUTname{i} '.metric'];
                fprintf('hem: %i  image: %i \n', h,i);
                caret_metriccalc({iname},oname,'(i1-0.25)/0.0765466','replaceNaNs',1);
            end;
        end;
    case 'surf_groupSmooth'           % Smooth the group maps
        %df1_imana('surf_groupSmooth', 'Smooth_iterations', 5)
        Smooth_iterations=[20 20 20 10 10 10];
        SPMname={'func_L_motor','func_R_motor',...
            'func_R_sens',...
            'dist_L_motor','dist_R_motor','dist_R_sens'};
        for h=1:2
            surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h}];
            cd(surfaceGroupDir)
            %----define name of coord and topology
            coordfile=[surfaceGroupDir filesep hem{h} '.FIDUCIAL.coord'];
            topofile=[surfaceGroupDir filesep hem{h} '.CLOSED.topo'];
            
            %----get the full directory name of the metric files and the smoothed metric files that we create below
            for i=1:length(SPMname);
                filenames{i}=[surfaceGroupDir filesep hem{h} '.' SPMname{i} '.metric']; % unsmoothed
                sfilenames=caret_smooth(filenames{i},'coord',coordfile,'topo',topofile,'iterations', Smooth_iterations(i));
            end;
            %----smooth the metric files and save them with the prefix 's'
        end;
    case 'surf_groupCon'              % Different group-GLM models  ....
        %sl1_pattern_behavior('surf_mapCorr',1:3)
        SPMname={'func_L_motor','func_R_motor',...
            'func_R_sens',...
            'dist_L_motor','dist_R_motor','dist_R_sens'};
        vararginoptions(varargin,{'SPMname'});
        contrasts = [1 0;0 1;-1 1];
        contrastName = {'Control','Dystonic','Dys>Con'};
        numContrasts = 3;
        numSPM=length(SPMname);
        numGroups = 2;
        for h=1:2
            surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
            cd(surfaceGroupDir);
            data=[];column_name={};
            outDir=[surfaceGroupDir filesep ];
            X=[];
            for i=1:numGroups
                X=[X double((subj_group==i))];
            end;
            for s=1:length(SPMname)
                %----get the full directory name of the metric files and the smoothed metric files that we create below
                filename=['s' hem{h} '.' SPMname{s} '.metric']; % smoothed
                outname=SPMname{s};
                cSPM=caret_getcSPM('regression',X,'data',filename,[1:length(subj_name)],'no_intercept','contrasts',contrasts);
                save([outDir, 'cSPM.' outname, '.mat'],'cSPM');
                caret_savecSPM([outDir, 'stats.',outname,'.metric'],cSPM);
                for g=1:numContrasts
                    data(:,(s-1)*numContrasts+g)=cSPM.con(g).con; % mean of group
                    data(:,(s-1)*numContrasts+g+numSPM*numContrasts)=cSPM.con(g).Z; % T
                    column_name{(s-1)*numContrasts+g}=['mean_' SPMname{s} '_' contrastName{g}];
                    column_name{(s-1)*numContrasts+g+numSPM*numContrasts}=['T_' SPMname{s} '_' contrastName{g}];
                end;
            end;
            C=caret_struct('metric','data',data,'column_name',column_name);
            caret_save([surfaceGroupDir '/' 'summary.metric'],C);
        end;
        
    case 'ROI_suit' %_____________Generates suit ROI nii file
        %df1_imana('suit_ROI', 1)
        subj=varargin{1};
        for sn=subj
            SUIT_num= {[3 5],[17 20],[4 7],[19 22]};
            V= spm_vol(fullfile(fileparts(which('suit_reslice' )), 'atlas', 'Cerebellum-SUIT.nii'));
            R_atlas= spm_read_vols(V);
            cd(fullfile(anatomicalDir,subj_name{sn}))
            R=zeros(size(R_atlas));
            for i= 1:4
                for j=1:length(SUIT_num{i})
                    R(R_atlas==SUIT_num{i}(j))=i;
                end;
            end
            V.fname= fullfile('suit', ['ROI_cerebellum.nii']);
            spm_write_vol(V,R);
            %----create the ROI image in subj space saved in the anatomical folder
            refImage= fullfile(baseDir, 'imaging_data',subj_name{sn}, ['rmask_noskull.nii']);
            defMat= fullfile(anatomicalDir,subj_name{sn},'suit',['mc_',subj_name{sn},'_anatomical_snc.mat']);
            source= fullfile(anatomicalDir,subj_name{sn},'suit', ['ROI_cerebellum.nii']);
            suit_reslice_inv(source, defMat,'reference',refImage,'prefix', '');
            movefile(fullfile(anatomicalDir,subj_name{sn},'suit', ['ROI_cerebellum.nii']), fullfile(anatomicalDir,subj_name{sn}, ['ROI_cerebellum.nii']))
        end;
    case 'ROI_define'                  % Defines the ROIs on a single-subject basis
        %df1_imana('ROI_define', 1:4)
        sn= varargin{1};
        
        for s=sn
            fprintf('%d \n',s);
            R=[];
            for h=1:2
                C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},['ROI.paint']));
                caretSubjDir=fullfile(caretDir,['x' subj_name{s}]);
                file=fullfile(glmDir,subj_name{s},'mask.img');
                for i=1:numregions_surf
                    R{i+(h-1)*numregions}.type='surf_nodes';
                    R{i+(h-1)*numregions}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                    R{i+(h-1)*numregions}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                    R{i+(h-1)*numregions}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                    R{i+(h-1)*numregions}.linedef=[5,0,1];
                    R{i+(h-1)*numregions}.image=file;
                    R{i+(h-1)*numregions}.name=[subj_name{s} '_' regname{i} '_' hem{h}];
                    R{i+(h-1)*numregions}.location=find(C.data(:,1)==i);
                end;
                % Cerebellum
                for i= 1:numregions_suit
                    R{i+numregions_surf+(h-1)*numregions}=region('roi_image',fullfile(anatomicalDir,subj_name{s}, ['ROI_cerebellum.nii']), i+(h-1)*numregions_suit);
                    R{i+numregions_surf+(h-1)*numregions}.name=[subj_name{s} '_' regname{i+numregions_surf} '_' hem{h}];
                    R{i+numregions_surf+(h-1)*numregions}.image=file;
                end;
                %                 end;
            end;
            
            R=region_calcregions(R,'exclude',[1 2;11 12],'exclude_thres',0.75);
            cd(regDir);
            save([subj_name{s} '_regions.mat'],'R');
        end;
    case 'ROI_data' %_____________DF______________This gives an ROI structure: relies on the fact that the regions
        % Are properly defined on the mask image, so should not include any
        % NaNs
        %df1_imana('ROI_data', 1:4)
        dataext='8';  % Standard 8 cortical sl1 regions
        sn=varargin{1};
        T=[];
        for s=sn
            fprintf('%s\n',subj_name{s});
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            SI=load('SPM_info.mat');    % Load information
            j=find(SI.regType>0); % Find all the regressors of interest
            
            P={};
            load(fullfile(regDir,[subj_name{s} '_regions.mat']));
            for i=1:length(j)
                P{i}=sprintf('beta_%4.4d.img',j(i));
            end;
            
            % Add a few extra images
            %----task against rest
            P{end+1}='psc_L_motor.nii'; %L_motor
            P{end+1}='psc_R_motor.nii'; %R_motor
            P{end+1}='psc_R_sens.nii'; %R_sens
            P{end+1}='ResMs.img';
            
            E=getrow(SI,j); % Retain only regressors of interest
            
            V=spm_vol(char(P));
            data = region_getdata(V,R);
            
            for i=1:length(R)
                if (~isempty(R{i}))
                    D=[];
                    vec=ones(size(data{i},2),1);
                    if s==17
                        D.beta=data{i}(1:length(j),:)';
                        sdiff = nan(size(D.beta,1),120-size(D.beta,2));
                        D.beta = [D.beta, sdiff];
                    else
                        D.beta=data{i}(1:length(j),:)';
                    end;
                    D.SN=vec*s;
                    D.regNum=vec*i;
                    D.xyz=R{i}.data;
                    for z=1:3 %----task against rest
                        D.psc(:,z)=data{i}(length(j)+z,:)';
                    end;
                    D.ResMs=data{i}(length(j)+4,:)';
                    T=addstruct(T,D);
                end;
            end;
        end;
        regDir=fullfile(baseDir,'RegionOfInterest');
        cd(regDir);
        save(fullfile(regDir,['reg_data_' dataext '.mat']),'-struct','T');
        varargout={T};
    case 'ROI_meanAct'
        reg = [1 2 11 12];
        T = load(fullfile(regDir,'reg_data_8.mat'));
        T = getrow(T,ismember(T.regNum,reg));
        
        S=[];
        for s=1:max(T.SN)
            E=load(fullfile(glmDir, subj_name{s},'SPM_info.mat'));
            for r=reg
                Er = E;
                Tr = getrow(T,T.regNum==r & T.SN==s);
                Er.meanAct = mean(Tr.beta,1)';
                Er = tapply(Er,{'sn','regType','digit','hand','stimType'},{'meanAct','nanmean'});
                Er.regSide = repmat(regSide(r),length(Er.sn),1);
                Er.regType = repmat(regType(r),length(Er.sn),1);
                Er.group   = repmat(subj_group(s),length(Er.sn),1);
                
                S = addstruct(S,Er);
            end;
        end;
        save(fullfile(regDir,'mean_betas.mat'),'-struct','S');
        
    case 'ROI_stats' %____________DF______________caluculates some statistic on the region...
        %D=df1_imana('ROI_stats', 'fcn','ROI_decomposition','selection', 'Acc_rightHand', 'prct', 80)
        %D=df1_imana('ROI_stats', 'fcn','ROI_decomposition','selection', 'none')
        %D=df1_imana('ROI_stats','fcn','accstats_subspace','selection','none','param',{80});save('reg_acc_8.mat','-struct','D');
        %
        selection='none';
        fcn='stats';
        prct=0;             % Cut off for selection criterion
        regions=[1:numregions*2];
        exclusion=1;
        data='beta';
        param={};
        dataext='8';
        
        vararginoptions(varargin,{'regions','data','dataext','selection','prct','standardize','exclusion','fcn','param','prct'});
        
        T=load(fullfile(regDir,['reg_data_' dataext '.mat']));
        
        S=[];
        for s=1:max(T.SN)
            E=load(fullfile(glmDir, subj_name{s},'SPM_info.mat'));
            for r=regions
                indicator=(T.SN==s & T.regNum==r);
                switch (selection)
                    case 'Favrg'
                        indx=find(indicator & T.Favrg>=prctile(T.Favrg(indicator,:),prct));
                    case 'Accavrg'
                        indx=find(indicator & T.accavrg>=prctile(T.accavrg(indicator,:),prct));
                    case 'Acc_rightHand'
                        acc_avrg= mean(T.volAcc(:,2:3),2);
                        indx=find(indicator & acc_avrg>=prctile(acc_avrg(indicator,:),prct));
                    case 'percIpsi'
                        indx=find(indicator & T.percIpsi>0.2);
                    case 'Fcontra_rightHand'
                        F_avrg= mean(T.Ftest(:,2:3),2);
                        indx=find(indicator & F_avrg>=prctile(F_avrg(indicator,:),prct));
                    otherwise
                        indx=find(indicator);
                end;
                
                y=T.(data)(indx,:);
                if (~isempty(indx))
                    D=feval(@df1_imana,fcn,y,E,getrow(T,indx),param{:});
                    fn=fieldnames(D);
                    vec=ones(size(D.(fn{1}),1),1);
                    D.SN=s*vec;
                    D.region=r*vec;
                    D.group=subj_group(s)*vec;
                    D.numvoxels=length(indx)*vec;
                    S=addstruct(S,D);
                end;
            end;
        end;
        Side=[ones(1,numregions) ones(1,numregions)*2];
        Type=[1:numregions 1:numregions];
        
        S.regSide=Side(S.region)';
        S.regType=Type(S.region)';
        
        varargout={S};
    case 'ROI_decomposition' %____DF______________
        y=varargin{1};
        I=varargin{2};
        T=varargin{3};
        %reorder y: first 8 block sensory finger (1:5) & 8 blocks motor fingers (1:5)
        i=1;
        for s=[0 1] %0=> of(motor) 1=> passive mov(sensory)
            for r=unique(I.run)'
                for d=unique(I.digit)'
                    idx(i)= find(I.hand==1 & I.stimType==s & I.run==r & I.digit==d);
                    i=i+1;
                end
            end
        end
        [C u res]= df1_decomp_struct(y(:,idx),I.digit(idx),I.stimType(idx),I.run(idx));
        
        %Calculate spatial kernels on
        if (~isempty(u))
            borders=[0 2 3 4 6 8 10 12 16 20 ];
            %                     % Different estimates
            indx={1,2,1:2,3:7,8:12,3:12,13:28};
            name={'cs','cm','c','fs','fm','f','b','e'};
            %
            %                     %got to mm (coords multiplied with the voxelsize)
            T.xyz= T.xyz;  % ask JD this is WRONG I think I need to add xyz_vox and how shall we do that we non cubic voxel  T.xyz_vox*2.3
            for k=1:8 % Different u's
                if (k<8)
                    U=u(indx{k},:);
                else
                    U=res;
                end;
                [C.(['Corrbin_' name{k}]),C.Distbin]=ss1_spatial_kernel(T.xyz,U',borders);
            end;
            
        end;
        varargout={C};
    case 'ROI_distance'                  % Do extraction of time series to LDA-t values
        % T=df1_imana('ROI_distance',[1:11]);
        % save('reg_distance_raw.mat','-struct','T');
        
        selection='none';
        fcn='stats';
        prct=0;
        regions=[1 2 11 12];
        param={};
        sn=varargin{1};
        T=[];
        
        for s=sn
            glmDirSubj=fullfile([glmDir], subj_name{s});
            cd(fullfile(glmDirSubj));
            load SPM;
            D=load(fullfile(glmDir, subj_name{s},'SPM_info.mat'));
            load(fullfile(regDir,[subj_name{s} '_regions.mat']));
            
            % Get and move the raw data files
            raw_data=SPM.xY.P; % Get regions
            Raw={};
            for i=1:size(raw_data,1)
                raw_data(i,:)=strrep(raw_data(i,:),'\','/');
                [d,name,ext,num]=spm_fileparts(raw_data(i,:));
                Raw{i}=fullfile(baseDir,'imaging_data',subj_name{s},[name ext num]);
            end;
            V=spm_vol(char(Raw));
            
            % Loop over the possile regions
            for r=regions
                
                % Get the data and prewhiten
                Y = region_getdata(V,R{r});  % Data is N x P
                P=size(Y,2);
                beta=mva_prewhiten_beta(Y,SPM);
                
                % Make the two contrast matrices:
                Z1=indicatorMatrix('identity_p',(double(D.hand==0 & D.stimType==0) .* D.digit));
                Z2=indicatorMatrix('identity_p',(double(D.hand==1 & D.stimType==0) .* D.digit));
                Z3=indicatorMatrix('identity_p',(double(D.hand==1 & D.stimType==1) .* D.digit));
                C=indicatorMatrix('allpairs',[1:5]);
                
                % Distance
                S.dist(1,:) = distance_ldc(beta,Z1,C,D.run);
                S.dist(2,:) = distance_ldc(beta,Z2,C,D.run);
                S.dist(3,:) = distance_ldc(beta,Z3,C,D.run);
                
                vec=[1;1;1];
                S.SN=s*vec;
                S.subj=repmat(subj_name(s),length(vec),1);
                S.region=r*vec;
                S.hand=[0;1;1];
                S.stimtype=[0;0;1];
                S.group=subj_group(s)*vec;
                
                T=addstruct(T,S);
                fprintf('%d %d\n',s,r);
            end;
        end;
        
        T.regSide=regSide(T.region)';
        T.regType=regType(T.region)';
        save(fullfile(regDir,'reg_distance_raw.mat'),'-struct','T');
        varargout={T};
    case 'ROI_distanceSplitHalf'                  % Do extraction of time series to LDA-t values
        selection='none';
        fcn='stats';
        prct=0;
        %         regions=[1 2 3 9 10    11 12 13 19 20];
        regions=[1 2 11 12];
        param={};
        sn=varargin{1};
        T=[];
        
        for s=sn
            glmDirSubj=fullfile([glmDir], subj_name{s});
            cd(fullfile(glmDirSubj));
            load SPM;
            D=load(fullfile(glmDir, subj_name{s},'SPM_info.mat'));
            load(fullfile(regDir,[subj_name{s} '_regions.mat']));
            
            % Get and move the raw data files
            raw_data=SPM.xY.P; % Get regions
            Raw={};
            for i=1:size(raw_data,1)
                raw_data(i,:)=strrep(raw_data(i,:),'\','/');
                [d,name,ext,num]=spm_fileparts(raw_data(i,:));
                Raw{i}=fullfile(baseDir,'imaging_data',subj_name{s},[name ext num]);
            end;
            V=spm_vol(char(Raw));
            
            % Loop over the possile regions
            for r=regions
                
                % Get the data and prewhiten
                Y = region_getdata(V,R{r});  % Data is N x P
                P=size(Y,2);
                beta=mva_prewhiten_beta(Y,SPM);
                
                % Make the two contrast matrices:
                Z1=indicatorMatrix('identity_p',(double(D.hand==0 & D.stimType==0) .* D.digit));
                Z2=indicatorMatrix('identity_p',(double(D.hand==1 & D.stimType==0) .* D.digit));
                Z3=indicatorMatrix('identity_p',(double(D.hand==1 & D.stimType==1) .* D.digit));
                C=indicatorMatrix('allpairs',[1:5]);
                
                % indexes for odd and even runs
                odd     = logical(mod(D.run,2));
                even    = logical(1-mod(D.run,2));
                
                % distance for odd runs
                S.dist1(1,:) = distance_ldc(beta(odd,:),Z1(odd,:),C,D.run(odd));
                S.dist1(2,:) = distance_ldc(beta(odd,:),Z2(odd,:),C,D.run(odd));
                S.dist1(3,:) = distance_ldc(beta(odd,:),Z3(odd,:),C,D.run(odd));
                
                % distance for even runs
                S.dist2(1,:) = distance_ldc(beta(even,:),Z1(even,:),C,D.run(even));
                S.dist2(2,:) = distance_ldc(beta(even,:),Z2(even,:),C,D.run(even));
                S.dist2(3,:) = distance_ldc(beta(even,:),Z3(even,:),C,D.run(even));
                
                vec=[1;1;1];
                S.SN=s*vec;
                S.subj=repmat(subj_name(s),length(vec),1);
                S.region=r*vec;
                S.hand=[0;1;1];
                S.stimtype=[0;0;1];
                S.group=subj_group(s)*vec;
                
                T=addstruct(T,S);
                fprintf('%d %d\n',s,r);
            end;
        end;
        
        T.regSide=regSide(T.region)';
        T.regType=regType(T.region)';
        save(fullfile(regDir,'reg_distance_raw_splithalf.mat'),'-struct','T');
        varargout={T};
    case 'ROI_spatVSfunc'                  % Do extraction of time series to LDA-t values
        % T=df1_imana('ROI_distance',[1:11]);
        % save('reg_distance_raw.mat','-struct','T');
        
        selection='none';
        fcn='stats';
        prct=0;
        % regions=[1 2 3 9 10    11 12 13 19 20];
        regions=[1 2];
        param={};
        Act = load(fullfile(regDir,'reg_data_8.mat'));
        sn=varargin{1};
        T=[];
        
        for s=sn
            glmDirSubj=fullfile([glmDir], subj_name{s});
            cd(fullfile(glmDirSubj));
            load SPM;
            D=load(fullfile(glmDir, subj_name{s},'SPM_info.mat'));
            load(fullfile(regDir,[subj_name{s} '_regions.mat']));
            
            % Get and move the raw data files
            raw_data=SPM.xY.P; % Get regions
            Raw={};
            for i=1:size(raw_data,1)
                raw_data(i,:)=strrep(raw_data(i,:),'\','/');
                [d,name,ext,num]=spm_fileparts(raw_data(i,:));
                Raw{i}=fullfile(baseDir,'imaging_data',subj_name{s},[name ext num]);
            end;
            V=spm_vol(char(Raw));
            
            % Loop over the possile regions
            for r=regions
                indx=(Act.SN==s & Act.regNum==r);
                
                % Get the 200 most activated voxels in each region
                % By magnitude of activation
                M=mean(Act.psc(indx,1:2),2)./sqrt(Act.ResMs(indx,:));
                [~,i]=sort(M,1,'descend');
                R{r}.data=R{r}.data(i(1:min(300,length(i))),:);
                
                A       = getrow(Act,indx);
                A       = getrow(A,i(1:min(300,length(i))));
                
                % calculate odd and even beta patterns
                e1   = (double(D.hand==0 & D.stimType==0) .* D.digit);
                e2   = (double(D.hand==1 & D.stimType==0) .* D.digit);
                e3   = (double(D.hand==1 & D.stimType==1) .* D.digit);
                
                e2(e2>0) = e2(e2>0) + 5;
                e3(e3>0) = e3(e3>0) + 10;
                % event   = sum([e1 e2 e3],2);
                event   = sum([e2 e3],2);
                
                Y   = region_getdata(V,R{r});  % Data is N x P
                
                for m=1:4
                    clear beta xyz
                    switch(m)
                        case 1      % use multivariate noise norm using avg condition vector
                            beta = mva_prewhiten_beta(Y,SPM)';
                            C = indicatorMatrix('identity',D.run);
                            beta = beta'-C*pinv(C)*beta';
                            
                            % mean run subtraction
                            % imagesc(indicatorMatrix('identity',D.run))
                            beta = pivottablerow(event,beta,'mean(x,1)','subset',event~=0);
                            xyz  = A.xyz;
                        case 2      % use univariate noise norm using avg condition vector
                            % subtract run mean
                            C = indicatorMatrix('identity',D.run);
                            
                            % univariate noise norm
                            beta = bsxfun(@rdivide,A.beta,A.ResMs);
                            beta = beta'-C*pinv(C)*beta';
                            beta = pivottablerow(event,beta,'mean(x,1)','subset',event~=0);
                            xyz  = A.xyz;
                        case 3      % use multivariate noise norm using full condition vector
                            beta = mva_prewhiten_beta(Y,SPM)';
                            C = indicatorMatrix('identity',D.run);
                            beta = beta'-C*pinv(C)*beta';
                            
                            % mean run subtraction
                            % imagesc(indicatorMatrix('identity',D.run))
                            beta = pivottablerow([D.run event],beta,'mean(x,1)','subset',event~=0);
                            xyz  = A.xyz;
                        case 4      % use univariate noise norm using full condition vector
                            % subtract run mean
                            C = indicatorMatrix('identity',D.run);
                            
                            % univariate noise norm
                            beta = bsxfun(@rdivide,A.beta,A.ResMs);
                            beta = beta'-C*pinv(C)*beta';
                            beta = pivottablerow([D.run event],beta,'mean(x,1)','subset',event~=0);
                            xyz  = A.xyz;
                    end;
                    
                    % estimate metrics
                    fDist = pdist(beta','euclidean');
                    cDist = pdist(beta','correlation');
                    sDist = pdist(xyz,'euclidean');
                    
                    e           = [-Inf 0:1:20 Inf];
                    [cnt, bin]  = histc(sDist',e);
                    fcnt        = pivottable(bin,[],fDist','mean');
                    ccnt        = pivottable(bin,[],cDist','mean');
                    eupdated    = e(unique(bin));
                    if length(eupdated)~=length(fcnt)
                        keyboard;
                    end;
                    
                    vec = ones(length(eupdated),1);
                    S.SN        = s*vec;
                    S.region    = r*vec;
                    S.hand      = 0*vec;
                    S.stimtype  = 0*vec;
                    S.method    = m*vec;
                    S.voxdist   = eupdated';
                    S.funcdist  = fcnt;
                    S.corrdist  = ccnt;
                    S.group     = subj_group(s)*vec;
                    
                    T=addstruct(T,S);
                    fprintf('%d %d %m\n',s,r);
                    disp(size(beta));
                end;
            end;
        end;
        
        T.regSide=regSide(T.region)';
        T.regType=regType(T.region)';
        save(fullfile(regDir,'reg_distance_spat_vs_func.mat'),'-struct','T');
        varargout={T};
    case 'ROI_activation'                  % Do extraction of time series to LDA-t values
        regions=[1 2 11 12];                % left right S1/M1
        sn=varargin{1};
        
        S = [];
        Si = [];
        for s=sn
            glmDirSubj=fullfile([glmDir], subj_name{s});
            cd(fullfile(glmDirSubj));
            load SPM;
            D=load(fullfile(glmDir, subj_name{s},'SPM_info.mat'));
            load(fullfile(regDir,[subj_name{s} '_regions.mat']));
            
            % Get and move the raw data files
            raw_data=SPM.xY.P; % Get regions
            Raw={};
            for i=1:size(raw_data,1)
                raw_data(i,:)=strrep(raw_data(i,:),'\','/');
                [d,name,ext,num]=spm_fileparts(raw_data(i,:));
                Raw{i}=fullfile(baseDir,'imaging_data',subj_name{s},[name ext num]);
            end;
            V=spm_vol(char(Raw));
            
            % Loop over the possile regions
            for r=regions
                % Get the data and prewhiten
                Y   = region_getdata(V,R{r});  % Data is N x P
                beta=mva_prewhiten_beta(Y,SPM);
                
                % Loop over the regions and get estimated data
                v           = ones(length(D.run),1);
                Si.SN       = s*v;
                Si.region   = r*v;
                Si.hand     = D.hand;
                Si.digit    = D.digit;
                Si.stimType = D.stimType;
                Si.act      = mean(beta,2);
                Si.group    = subj_group(s)*v;
                S           = addstruct(S,Si);
                
                fprintf('%d %d\n',s,r);
            end;
        end;
        Side = [1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2];
        S.regSide = Side(S.region)';
        S = getrow(S,S.hand~=S.regSide);
        save(fullfile(regDir,'reg_activations.mat'),'-struct','S');
        varargout={S};
        
    case 'ROI_distance_plot'                  % Plot distance values
        % df1_imana('ROI_distance_plot',1,'ldc');
        color={[0 0 1],[1 0 0],[0 0 1],[1 0 0]};
        style_s={'-','-',':',':'};
        D=load(fullfile(regDir,'reg_distance_raw.mat'));
        regType=varargin{1};
        type = varargin{2};
        HAND=[0 1 1];
        STIMTYPE=[0 0 1];
        plotname={'Left Motor','Right Motor','Right Sensory'};
        
        D.ndist=bsxfun(@rdivide,D.dist,sqrt(sum(D.dist.^2,2)));
        
        set(gcf,'Name',regname{regType})
        D=getrow(D,D.regType==regType);
        maxim=max(max(D.(type)));
        
        for h=1:3
            for s=0:1
                
                subplot(2,3,h+s*3);
                traceplot([1:10],D.(type),'split',[D.group],'leg',{'control','patient'},...
                    'subset',D.regType==regType & D.hand==HAND(h) & D.stimtype==STIMTYPE(h) & D.regSide==s,...
                    'linecolor',color,'patchcolor',color,'linestyle',style_s,'errorfcn','stderr');
                
                set(gca,'YLim',[-maxim/10 maxim*0.8]);
                if(s==1)
                    title(plotname{h});
                end;
                if(h==1 & s==0)
                    ylabel('Left Hemisphere');
                end;
                if(h==1 & s==1)
                    ylabel('Right Hemisphere');
                end;
                
            end;
        end;
    case 'ROI_Gstats'                  % Do extraction of time series to LDA-t values
        % T=df1_imana('ROI_distance',[1:11]);
        % save('reg_distance_raw.mat','-struct','T');
        
        selection='none';
        fcn='stats';
        prct=0;
        % regions=[1 2 3 9 10    11 12 13 19 20];
        % regions=[1 2 3 9 10 11 12 13 19 20];
        regions=[1 2 11 12];
        param={};
        Act = load(fullfile(regDir,'reg_data_8.mat'));
        sn=varargin{1};
        T=[];
        
        for s=sn
            glmDirSubj=fullfile([glmDir], subj_name{s});
            cd(fullfile(glmDirSubj));
            load SPM;
            D=load(fullfile(glmDir, subj_name{s},'SPM_info.mat'));
            load(fullfile(regDir,[subj_name{s} '_regions.mat']));
            
            % Get and move the raw data files
            raw_data=SPM.xY.P; % Get regions
            Raw={};
            for i=1:size(raw_data,1)
                raw_data(i,:)=strrep(raw_data(i,:),'\','/');
                [d,name,ext,num]=spm_fileparts(raw_data(i,:));
                Raw{i}=fullfile(baseDir,'imaging_data',subj_name{s},[name ext num]);
            end;
            V=spm_vol(char(Raw));
            
            
            % Make design matrix
            
            % Loop over the possile regions
            for r=regions
                indx=(Act.SN==s & Act.regNum==r);
                
                % Get the 200 most activated voxels in each region
                % By magnitude of activation
                M=mean(Act.psc(indx,1:2),2)./sqrt(Act.ResMs(indx,:));
                [~,i]=sort(M,1,'descend');
                R{r}.data=R{r}.data(i(1:min(200,length(i))),:);
                
                % Get the data and prewhiten
                Y = region_getdata(V,R{r});  % Data is N x P
                P=size(Y,2);
                beta=mva_prewhiten_beta(Y,SPM);
                
                % Make the Design matrix:
                D.cond=(D.stimType+D.hand)+1;
                D.trialType=(D.cond-1)*5+D.digit;
                Z=indicatorMatrix('identity',D.trialType);
                
                % Generate contrast matrices on the G-matrix
                con=kron([1:3],ones(1,5));
                digit=kron(ones(1,3),[1:5]);
                H={};
                for i=1:3
                    H{end+1}=double(con==i)'*double(con==i);
                end;
                for i=1:3
                    for j=i+1:3
                        H{end+1}=double(con==i)'*double(con==j)+double(con==j)'*double(con==i);
                    end;
                end;
                DI=repmat(digit,15,1)==repmat(digit',1,15);
                for i=1:6
                    H{i+6}=H{i}.*DI;
                end;
                for i=1:12
                    HX(:,i)=H{i}(:);
                end;
                
                % Make the contrast matrix for possible pairwise
                C1=[indicatorMatrix('allpairs',[1:5]) zeros(10,10)];
                C2=[zeros(10,5) indicatorMatrix('allpairs',[1:5]) zeros(10,5)];
                C3=[zeros(10,10) indicatorMatrix('allpairs',[1:5])];
                C=[C1;C2;C3];
                
                [~,G]=crossval_estG(beta,Z,D.run);  % removing negative eigenvalues
                
                b=pinv(HX)*G(:);
                
                S.SN=s;
                S.region=r;
                S.psc=[mean(Act.psc(indx,1)) mean(Act.psc(indx,2)) mean(Act.psc(indx,3))];
                S.avrgdist=b(7:9)';
                S.crossdist=b(10:12)';
                S.G=G(:)';
                S.dist=diag(C*G*C')';
                S.group=subj_group(s);
                
                A=diag(S.avrgdist)+squareform(S.crossdist); % covariance matrix
                Ai=inv(sqrt(diag(diag(A))));                % D=sqrt(diag(C)); Corr = inv(D)*C*inv(D)
                CorrA=Ai*A*Ai;
                S.corr=1-squareform((1-CorrA).*(1-eye(3))); % 1-eye(3) is for numerical precision of diagonal not being exactly zero
                
                T=addstruct(T,S);
                fprintf('%d %d\n',s,r);
            end;
        end;
        
        T.regSide=regSide(T.region)';
        T.regType=regType(T.region)';
        save(fullfile(regDir,'reg_Gstats_raw.mat'),'-struct','T');
        varargout={T};
    case 'ROI_Gstat_plot'
        reg=[1 2 11 12];
        group=[1 2];
        T=load(fullfile(regDir,'reg_Gstats_raw.mat'));
        
        vararginoptions(varargin,{'reg','group'});
        T=getrow(T,ismember(T.group,group));
        
        figure(1);
        subplot(2,1,1)
        barplot(T.region,[T.avrgdist(:,1:2) T.crossdist(:,1)],'subset',ismember(T.region,reg));
        title('left-right motor');
        subplot(2,1,2)
        barplot(T.region,[T.avrgdist(:,2:3) T.crossdist(:,3)],'subset',ismember(T.region,reg));
        title('right sens-motor');
        
        figure(2);
        numregions=length(reg);
        for r=1:numregions
            ridx=find(T.region==reg(r));
            for i=1:3
                subplot(numregions,3,(r-1)*3+i);
                indx=[1:10]+(i-1)*10;
                ndist=bsxfun(@rdivide,T.dist(ridx,indx),sqrt(sum(T.dist(ridx,indx).^2,2)));
                A=1-pdist(ndist,'correlation');
                imagesc(squareform(mean(ssqrt(ndist))));
                title(sprintf('Rel: %2.3f',mean(A)));
            end;
        end;
        
        figure(3);
        barplot(T.region,T.corr(:,3),'subset',ismember(T.region,reg),'split',T.group);
        title('corr');
        
    case 'ROI_psc_plot'                  % Plot distance values
        % df1_imana('ROI_distance_plot',1,'ldc');
        D=load(fullfile(regDir,'reg_distance_raw.mat'));
        regType=varargin{1};
        
        HAND=[0 1 1];
        STIMTYPE=[0 0 1];
        plotname={'Left Motor','Right Motor','Right Sensory'};
        
        D=getrow(D,D.regType==regType);
        
        myboxplot([D.hand D.stimtype],D.psc,'split',[D.regSide D.group]);
    case 'distance_compare'         % Compare motor distances across studies
        
        color={[0 0 1],[1 0 0],[0 1 0]};
        normalize=0;
        var='dist';
        reg=[1 2];
        vararginoptions(varargin,{'normalize','var','reg'});
        
        field_oth={'SN','hand','region','regSide','regType'};
        
        
        TD=load(fullfile(regDir,'reg_distance_raw.mat'));
        TD.hand=TD.hand+1; % Make hands 1,2: DO THIS IN GENERAL??
        TD.regSide=TD.regSide+1;
        TD=getrow(TD,TD.hand~=TD.regSide & TD.stimtype==0 & ismember(TD.regType,reg));
        
        T1=load(fullfile('/Users/joern/Projects/FingerPattern/tendigit1/RegionOfInterest','reg_distance_raw.mat'));
        T1=getrow(T1,T1.hand~=T1.regSide & T1.method==6 & ismember(T1.regType,reg));
        
        
        % wrangle the data into a new format: In general it may be useful
        % to have this type of structure
        TDa.ndist=bsxfun(@rdivide,TD.dist,sqrt(sum(TD.dist.^2,2)))';TDa.ndist=TDa.ndist(:);
        T1a.ndist=bsxfun(@rdivide,T1.dist,sqrt(sum(T1.dist.^2,2)))';T1a.ndist=T1a.ndist(:);
        
        TDa.dist=TD.dist';TDa.dist=TDa.dist(:);
        T1a.dist=T1.dist';T1a.dist=T1a.dist(:);
        
        for i=1:length(field_oth);
            TDa.(field_oth{i})=kron(TD.(field_oth{i}),ones(10,1));
            T1a.(field_oth{i})=kron(T1.(field_oth{i}),ones(10,1));
        end;
        TDa.pair=kron(ones(length(TD.SN),1),[1:10]');
        T1a.pair=kron(ones(length(T1.SN),1),[1:10]');
        
        T1a.group=ones(length(T1a.SN),1)*1; % Control from tendigit1
        TDa.group=subj_group(TDa.SN)+1; % 2:Pianists 3: Dystonics
        
        T=addstruct(TDa,T1a);
        
        for i=1:2
            subplot(1,2,i);
            lineplot(T.pair,T.(var),'split',[T.group],...
                'catcol',1,...
                'leg',{'TD1','musicians','dystonia'},'style_thickline','subset',T.hand==i);
        end;
        
        set(gca,'XTickLabel',fingerPairs);
        anovaMixed(T.ndist,T.SN+T.group*20,'within',[T.pair T.hand],{'pair','hand'},'between',T.group,{'group'});
        
        varargout={T};
    case 'fit_spatial_figure' % This is Figure X
        %df1_imana('fit_spatial_figure', D)
        T=varargin{1};
        reg= varargin{2};
        eff={'f'};
        %eff={'c','f','b','e'};
        lcolor= [ 0 0 0; 0 0 1; 1 0 0; 0 1 0];
        c=1;
        D=[]
        for reg=reg %=unique(T.region)'
            A=[];
            for e=1:length(eff)
                for sn=unique(T.SN)'
                    i=find(T.SN==sn & T.region==reg);
                    C.db=T.Distbin(i,:);
                    C.r=T.(['Corrbin_' eff{e}])(i,:);
                    C.SN= sn;
                    C.group= subj_group(sn);
                    % With new exponential fit:
                    C.(['SigmaE_' eff{e}])=fminsearch(@exp_error,2,[],C.db,C.r);
                    C.(['FWHME_' eff{e}])=2*sqrt(log(2))*C.(['SigmaE_' eff{e}]);
                    A=addstruct(A,C);
                end;
                %subplot(2,2,e);
                traceplot(mean(A.db),A.r,'errorfcn','stderr', 'leg', 'auto', 'linecolor',lcolor(c,:)); %'split',[D.group],
                c=c+1;
                hold on
                D=addstruct(D,A);
                %waitforbuttonpress
            end;
            %figure
            %keyboard;
        end;
        varargout={D};
    case 'force_scanner'
        sn=varargin{1};
        cd(fullfile(behaviourDir,subj_name{sn}))
        D=dload(['DF1_', subj_name{sn}, '.dat']);
        barplot(D.hand, D.maxForce, 'split', [D.stimType D.digit], 'subset', [D.announce==0& D.digit~=0], 'facecolor',[0.5 0.5 0.5]);
        set(gcf,'PaperPosition',[1 1 6 2]);
        wysiwyg;
        saveas(gcf, 'maxForce', 'jpg')
    case 'searchlight_make_nii'        % Makes an nii file for a location in the image
        sn=varargin{1};
        coord=varargin{2};
        wdir = fullfile(baseDir,'glm_firstlevel_1', subj_name{sn});
        S=load(fullfile(wdir,'vol_roi_80vox.mat'));
        V=spm_vol('mask.img');
        [i,j,k]=spmj_affine_transform(coord(1),coord(2),coord(3),inv(V.mat));
        I=surfing_coords2linvoxelidxs(coord',V);
        i=find(S.voxel==I);
        if (isempty(i))
            error('voxel not found');
        end;
        X=zeros(V.dim);
        X(S.LI{i})=1;
        V.fname=fullfile(wdir,sprintf('searchlight_%d_%d_%d.nii',coord(1),coord(2),coord(3)));
        spm_write_vol(V,X);
        
    case 'plot_fingerpatterns_ci'               % Makes a plot of all finger patterns
        sn=varargin{1};
        h=varargin{2};
        atlas=2;
        vararginoptions({varargin{3:end}},{'atlas'});
        
        groupDir=[caretDir filesep atlasname{atlas}  filesep hemName{h} ];
        cd(groupDir);
        border=fullfile(caretDir,atlasname{atlas},hemName{h},['CS.border']);
        switch(h)
            case 2
                coord='rh.FLAT.coord';
                topo='rh.CUT.topo';
                data='rh.surface_shape';
                xlims=[-10 5];
                ylims=[-5 10];
            case 1
                coord='lh.FLAT.coord';
                topo='lh.CUT.topo';
                data='lh.surface_shape';
                xlims=[-5 10];
                ylims=[-5 10];
        end;
        B=caret_load(border);
        
        data=fullfile(caretDir,[atlasA{atlas} subj_name{sn}],hemName{h},[subj_name{sn} '_finger.metric']);
        sshape=fullfile(caretDir,atlasname{atlas},hemName{h},[hem{h} '.surface_shape']);
        
        subplot(4,5,16);
        M=caret_plotflatmap('col',2,'data',sshape,'col',1,'border',B.Border,'topo',topo,'coord',coord,'xlims',xlims,'ylims',ylims);
        
        for i=1:15
            subplot(4,5,i);
            [M,d]=caret_plotflatmap('M',M,'col',i,'data',data,'cscale',[-6 6],...
                'border',B.Border,'topo',topo,'coord',coord);
            maxT(i)=max(d(:));
        end;
        mm=max(maxT);
        for i=1:15
            subplot(4,5,i);
            caxis([-mm/4 mm]);
        end;
        set(gcf,'PaperPosition',[1 1 9 5]);
        wysiwyg;
    case 'plot_fingerpatterns_sulculdepth' 		% make a plot of finger maps for each individual, left hem, passive condition
        sn=varargin{1};
        h = 1;
        atlas=2;
        vararginoptions({varargin{3:end}},{'atlas'});
        
        groupDir=[caretDir filesep atlasname{atlas}  filesep hemName{h} ];
        cd(groupDir);
        border=fullfile(caretDir,atlasname{atlas},hemName{h},['CS.border']);
        switch(h)
            case 2
                coord='rh.FLAT.coord';
                topo='rh.CUT.topo';
                data='rh.surface_shape';
                xlims=[-10 5];
                ylims=[-5 10];
            case 1
                coord='lh.FLAT.coord';
                topo='lh.CUT.topo';
                data='lh.surface_shape';
                xlims=[-5 10];
                ylims=[-5 10];
        end;
        B=caret_load(border);
        
        data=fullfile(caretDir,[atlasA{atlas} subj_name{sn}],hemName{h},[subj_name{sn} '_finger.metric']);
        sshape=fullfile(caretDir,atlasname{atlas},hemName{h},[hem{h} '.surface_shape']);
        
        M=caret_plotflatmap('col',2,'data',sshape,'col',1,'border',B.Border,'topo',topo,'coord',coord,'xlims',xlims,'ylims',ylims);
        plt.subplot(4,5,1);
        plt.set('xtick',[],'ytick',[],'ratio','square');
        
    case 'plot_fingerpatterns_passive' 		% make a plot of finger maps for each individual, left hem, passive condition
        sn=varargin{1};
        offset=varargin{2};
        h = 1;
        atlas=2;
        vararginoptions({varargin{3:end}},{'atlas'});
        
        groupDir=[caretDir filesep atlasname{atlas}  filesep hemName{h} ];
        cd(groupDir);
        border=fullfile(caretDir,atlasname{atlas},hemName{h},['CS.border']);
        switch(h)
            case 2
                coord='rh.FLAT.coord';
                topo='rh.CUT.topo';
                data='rh.surface_shape';
                xlims=[-10 5];
                ylims=[-5 10];
            case 1
                coord='lh.FLAT.coord';
                topo='lh.CUT.topo';
                data='lh.surface_shape';
                xlims=[-5 10];
                ylims=[-5 10];
        end;
        B=caret_load(border);
        
        data=fullfile(caretDir,[atlasA{atlas} subj_name{sn}],hemName{h},[subj_name{sn} '_finger.metric']);
        sshape=fullfile(caretDir,atlasname{atlas},hemName{h},[hem{h} '.surface_shape']);
        
        plt.subplot(5,5,1);
        M=caret_plotflatmap('col',2,'data',sshape,'col',1,'border',B.Border,'topo',topo,'coord',coord,'xlims',xlims,'ylims',ylims);
        plt.set('xtick',[],'ytick',[],'ratio','square');
        
        for i=1:5
            plt.subplot(5,5,5*offset + i);
            [M,d]=caret_plotflatmap('M',M,'col',10+i,'data',data,'cscale',[-6 12],...
                'border',B.Border,'topo',topo,'coord',coord);
            maxT(i)=max(d(:));
        end;
        mm=max(maxT);
        for i=1:5
            plt.subplot(5,5,5*offset + i);
            caxis([-mm/4 mm]);
            
            plt.set('xtick',[],'ytick',[],'ratio','square');
        end;
    case 'FIG_group_FingerPatterns_Dystonics'         % Finger patterns for the dystonic group
        atlas       = 2;
        hand        = 2;    % right hand
        hemi        = 1;    % left hemisphere
        reg         = 1;    % S1
        metric      = 3;    % max
        condition   = 2;
        vararginoptions(varargin,{'hemi','condition','hand'});
        
        cData   = zeros(163842,5,length(subj_name));
        handIdx = 5*(hand-1)+[1:5];
        x = zeros(5,length(subj_name));
        y = zeros(5,length(subj_name));
        
        % 1. Getting activation patterns for all subjects/fingers
        D = load(fullfile(regDir,'spatial_CoG.mat'));
        for sn=1:length(subj_name)
            fprintf('SN %d, REG %d, HEM %d\n',sn,reg,hemi);
            
            dataFile    = fullfile(caretDir,[atlasA{atlas} subj_name{sn}],hemName{hemi},[subj_name{sn} '_finger.metric']);
            data        = caret_load(dataFile);
            data        = data.data(:,handIdx);
            
            % storing data
            cData(:,:,sn) = data;
            
            % get CoG measurements
            Di = getrow(D,strcmp(D.subjname,subj_name{sn}) & D.metric==metric & D.region==reg & D.hand==hand & D.condition==condition & D.regSide==hemi);
            x(:,sn) = Di.x;
            y(:,sn) = Di.y;
        end;
        
        % 2. Getting means for finger patterns across groups
        fp_cont.data    = nanmean(cData(:,:,subj_group==1),3);
        fp_dyst.data    = nanmean(cData(:,:,subj_group==2),3);
        
        % 3. Getting group surface
        groupDir    = [caretDir filesep atlasname{atlas}  filesep hemName{hemi} ];
        border      = fullfile(caretDir,atlasname{atlas},hemName{hemi},['CS.border']);
        paint       = fullfile(caretDir,atlasname{atlas},hemName{hemi},['ROI.paint']);
        shape       = fullfile(caretDir,atlasname{atlas},hemName{hemi},[hem{hemi} '.surface_shape']);
        cd(groupDir);
        
        switch(hemi)
            case 2
                coord='rh.FLAT.coord';
                topo='rh.CUT.topo';
                xlims=[-10 5];
                ylims=[-5 10];
            case 1
                coord='lh.FLAT.coord';
                topo='lh.CUT.topo';
                xlims=[-5 10];
                ylims=[-5 10];
        end;
        B       = caret_load(border);
        P       = caret_load(paint);
        sshape  = caret_load(shape);
        
        % % TO DELETE!!!!
        % idxROI  = P.data==reg;
        % fp_cont.data(~idxROI,:)=nan;
        
        figure;
        subplot(5,5,1);
        colormap(gray);
        M       = caret_plotflatmap('col',2,'data',sshape,'col',1,'border',B.Border,...
            'topo',topo,'coord',coord,'xlims',xlims,'ylims',ylims);
        axis square;
        set(gca,'XTickLabel',[],'YTickLabel',[]);
        
        % 4. Plotting finger patterns
        for i=1:length(handIdx)
            plt.subplot(5,5,5*offset + i);
            [M,d]=caret_plotflatmap('M',M,'col',i,'data',fp_dyst,'cscale',[-6 12],...
                'border',B.Border,'topo',topo,'coord',coord);
            maxT(i+5)=max(d(:));
            axis square;
            set(gca,'XTickLabel',[],'YTickLabel',[]);
            hold on;
            plot(x(i,subj_group==2),y(i,subj_group==2),'+r','MarkerSize',4);
            hold off;
            
        end;
        
        for i=1:5
            plt.subplot(5,5,5*offset + i);
            caxis([-mm/4 mm]);
            
            plt.set('xtick',[],'ytick',[],'ratio','square');
        end;
        
        % saving results
        % save_figure(gcf,fullfile(figureDir,sprintf('%s.pdf',what)),'style','brain_3row');
        
    case 'FIG_fingermaps_passive'  % Maps for Figure 4 in paper
        df1_imana('plot_fingerpatterns_passive',6,1);
        df1_imana('plot_fingerpatterns_passive',8,2);
        df1_imana('plot_fingerpatterns_passive',14,3);
        df1_imana('plot_fingerpatterns_passive',15,4);
        f = gcf;
        % saveas(f,'Fingermaps_passive.tif');
        
    case 'plot_fingerpatterns_ipsicontra'               % Makes a plot of all finger patterns
        sn=varargin{1};
        atlas=2;
        vararginoptions({varargin{3:end}},{'atlas'});
        
        loc = [1 1 1 1];
        dig = [7 9 2 4];
        sc  = {[-6 12],[-6 12],[-6 12],[-6 12]};
        vec = [];
        S = [];
        
        for i=1:length(loc)
            h = loc(i);
            groupDir=[caretDir filesep atlasname{atlas}  filesep hemName{h} ];
            cd(groupDir);
            border=fullfile(caretDir,atlasname{atlas},hemName{h},['CS.border']);
            switch(h)
                case 2
                    coord='rh.FLAT.coord';
                    topo='rh.CUT.topo';
                    data='rh.surface_shape';
                    % xlims=[-10 20];
                    % ylims=[-15 30];
                    xlims=[-5 15];
                    ylims=[-5 15];
                case 1
                    coord='lh.FLAT.coord';
                    topo='lh.CUT.topo';
                    data='lh.surface_shape';
                    % xlims=[-20 10];
                    % ylims=[-15 30];
                    xlims=[-15 5];
                    ylims=[-5 15];
            end;
            B=caret_load(border);
            
            
            data=fullfile(caretDir,[atlasA{atlas} subj_name{sn}],hemName{h},[subj_name{sn} '_finger.metric']);
            sshape=fullfile(caretDir,atlasname{atlas},hemName{h},[hem{h} '.surface_shape']);
            
            subplot(2,4,i);
            M=caret_plotflatmap('col',2,'data',sshape,'col',1,'border',B.Border,'topo',topo,'coord',coord,'xlims',xlims,'ylims',ylims);
            
            subplot(2,4,4+i);
            % x           = caret_load(data);
            % x           = x.data(:,dig(i));
            % x           = (x-nanmean(x))/nanstd(x);
            
            [M,d]=caret_plotflatmap('M',M,'col',dig(i),'data',data,'cscale',sc{i},...
                'border',B.Border,'topo',topo,'coord',coord);
            maxT(i)=max(d(:));
            
            x           = caret_load(data);
            Si.vec      = x.data(:,dig(i));
            Si.dig      = i * ones(length(Si.vec),1);
            S           = addstruct(S,Si);
            
            vec(:,i) = Si.vec;
        end;
        
        mm=max(maxT);
        for i=1:4
            subplot(2,4,4+i);
            caxis([-mm/2 mm]);
        end;
        
        subplot(2,4,1:4)
        vec = vec(~any(isnan(vec),2),:);
        vec = vec(~any(isnan(vec),2),:);
        plot(vec);
        disp(corr(vec,'type','pearson'));
        varargout = {S};
    case 'SPAT_calcmetric'  % This calculates the various spatial statistics
        dat = varargin{1};
        x = varargin{2};
        y = varargin{3};
        m = varargin{4};
        k_const     = [0 0.05 0.1 0.2 0.4 0.8 1.6 3.2 10];
        k       = k_const(m);
        
        switch(m)
            case 1          % weighted mean
                dat(dat<0)=0;
                w       = dat/sum(dat);
                X       = sum(w.*x);
                Y       = sum(w.*y);
            case {2,3,4,5,6,7,8}   % Softmax 0-1
                w       = exp(k*dat);
                % w(dat<0)= 0;
                w       = w/sum(w);
                X       = sum(w.*x);
                Y       = sum(w.*y);
            case 9
                [~,idx] = max(dat); 
                X       = x(idx);
                Y       = y(idx);
        end;
        varargout={X,Y,k};
    case 'SPAT_metrics2d'
        atlas       = 2;
        roiType     = [1 2]; % Use subsection of M1 and S1 
        isplot      = 1;
        symbol      = {'kx','go','ko','ko','ko','ko','ko','ro','k*'};
        vararginoptions(varargin,{'atlas','roiType','isplot'});
        
        S=[];
        hand  = [ones(1,5) 2*ones(1,10)];
        digit = repmat(1:5,1,3);
        condition  = [ones(1,10) 2*ones(1,5)];
        
        for h=1:2
            switch(h)
                case 2
                    coord='rh.FLAT.coord';
                    topo='rh.CUT.topo';
                    xlims=[-15 15];
                    ylims=[-10 20];
                case 1
                    coord='lh.FLAT.coord';
                    topo='lh.CUT.topo';
                    xlims=[-15 15];
                    ylims=[-10 20];
            end
            groupDir    = [caretDir filesep atlasname{atlas}  filesep hemName{h} ];
            border      = fullfile(caretDir,atlasname{atlas},hemName{h},['CS.border']);
            paint       = fullfile(caretDir,atlasname{atlas},hemName{h},['ROI.paint']); % Use subsections of M1 / S1 
            cd(groupDir);
            C = caret_load(coord);
            B=caret_load(border);
            P=caret_load(paint);
            
            % Plot surface shape for the area
            if (isplot)
                subplot(4,5,16);
                sshape=fullfile(caretDir,atlasname{atlas},hemName{h},[hem{h} '.surface_shape']);
                M=caret_plotflatmap('data',sshape,'col',1,'border',B.Border,'topo',topo,'coord',coord,'xlims',xlims,'ylims',ylims);
            end
            
            for sn=1:length(subj_name)
                dataFile=fullfile(caretDir,[atlasA{atlas} subj_name{sn}],hemName{h},[subj_name{sn} '_finger.metric']);
                data    = caret_load(dataFile);
                for i=1:15
                    % Plot the flatmap activation for all finger maps
                    if (isplot)
                        subplot(4,5,i);
                        caret_plotflatmap('M',M,'col',i,'data',data,'cscale',[-6 12],...
                            'border',B.Border,'topo',topo,'coord',coord);
                    end
                end
                for r=roiType
                    fprintf('SN %d, REG %d, HEM %d\n',sn,r,h);
                    % Determine the vertices in the region of question
                    idxROI  = P.data==r & ~isnan(sum(data.data,2));
                    pNans = sum(P.data==r & isnan(sum(data.data,2)))/sum(P.data==r);  % Proportion of Vertices with Nans
                    
                    for i=1:15
                        % Determine the activations and the coordinates
                        dat = data.data(idxROI,i);
                        group = floor((i-1)/5); 
                        datN = dat-mean(data.data(idxROI,group*5+1:group*5+5),2)
                        x = C.data(idxROI,1);
                        y = C.data(idxROI,2);
                        % Now try out different method to determine the
                        % COG:
                        for m=1:10
                            if (m==10) 
                                [X,Y,k]=df1_imana('SPAT_calcmetric',datN,x,y,1);
                            else 
                                [X,Y,k]=df1_imana('SPAT_calcmetric',dat,x,y,m);
                            end; 
                            % Add the symbol
                            if (isplot)
                                subplot(4,5,i);
                                hold on;
                                plot(X,Y,symbol{m});
                                hold off;
                            end
                            
                            % Saving data for each subject/hemisphere/hand
                            Si.x        = X;
                            Si.y        = Y;
                            Si.regSide  = h;
                            Si.region   = r;
                            Si.hand     = hand(i);
                            Si.digit    = digit(i);
                            Si.condition= condition(i);
                            Si.metric   = m;
                            Si.k        = k;
                            Si.pNans    = pNans;
                            Si.sn       = sn;
                            Si.group    = subj_group(sn);
                            Si.subjname = subj_name(sn);
                            S           = addstruct(S,Si);
                        end
                    end
                end
                if (isplot)
                    keyboard;
                end;
            end
        end
        
        varargout = {S};
        save(fullfile(regDir,'spatial_CoG.mat'),'-struct','S');
    case 'SPAT_metrics2d_splithalf'
        atlas       = 2;
        roiType     = [1 2];
        isplot      = 1;
        symbol      = {'kx','go','ko','ko','ko','ko','ko','ro','k*'};
        vararginoptions(varargin,{'atlas','roiType','isplot'});
        
        S=[];
        hand  = [ones(1,5) 2*ones(1,10)];
        digit = repmat(1:5,1,3);
        condition  = [ones(1,10) 2*ones(1,5)];
        
        type = {'odd','even'};
        
        for h=1:2
            switch(h)
                case 2
                    coord='rh.FLAT.coord';
                    topo='rh.CUT.topo';
                    xlims=[-15 15];
                    ylims=[-10 20];
                case 1
                    coord='lh.FLAT.coord';
                    topo='lh.CUT.topo';
                    xlims=[-15 15];
                    ylims=[-10 20];
            end
            groupDir    = [caretDir filesep atlasname{atlas}  filesep hemName{h} ];
            border      = fullfile(caretDir,atlasname{atlas},hemName{h},['CS.border']);
            paint       = fullfile(caretDir,atlasname{atlas},hemName{h},['ROI.paint']);
            cd(groupDir);
            C = caret_load(coord);
            B=caret_load(border);
            P=caret_load(paint);
            
            % Plot surface shape for the area
            if (isplot)
                subplot(4,5,16);
                sshape=fullfile(caretDir,atlasname{atlas},hemName{h},[hem{h} '.surface_shape']);
                M=caret_plotflatmap('col',2,'data',sshape,'col',1,'border',B.Border,'topo',topo,'coord',coord,'xlims',xlims,'ylims',ylims);
            end
            
            for sn=1:length(subj_name)
                for ty=1:2
                    dataFile=fullfile(caretDir,[atlasA{atlas} subj_name{sn}],hemName{h},[subj_name{sn} '_finger_' type{ty} '.metric']);
                    data    = caret_load(dataFile);
                    for i=1:15
                        % Plot the flatmap activation for all finger maps
                        if (isplot)
                            subplot(4,5,i);
                            caret_plotflatmap('M',M,'col',i,'data',data,'cscale',[-6 12],...
                                'border',B.Border,'topo',topo,'coord',coord);
                        end
                    end
                    for r=roiType
                        fprintf('Type %d, SN %d, REG %d, HEM %d\n',ty,sn,r,h);
                        % Determine the vertices in the region of question
                        idxROI  = P.data==r & ~isnan(sum(data.data,2));
                        pNans = sum(P.data==r & isnan(sum(data.data,2)))/sum(P.data==r);  % Proportion of Vertices with Nans
                        
                        for i=1:15
                            % Determine the activations and the coordinates
                            dat = data.data(idxROI,i);
                            group = floor((i-1)/5); 
                            datN = dat-mean(data.data(idxROI,group*5+1:group*5+5),2);
                            x = C.data(idxROI,1);
                            y = C.data(idxROI,2);
                            % Now try out different method to determine the
                            % COG:
                            for m=1:10
                                if (m==10) 
                                    [X,Y,k]=df1_imana('SPAT_calcmetric',datN,x,y,1);
                                else 
                                    [X,Y,k]=df1_imana('SPAT_calcmetric',dat,x,y,m);
                                end; 
                                % Add the symbol
                                if (isplot)
                                    subplot(4,5,i);
                                    hold on;
                                    plot(X,Y,symbol{m});
                                    hold off;
                                end
                                
                                % Saving data for each subject/hemisphere/hand
                                Si.x        = X;
                                Si.y        = Y;
                                Si.regSide  = h;
                                Si.region   = r;
                                Si.hand     = hand(i);
                                Si.digit    = digit(i);
                                Si.condition= condition(i);
                                Si.metric   = m;
                                Si.k        = k;
                                Si.type     = ty;
                                Si.pNans    = pNans;
                                Si.sn       = sn;
                                Si.group    = subj_group(sn);
                                Si.subjname = subj_name(sn);
                                S           = addstruct(S,Si);
                            end
                        end
                        if (isplot)
                            keyboard;
                        end;
                    end
                end
            end
        end
        
        varargout = {S};
        save(fullfile(regDir,['spatial_CoG_splithalf.mat']),'-struct','S');
    case 'SPAT_metrics2d_plot'
        atlas       = 2;
        roiType     = 3;  % BA 3 
        h           = 1;
        sn          = 7;
        condition   = 2;
        hand        = 2;
        symbol      = {'kx','go','ko','ko','ko','ko','ko','ro','k*'};
        vararginoptions(varargin,{'atlas','roiType','h','condition','hand','sn'});
        
        figure(1);
        set(gcf,'PaperPosition',[2 2 12 2]);
        wysiwyg;
        S=[];
        handV  = [ones(1,5) 2*ones(1,10)];
        digitV = repmat(1:5,1,3);
        conditionV  = [ones(1,10) 2*ones(1,5)];
        switch(h)
            case 2
                coord='rh.FLAT.coord';
                topo='rh.CUT.topo';
                xlims=[-10 5];
                ylims=[-5 10];
            case 1
                coord='lh.FLAT.coord';
                topo='lh.CUT.topo';
                xlims=[-5 10];
                ylims=[-5 10];
        end
        groupDir    = [caretDir filesep atlasname{atlas}  filesep hemName{h} ];
        border      = fullfile(caretDir,atlasname{atlas},hemName{h},['CS.border']);
        paint       = fullfile(caretDir,atlasname{atlas},hemName{h},['ROIsm.paint']);
        cd(groupDir);
        C = caret_load(coord);
        B=caret_load(border);
        P=caret_load(paint);
        
        dataFile=fullfile(caretDir,[atlasA{atlas} subj_name{sn}],hemName{h},[subj_name{sn} '_finger.metric']);
        data    = caret_load(dataFile);
        axis equal;
        
        M = caret_plotflatmap('coord',coord,'topo',topo,'col',1,'data',data,'cscale',[-6 12],...
            'border',B.Border,'xlims',xlims,'ylims',ylims);
        % Determine the vertices in the region of question
        idxROI  = P.data==roiType & C.data(:,1)>xlims(1) & C.data(:,1)<xlims(2)& C.data(:,2)>ylims(1) & C.data(:,2)<ylims(2) & ~isnan(sum(data.data,2));
        idxROI  = P.data==roiType & ~isnan(sum(data.data,2));
        
        for digit=1:5
            subplot(1,5,digit);
            i = find(handV==hand & digitV==digit & conditionV==condition);
            caret_plotflatmap('M',M,'col',i,'data',data,'cscale',[-6 12],'border',B.Border,'coord',coord,'topo',topo);
            set(gca,'XTick',[],'YTick',[]);
            axis equal; 
            
            
            % Determine the activations and the coordinates
            dat = data.data(idxROI,i);
            x = C.data(idxROI,1);
            y = C.data(idxROI,2);
            % Now try out different method to determine the
            
            for m=1:9
                [X,Y,k]=df1_imana('SPAT_calcmetric',dat,x,y,m);
                % Add the symbol
                hold on;
                plot(X,Y,symbol{m});
                hold off;
            end
        end;
    case 'SPAT_tuningFunctionDistribution'      % shows the distribution of tuning functions within an individuals
        atlas       = 2;
        hemi        = 1;     % 1. left/2. right
        reg         = 2;     % S1 & M1 only
        condition   = 1;     % only looking at active motor condition
        vararginoptions(varargin,{'atlas','hemi','reg'});
        
        S               = [];
        handLab         = [ones(1,5) 2*ones(1,10)];
        digitLab        = repmat(1:5,1,3);
        conditionLab    = [ones(1,10) 2*ones(1,5)];
        dataIdx         = find(conditionLab==condition & hemi~=handLab);  % contra-lat hand, active motor condition
        for sn=1:length(subj_name)
            fprintf('SN %d, HEM %d\n',sn,hemi);
            
            groupDir    = [caretDir filesep atlasname{atlas}  filesep hemName{hemi} ];
            border      = fullfile(caretDir,atlasname{atlas},hemName{hemi},['CS.border']);
            paint       = fullfile(caretDir,atlasname{atlas},hemName{hemi},['ROI.paint']);
            cd(groupDir);
            
            switch(hemi)
                case 2
                    coord='rh.FLAT.coord';
                    topo='rh.CUT.topo';
                    xlims=[-10 20];
                    ylims=[-10 20];
                case 1
                    coord='lh.FLAT.coord';
                    topo='lh.CUT.topo';
                    xlims=[-20 10];
                    ylims=[-10 20];
            end;
            B=caret_load(border);
            P=caret_load(paint);
            
            dataFile=fullfile(caretDir,[atlasA{atlas} subj_name{sn}],hemName{hemi},[subj_name{sn} '_finger.metric']);
            sshape=fullfile(caretDir,atlasname{atlas},hemName{hemi},[hem{hemi} '.surface_shape']);
            
            % selecting only vertices that belong to the defined region
            % of interest
            data                    = caret_load(dataFile);
            idxROI                  = ismember(P.data,reg);
            data.data(~idxROI,:)    = nan;
            [~,tf]                  = max(data.data(:,dataIdx),[],2);
            tf(~idxROI,:)           = nan;
            
            M=caret_plotflatmap('col',2,'data',sshape,'col',1,'border',B.Border,'topo',topo,'coord',coord,'xlims',xlims,'ylims',ylims);
            caret_plotflatmap('M',M,'data',tf,'border',B.Border,'topo',topo,'coord',coord);
            colormap(jet);
            h = colorbar;
            % set(h,'Ticks',0:5);
            
            keyboard;
        end;
    case 'SPAT_checkSomatotopy'
        condition   = 2;    % sensory
        hand        = 2;    % right hand
        group       = [1];    % healthy participants
        reg         = 1;    % S1
        vararginoptions(varargin,{'reg','hand','condition','metric','group'});
        grpLabel    = {'controls','dystonic'};
        
        % Run Manova on healthy group to check somatotopy for the different
        % softmax values
        D = load(fullfile(regDir,'spatial_CoG.mat'));
        D = getrow(D,ismember(D.group,group) & D.condition==condition & ...
            D.region==reg & D.hand==hand & D.hand~=D.regSide);
        
        for m=6
            fprintf('metric=%d\n',m);
            [f1,r]=pivottable([D.sn D.digit],[],D.x,'mean','subset',D.metric==m);
            [f2,r]=pivottable([D.sn D.digit],[],D.y,'mean','subset',D.metric==m);
            Y       = [f1 f2];
            R       = r(:,1);
            F       = r(:,2);
            T=MANOVArp(R,F,Y);
            fprintf('\n\n');
        end;
        
        % Visualizing somatotopic gradient for softmax k=0.2
        CAT.markercolor={[0 0 0.7] [0.5 0 0.5] [0.7 0 0] [0.5 0.5 0],[0 0.7 0]};
        CAT.markerfill={[0 0 0.7] [0.5 0 0.5] [0.7 0 0] [0.5 0.5 0],[0 0.7 0]};
        CAT.markertype='o';
        CAT.markersize=10;
        CAT.errorcolor={[0 0 0.7] [0.5 0 0.5] [0.7 0 0] [0.5 0.5 0],[0 0.7 0]};
        CAT.errorwidth=2;
        
        h = plt.figure;
        [x,y]=xyplot(D.x,D.y,D.digit,'subset',D.metric==4,'split',D.digit,...
            'CAT',CAT,'leg','none');%,'leglocation','southwest','leg','auto');
        % plt.set(h,'ylim',[2 6],'xlim',[3 5]);
        % plt.labels('anterior -> posterior','ventral -> dorsal');
        
        % saving results
        % save_figure(gcf,fullfile(figureDir,sprintf('%s_%s.pdf',what,grpLabel{group})),'style','brain_1col');
        varargout={T}; 
    case 'SPAT_estimateDistances'           % estimate distances based on the calculate x-y coordinates
        metric      = 3;    % cog, without negative values
        vararginoptions(varargin,{'hand','condition','metric'});
        
        D = load(fullfile(regDir,'spatial_CoG.mat'));
        % D = getrow(D,D.metric==metric & D.hand~=D.regSide); % hand is always contralateral and metric is softmax
        D = getrow(D,D.hand~=D.regSide); % hand is always contralateral
        
        % Looping over subjects/regions/conditions/hand
        S = [];
        for sn = unique(D.sn)'
            Ds = getrow(D,D.sn==sn);
            
            for reg = unique(Ds.region)'
                Dr = getrow(Ds,Ds.region==reg);
                
                for cond = unique(Dr.condition)'
                    Dc = getrow(Dr,Dr.condition==cond);
                    
                    for h = unique(Dc.hand)'
                        Dh  = getrow(Dc,Dc.hand==h);
                        
                        for m = unique(Dh.metric)'
                            Dm  = getrow(Dh,Dh.metric==m);
                            x   = pivottablerow(Dm.digit,[Dm.x Dm.y],'nanmean(x,1)');
                            d   = pdist(x,'euclidean');
                            
                            % Saving results
                            Si.sn           = sn;
                            Si.metric       = m;
                            Si.region       = reg;
                            Si.condition    = cond;
                            Si.hand         = h;
                            Si.group        = mean(Dm.group);
                            Si.dist         = d;
                            Si.k            = mean(Dm.k);
                            S               = addstruct(S,Si);
                            
                            fprintf('SN: %d, Reg: %d, Cond: %d, Hand: %d, Metric: %d\n',sn,reg,cond,h,m);
                        end;
                    end;
                end;
            end;
        end;
        varargout = {S};
        save(fullfile(regDir,'spatial_distances.mat'),'-struct','S');
    case 'SPAT_estimateDistances_splithalf'           % estimate distances based on the calculate x-y coordinates
        D = load(fullfile(regDir,'spatial_CoG_splithalf.mat'));
        D = getrow(D,D.hand~=D.regSide); % hand is always contralateral
        
        % Looping over subjects/regions/conditions/hand
        S = [];
        for ty=1:2
            Dt = getrow(D,D.type==ty);
            
            for sn = unique(D.sn)'
                Ds = getrow(Dt,Dt.sn==sn);
                
                for reg = unique(Ds.region)'
                    Dr = getrow(Ds,Ds.region==reg);
                    
                    for cond = unique(Dr.condition)'
                        Dc = getrow(Dr,Dr.condition==cond);
                        
                        for h = unique(Dc.hand)'
                            Dh  = getrow(Dc,Dc.hand==h);
                            
                            for m = unique(Dh.metric)'
                                Dm     	= getrow(Dh,Dh.metric==m);
                                [x,dig] = pivottablerow(Dm.digit,[Dm.x Dm.y],'nanmean(x,1)');
                                d       = pdist(x,'euclidean');
                                
                                % Saving results
                                Si.sn           = sn;
                                Si.metric       = m;
                                Si.region       = reg;
                                Si.condition    = cond;
                                Si.hand         = h;
                                Si.group        = mean(Dm.group);
                                Si.type         = ty;
                                Si.k            = mean(Dm.k);
                                Si.dist         = d;
                                Si.x            = x(:,1)';
                                Si.y            = x(:,2)';
                                Si.digit        = dig';
                                S               = addstruct(S,Si);
                                
                                fprintf('Type: %d, SN: %d, Reg: %d, Cond: %d, Hand: %d, Metric: %d\n',ty,sn,reg,cond,h,m);
                            end;
                        end;
                    end;
                end;
            end;
        end;
        
        varargout = {S};
        save(fullfile(regDir,'spatial_distances_splithalf.mat'),'-struct','S');
    case 'SPAT_compareDistances'            % compare distances across
        % Compare distances only for sensory condition only
        D = load(fullfile(regDir,'spatial_distances.mat'));
        D = getrow(D,D.condition==2 & D.metric==3);
        
        
        s   = style_sheet(sty_grp,'leg',{'controls','dystonic'},'leglocation','northeast','errorcolor','match','fillcolor','match');
        h   = make_figure;
        myboxplot(D.region,mean(D.dist,2),'split',D.group,'style_tukey','plotall',0','xtickoff','CAT',s(1).CAT,s(1).PLOT{:});
        set_graphics(h,'ylabel','avg. distance (mm)','xticklabel',{'S1','S1','M1','M1'},'ytick',0:0.2:1.4,'ylim',[0 1.3]);
        
        % saving results
        save_figure(gcf,fullfile(figureDir,sprintf('%s.pdf',what)),'style','brain_1col');
    case 'SPAT_compareDistancesStats'       % compare distances across STATS
        reg = 1;
        vararginoptions(varargin,{'reg'});
        
        % Compare distances only for sensory condition only
        D = load(fullfile(regDir,'spatial_distances.mat'));
        D = getrow(D,D.condition==2 & D.metric~=2 & D.region==reg);
        
        for m=unique(D.metric)'
            Di = getrow(D,D.metric==m);
            disp(['Metric ' num2str(m)]);
            MANOVA1(Di.group,Di.dist);
        end;
    case 'FIG_group_FingerPatterns'         % Finger patterns for the two groups (healthy vs dystonic)
        atlas       = 2;
        hand        = 2;    % right hand
        hemi        = 1;    % left hemisphere
        reg         = 1;    % S1
        metric      = 3;    % max
        condition   = 2;
        vararginoptions(varargin,{'hemi','condition','hand'});
        
        cData   = zeros(163842,5,length(subj_name));
        handIdx = 5*(hand-1)+[1:5];
        x = zeros(5,length(subj_name));
        y = zeros(5,length(subj_name));
        
        % 1. Getting activation patterns for all subjects/fingers
        D = load(fullfile(regDir,'spatial_CoG.mat'));
        for sn=1:length(subj_name)
            fprintf('SN %d, REG %d, HEM %d\n',sn,reg,hemi);
            
            dataFile    = fullfile(caretDir,[atlasA{atlas} subj_name{sn}],hemName{hemi},[subj_name{sn} '_finger.metric']);
            data        = caret_load(dataFile);
            data        = data.data(:,handIdx);
            
            % storing data
            cData(:,:,sn) = data;
            
            % get CoG measurements
            Di = getrow(D,strcmp(D.subjname,subj_name{sn}) & D.metric==metric & D.region==reg & D.hand==hand & D.condition==condition & D.regSide==hemi);
            x(:,sn) = Di.x;
            y(:,sn) = Di.y;
        end;
        
        % 2. Getting means for finger patterns across groups
        fp_cont.data    = nanmean(cData(:,:,subj_group==1),3);
        fp_dyst.data    = nanmean(cData(:,:,subj_group==2),3);
        
        % 3. Getting group surface
        groupDir    = [caretDir filesep atlasname{atlas}  filesep hemName{hemi} ];
        border      = fullfile(caretDir,atlasname{atlas},hemName{hemi},['CS.border']);
        paint       = fullfile(caretDir,atlasname{atlas},hemName{hemi},['ROI.paint']);
        shape       = fullfile(caretDir,atlasname{atlas},hemName{hemi},[hem{hemi} '.surface_shape']);
        cd(groupDir);
        
        switch(hemi)
            case 2
                coord='rh.FLAT.coord';
                topo='rh.CUT.topo';
                xlims=[-10 5];
                ylims=[-5 10];
            case 1
                coord='lh.FLAT.coord';
                topo='lh.CUT.topo';
                xlims=[-5 10];
                ylims=[-5 10];
        end;
        B       = caret_load(border);
        P       = caret_load(paint);
        sshape  = caret_load(shape);
        
        % % TO DELETE!!!!
        % idxROI  = P.data==reg;
        fp_cont.data(~idxROI,:)=nan;
        
        figure;
        subplot(3,5,11);
        colormap(gray);
        M       = caret_plotflatmap('col',2,'data',sshape,'col',1,'border',B.Border,...
            'topo',topo,'coord',coord,'xlims',xlims,'ylims',ylims);
        axis square;
        set(gca,'XTickLabel',[],'YTickLabel',[]);
        
        % 4. Plotting finger patterns
        for i=1:length(handIdx)
            subplot(3,5,i);
            colormap('default');
            [M,d]=caret_plotflatmap('M',M,'col',i,'data',fp_cont,'cscale',[-6 12],...
                'border',B.Border,'topo',topo,'coord',coord);
            maxT(i)=max(d(:));
            axis square;
            set(gca,'XTickLabel',[],'YTickLabel',[]);
            hold on;
            plot(x(i,subj_group==1),y(i,subj_group==1),'+r','MarkerSize',4);
            hold off;
            
            subplot(3,5,i+5);
            [M,d]=caret_plotflatmap('M',M,'col',i,'data',fp_dyst,'cscale',[-6 12],...
                'border',B.Border,'topo',topo,'coord',coord);
            maxT(i+5)=max(d(:));
            axis square;
            set(gca,'XTickLabel',[],'YTickLabel',[]);
            hold on;
            plot(x(i,subj_group==2),y(i,subj_group==2),'+r','MarkerSize',4);
            hold off;
            
        end;
        
        % Re-scaling
        mm=max(maxT);
        for i=1:10
            subplot(3,5,i);
            caxis([-mm/2 mm]);
        end;
        
        % saving results
        % save_figure(gcf,fullfile(figureDir,sprintf('%s.pdf',what)),'style','brain_3row');
        
    case 'FIG_group_FingerPatternsThumb'         % Finger patterns for the two groups (healthy vs dystonic)
        atlas       = 2;
        hand        = 2;    % right hand
        hemi        = 1;    % left hemisphere
        reg         = 1;    % S1
        metric      = 3;    % max
        condition   = 2;
        sn1         = 6;
        sn2         = 12;
        vararginoptions(varargin,{'hemi','condition','hand','sn1','sn2'});
        
        cData   = zeros(163842,5,length(subj_name));
        handIdx = 5*(hand-1)+[1:5];
        x = zeros(5,length(subj_name));
        y = zeros(5,length(subj_name));
        xmax = zeros(5,length(subj_name));
        ymax = zeros(5,length(subj_name));
        
        % 1. Getting activation patterns for all subjects/fingers
        D = load(fullfile(regDir,'spatial_CoG.mat'));
        for sn=1:length(subj_name)
            fprintf('SN %d, REG %d, HEM %d\n',sn,reg,hemi);
            
            dataFile    = fullfile(caretDir,[atlasA{atlas} subj_name{sn}],hemName{hemi},[subj_name{sn} '_finger.metric']);
            data        = caret_load(dataFile);
            data        = data.data(:,handIdx);
            
            % storing data
            cData(:,:,sn) = data;
            
            % get CoG measurements
            Di = getrow(D,strcmp(D.subjname,subj_name{sn}) & D.metric==metric & D.region==reg & D.hand==hand & D.condition==condition & D.regSide==hemi);
            x(:,sn) = Di.x;
            y(:,sn) = Di.y;
            
            Di = getrow(D,strcmp(D.subjname,subj_name{sn}) & D.metric==1 & D.region==reg & D.hand==hand & D.condition==condition & D.regSide==hemi);
            xmax(:,sn) = Di.x;
            ymax(:,sn) = Di.y;
        end;
        
        % 2. Getting means for finger patterns across groups
        % sn1 = 5; sn2 = 6;
        fp_cont.data    = cData(:,:,sn1);
        fp_dyst.data    = cData(:,:,sn2);
        
        % 3. Getting group surface
        groupDir    = [caretDir filesep atlasname{atlas}  filesep hemName{hemi} ];
        border      = fullfile(caretDir,atlasname{atlas},hemName{hemi},['CS.border']);
        paint       = fullfile(caretDir,atlasname{atlas},hemName{hemi},['ROI.paint']);
        shape       = fullfile(caretDir,atlasname{atlas},hemName{hemi},[hem{hemi} '.surface_shape']);
        cd(groupDir);
        
        switch(hemi)
            case 2
                coord='rh.FLAT.coord';
                topo='rh.CUT.topo';
                xlims=[-10 5];
                ylims=[-5 10];
            case 1
                coord='lh.FLAT.coord';
                topo='lh.CUT.topo';
                xlims=[-5 10];
                ylims=[-5 10];
        end;
        B       = caret_load(border);
        P       = caret_load(paint);
        sshape  = caret_load(shape);
        
        % % TO DELETE!!!!
        % idxROI  = P.data==reg;
        % fp_cont.data(~idxROI,:)=nan;
        
        figure;
        subplot(3,5,11);
        colormap(gray);
        M       = caret_plotflatmap('col',2,'data',sshape,'col',1,'border',B.Border,...
            'topo',topo,'coord',coord,'xlims',xlims,'ylims',ylims);
        axis square;
        set(gca,'XTickLabel',[],'YTickLabel',[]);
        
        % 4. Plotting finger patterns
        % for i=1:length(handIdx)
        for i=1
            subplot(121);
            colormap('default');
            [M,d]=caret_plotflatmap('M',M,'col',i,'data',fp_cont,...
                'border',B.Border,'topo',topo,'coord',coord);
            maxT(i)=max(d(:));
            axis square;
            set(gca,'XTickLabel',[],'YTickLabel',[]);
            hold on;
            plot(x(i,sn1),y(i,sn1),'+r','MarkerSize',15);
            plot(xmax(i,sn1),ymax(i,sn1),'or','MarkerSize',15);
            hold off;
            
            subplot(122);
            [M,d]=caret_plotflatmap('M',M,'col',i,'data',fp_dyst,...
                'border',B.Border,'topo',topo,'coord',coord);
            maxT(i+5)=max(d(:));
            axis square;
            set(gca,'XTickLabel',[],'YTickLabel',[]);
            hold on;
            plot(x(i,sn2),y(i,sn2),'+r','MarkerSize',15);
            plot(xmax(i,sn2),ymax(i,sn2),'or','MarkerSize',15);
            hold off;
            
        end;
        
        % Re-scaling
        % mm=max(maxT);
        % for i=1:10
        %     subplot(3,5,i);
        %     caxis([-mm/2 mm]);
        % end;
        %
        % saving results
        plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'1x2');
    case 'FIG_group_FingerPatternsSomatotopy'         % Finger patterns for the two groups (healthy vs dystonic)
        atlas       = 2;
        hand        = 2;    % right hand
        hemi        = 1;    % left hemisphere
        reg         = 1;    % S1
        metric      = 3;    % max
        condition   = 2;
        vararginoptions(varargin,{'hemi'});
        
        cData   = zeros(163842,5,length(subj_name));
        handIdx = 5*(hand-1)+[1:5];
        x = zeros(5,length(subj_name));
        y = zeros(5,length(subj_name));
        
        % 1. Getting activation patterns for all subjects/fingers
        D = load(fullfile(regDir,'spatial_CoG.mat'));
        D = getrow(D,D.metric==metric & D.region==reg & D.hand==hand & D.condition==condition & D.regSide==hemi);
        for sn=1:length(subj_name)
            fprintf('SN %d, REG %d, HEM %d\n',sn,reg,hemi);
            
            dataFile    = fullfile(caretDir,[atlasA{atlas} subj_name{sn}],hemName{hemi},[subj_name{sn} '_finger.metric']);
            data        = caret_load(dataFile);
            data        = data.data(:,handIdx);
            
            % storing data
            cData(:,:,sn) = data;
            
            % get CoG measurements
            Di = getrow(D,strcmp(D.subjname,subj_name{sn}) & D.metric==metric & D.region==reg & D.hand==hand & D.condition==condition & D.regSide==hemi);
            x(:,sn) = Di.x;
            y(:,sn) = Di.y;
        end;
        
        % 3. Getting group surface
        groupDir    = [caretDir filesep atlasname{atlas}  filesep hemName{hemi} ];
        border      = fullfile(caretDir,atlasname{atlas},hemName{hemi},['CS.border']);
        paint       = fullfile(caretDir,atlasname{atlas},hemName{hemi},['ROI.paint']);
        shape       = fullfile(caretDir,atlasname{atlas},hemName{hemi},[hem{hemi} '.surface_shape']);
        cd(groupDir);
        
        switch(hemi)
            case 2
                coord='rh.FLAT.coord';
                topo='rh.CUT.topo';
                xlims=[-10 5];
                ylims=[-5 10];
            case 1
                coord='lh.FLAT.coord';
                topo='lh.CUT.topo';
                xlims=[-5 10];
                ylims=[-5 10];
        end;
        B       = caret_load(border);
        P       = caret_load(paint);
        sshape  = caret_load(shape);
        
        % Visualizing somatotopic gradient
        CAT.markercolor={[0 0 0.7] [0.5 0 0.5] [0.7 0 0] [0.5 0.5 0],[0 0.7 0]};
        CAT.markerfill={[0 0 0.7] [0.5 0 0.5] [0.7 0 0] [0.5 0.5 0],[0 0.7 0]};
        CAT.markertype='o';
        CAT.markersize=10;
        CAT.errorcolor={[0 0 0.7] [0.5 0 0.5] [0.7 0 0] [0.5 0.5 0],[0 0.7 0]};
        CAT.errorbars='none';
        
        plt.figure;
        for i=1:2
            subplot(1,2,i);
            [M,d]       = caret_plotflatmap('col',2,'data',sshape,'col',1,'border',B.Border,...
                'topo',topo,'coord',coord,'xlims',xlims,'ylims',ylims);
            colormap(gray);
            axis square;
            set(gca,'XTickLabel',[],'YTickLabel',[]);
        end;
        
        % Re-scaling
        mm = 1.2;
        for i=1:2
            subplot(1,2,i);
            caxis([-mm/2 mm]);
            
            hold on;
            [x,y]=xyplot(D.x,D.y,D.digit,'subset',D.metric==3,'split',D.digit,'subset',D.group==i,...
                'CAT',CAT,'leg','none');
            hold off;
        end;
        
        % saving results
        plt.figure(gcf,fullfile(figureDir,sprintf('%s.pdf',what)),'style','1x2');
    case 'FIG_group_TuningFunctions'        % Tuning functions for each voxel for the two groups (healthy vs dystonic)
        atlas       = 2;
        hand        = 2;    % right hand
        hemi        = 1;    % left hemisphere
        reg         = 1;    % S1
        vararginoptions(varargin,{'hemi'});
        
        % 1. Getting group surface
        groupDir    = [caretDir filesep atlasname{atlas}  filesep hemName{hemi} ];
        border      = fullfile(caretDir,atlasname{atlas},hemName{hemi},['CS.border']);
        paint       = fullfile(caretDir,atlasname{atlas},hemName{hemi},['ROI.paint']);
        shape       = fullfile(caretDir,atlasname{atlas},hemName{hemi},[hem{hemi} '.surface_shape']);
        cd(groupDir);
        
        switch(hemi)
            case 2
                coord='rh.FLAT.coord';
                topo='rh.CUT.topo';
                xlims=[-10 5];
                ylims=[-5 10];
            case 1
                coord='lh.FLAT.coord';
                topo='lh.CUT.topo';
                xlims=[-5 10];
                ylims=[-5 10];
        end;
        B       = caret_load(border);
        P       = caret_load(paint);
        sshape  = caret_load(shape);
        
        % 2. Getting activation patterns for all subjects/fingers
        cData   = zeros(163842,5,length(subj_name));
        handIdx = 5*(hand-1)+[1:5];
        x = zeros(5,length(subj_name));
        y = zeros(5,length(subj_name));
        
        S = [];
        for sn=1:length(subj_name)
            fprintf('SN %d, REG %d, HEM %d\n',sn,reg,hemi);
            
            dataFile    = fullfile(caretDir,[atlasA{atlas} subj_name{sn}],hemName{hemi},[subj_name{sn} '_finger.metric']);
            data        = caret_load(dataFile);
            data        = data.data(:,handIdx);
            
            % storing data
            cData(:,:,sn) = data;
            
            % pick up ROI data
            idxROI  = P.data==reg;
            data    = data(idxROI,:);
            
            [~,tf]      = nanmax(data,[],2);
            Si.sn       = repmat(sn,5,1);
            Si.group    = repmat(subj_group(sn),5,1);
            Si.digit    = [1 2 3 4 5]';
            Si.count    = pivottable(ones(length(tf),1),tf,tf,'length')';
            S           = addstruct(S,Si);
        end;
        
        % 3. Getting means for finger patterns across groups
        %       - get winner take all tuning function
        fp_cont.data    = nanmean(cData(:,:,subj_group==1),3);
        fp_dyst.data    = nanmean(cData(:,:,subj_group==2),3);
        
        [~,fp_cont.data] = nanmax(fp_cont.data,[],2);
        [~,fp_dyst.data] = nanmax(fp_dyst.data,[],2);
        
        % 4. Plotting
        figure;
        subplot(1,3,1);
        M       = caret_plotflatmap('col',2,'data',sshape,'col',1,'border',B.Border,...
            'topo',topo,'coord',coord,'xlims',xlims,'ylims',ylims);
        axis square;
        set(gca,'XTickLabel',[],'YTickLabel',[]);
        
        % 4. Plotting finger patterns
        subplot(1,3,2);
        [M,d]=caret_plotflatmap('M',M,'col',1,'data',fp_cont,'cscale',[1 5],...
            'border',B.Border,'topo',topo,'coord',coord,'interpol','mode');
        axis square;
        set(gca,'XTickLabel',[],'YTickLabel',[]);
        
        subplot(1,3,3);
        [M,d]=caret_plotflatmap('M',M,'col',1,'data',fp_dyst,'cscale',[1 5],...
            'border',B.Border,'topo',topo,'coord',coord,'interpol','mode');
        axis square;
        set(gca,'XTickLabel',[],'YTickLabel',[]);
        
        map = [0 0 0.7; 0.5 0 0.5; 0.7 0 0; 0.5 0.5 0; 0 0.7 0];
        colormap(map);
        
        % saving results
        save_figure(gcf,fullfile(figureDir,sprintf('%s.pdf',what)),'style','brain_1row');
    case 'FIG_group_Distance'               % distances between groups
        D = load(fullfile(regDir,'spatial_distances.mat'));
        D = getrow(D,D.condition==2 & D.hand==2); % sensory condition, S1, right hand
        
        barplot(D.region,mean(D.dist,2),'split',D.group,...
            'leg',{'musician','dystonic'},'leglocation','northeast');
        title('Right hand (symptomatic)');
        ylabel('interdigit CoG distance (mm)');
        set(gca,'XTickLabel',{'S1','S1','M1','M1'},'fontsize',14);
        
        % mean distance in S1
        fprintf('S1\n--\n');
        F = pivottable(D.sn,D.group,D.dist,'mean(x,2)','subset',D.region==1);
        ttest(F(:,1),F(:,2),2,'independent');
        fprintf('\n\n');
        
        % mean distance in M1
        fprintf('M1\n--\n');
        F = pivottable(D.sn,D.group,D.dist,'mean(x,2)','subset',D.region==2);
        ttest(F(:,1),F(:,2),2,'independent');
    case 'FIG_dist_Rel'               % reliability of CoG/RSA distances across groups
        D = load(fullfile(regDir,'spatial_distances.mat'));
        D = getrow(D,D.condition==1 & D.hand==2 & ismember(D.region,[1 2]) & D.group==1); % motor condition, S1, right hand
        
        T = load(fullfile(regDir,'reg_distance_raw.mat'));
        T = getrow(T,T.stimtype==0 & ismember(T.region,[1 2]) & T.hand~=T.regSide & T.group==1); % motor condition, S1, right hand
        
        % Loop over each subj to calculate reliability
        S = [];
        for reg=[1 2]
            Ti = getrow(T,T.region==reg);
            Di = getrow(D,D.region==reg);
            
            for sn=unique(T.SN)'
                Si.sn   = sn;
                Si.reg  = reg;
                
                % CoG distances
                idx         = Di.sn==sn;
                Si.r        = mean(corr(Di.dist(idx,:)',Di.dist(~idx,:)'));
                Si.metric   = 1;
                S           = addstruct(S,Si);
                
                % cortical distances
                idx         = Ti.SN==sn;
                Si.r        = mean(corr(Ti.dist(idx,:)',Ti.dist(~idx,:)'));
                Si.metric   = 2;
                S           = addstruct(S,Si);
            end;
        end;
        
        % Figure
        barplot(S.reg,S.r,'split',S.metric,...
            'leg',{'CoG','RSA'},'leglocation','northwest');
        title('Right hand (symptomatic)');
        ylabel('Pearsons r');
        set(gca,'XTickLabel',{'S1','S1','M1','M1'},'fontsize',14);
        
        % mean distance in S1
        fprintf('S1\n--\n');
        F = pivottable(S.sn,S.metric,S.r,'mean(x,2)','subset',S.reg==1);
        ttest(F(:,1),F(:,2),2,'independent');
        fprintf('\n\n');
        
        % mean distance in M1
        fprintf('M1\n--\n');
        F = pivottable(S.sn,S.metric,S.r,'mean(x,2)','subset',S.reg==2);
        ttest(F(:,1),F(:,2),2,'independent');
    case 'FIG_meanDistance'                 % mean distance across S1 & M1'
        D = load(fullfile(regDir,'reg_distance_Raw.mat'));
        D = rmfield(D,'subj');
        D = getrow(D,D.stimtype==1 & ismember(D.region,[1 2]) & D.hand~=D.regSide);
        
        s   = style_sheet(sty_grp,'leg',{'controls','dystonic'},'leglocation','northoutside','errorcolor','match','fillcolor','match');
        h   = make_figure;
        barplot([D.region],mean(D.dist,2),'split',D.group,'CAT',s(1).CAT,s(1).PLOT{:});
        set_graphics(h,'ylabel','avg. distance (a.u.)','xticklabel',{'S1','S1','M1','M1'});
        
        % saving results
        save_figure(gcf,fullfile(figureDir,sprintf('%s.pdf',what)),'style','brain_2col');
    case 'FIG_patternStructureS1'           % avg. pattern structure in controls
        % Compare distances only for sensory condition only
        D   = load(fullfile(regDir,'reg_distance_raw.mat'));
        D   = rmfield(D,'subj');
        D1  = getrow(D,ismember(D.region,[1]) & D.stimtype==1 & D.hand~=D.regSide);
        D2  = getrow(D,ismember(D.region,[2]) & D.stimtype==1 & D.hand~=D.regSide);
        
        % make image of mean distance pattern
        x   = D1.dist(D1.group,:);
        h   = make_figure;
        imagesc_rectangle(squareform(mean(x,1)),'YDir','reverse','MAP',colormap(hot),'scale',[0 1.5]);
        set_graphics(h,'xtick',[],'ytick',[]);
        
        % saving results
        save_figure(gcf,fullfile(figureDir,sprintf('%s.pdf',what)),'style','brain_1col');
    case 'FIG_patternStructureGroup'        % avg. pattern structure across groups
        % Compare distances only for sensory condition only
        D   = load(fullfile(regDir,'reg_distance_raw.mat'));
        D   = rmfield(D,'subj');
        D1  = getrow(D,ismember(D.region,[1]) & D.stimtype==1 & D.hand~=D.regSide);
        D2  = getrow(D,ismember(D.region,[2]) & D.stimtype==1 & D.hand~=D.regSide);
        
        h   = make_figure;
        subplot(221);
        naturalglove_analyze('MDS_plot',abs(D1.dist(D1.group==1,:)));
        set_graphics(gca,'ylim',[-0.6 0.6],'xlim',[-0.58 1.2],'title','S1 - controls',...
            'ytick',[],'xtick',[]);
        subplot(222);
        naturalglove_analyze('MDS_plot',abs(D1.dist(D1.group==2,:)));
        set_graphics(gca,'ylim',[-0.6 0.6],'xlim',[-0.58 1.2],'title','S1 - dystonics',...
            'ytick',[],'xtick',[]);
        subplot(223);
        naturalglove_analyze('MDS_plot',abs(D2.dist(D2.group==1,:)));
        set_graphics(gca,'ylim',[-0.3 0.3],'xlim',[-0.55 0.7],'title','M1 - controls',...
            'ytick',[],'xtick',[]);
        subplot(224);
        naturalglove_analyze('MDS_plot',abs(D2.dist(D2.group==2,:)));
        set_graphics(gca,'ylim',[-0.3 0.3],'xlim',[-0.55 0.7],'title','M1 - dystonics',...
            'ytick',[],'xtick',[]);
        
        
        % saving results
        save_figure(gcf,fullfile(figureDir,sprintf('%s.pdf',what)),'style','brain_2row');
    case 'FIG_meanDistanceMotor'          % mean distance across regions and conditions
        D   = load(fullfile(regDir,'reg_distance_Raw.mat'));
        D   = rmfield(D,'subj');
        D1  = getrow(D,D.stimtype==0 & ismember(D.region,[1 2]) & D.hand~=D.regSide);
        D2  = getrow(D,D.stimtype==0 & ismember(D.region,[1 2]) & D.hand==D.regSide);
        
        s   = style_sheet(sty_grp,'leg',{'controls','dystonic'},'leglocation','northoutside','errorcolor','match','fillcolor','match');
        h   = make_figure;
        subplot(121);
        barplot(D1.region,mean(D1.dist,2),'split',D1.group,'CAT',s(1).CAT,s(1).PLOT{:});
        title('contralateral');
        subplot(122);
        barplot(D2.region,mean(D2.dist,2),'split',D2.group,'CAT',s(1).CAT,s(1).PLOT{:});
        title('ipsilateral');
        set_graphics(h,'ylabel','avg. distance (a.u.)','xticklabel',{'S1','S1','M1','M1'},...
            'match','ylim','ytick',[0:0.1:0.3]);
        
        % saving results
        save_figure(gcf,fullfile(figureDir,sprintf('%s.pdf',what)),'style','brain_2col');
    case 'SPAT_REL_betweenSubj'  % corr distance measures between subjects
        reg      =1;
        metric   =2;
        group    =1;
        cond     =2;
        vararginoptions(varargin,{'reg','metric','group','cond'})
        
        % 1. Get subset of data
        D = load(fullfile(regDir,'spatial_distances.mat'));
        d = getrow(D,D.region==reg & D.metric==metric & D.condition==cond & D.group==group);
        
        d = getrow(d,d.hand==2);
        
        % 2. Get correlation between subjects
        for h=unique(d.hand)'
            i = 1;
            r = zeros(length(d.sn),1);
            for s=unique(d.sn)'
                d1      = d.dist(d.sn==s & d.hand==h,:);
                d2      = d.dist(d.sn~=s & d.hand==h,:);
                r(i)    = nanmean(corr(d1',d2'))';
                i       = i+1;
            end;
            
            m=nanmean(fisherz(r));
            SE=stderr(fisherz(r));
            fprintf('%1.2f (%1.2f - %1.2f)\n',fisherinv(m),fisherinv(m-1.96*SE),fisherinv(m+1.96*SE));
        end;
        varargout = {r};
    case 'TUNING_REL_betweenSubj'  % corr distance measures between subjects
        reg     = [1 2];
        stimType    = 1;
        group       = 2;
        vararginoptions(varargin,{'reg','stimType','group'})
        
        
        % 1. Get subset of data
        D = load(fullfile(regDir,'reg_distance_raw.mat'));
        % D = rmfield(D,'subj');
        D = getrow(D,ismember(D.region,reg) & D.stimtype==stimType & ...
            D.group==group & D.hand~=D.regSide);
        
        % 2. Get one distance structure for each subject
        T = tapply(D,{'SN'},{'dist','mean(x,1)'});
        
        % 3. Get correlation between subjects
        i = 1;
        r = zeros(length(T.SN),1);
        for s=unique(T.SN)'
            d1      = T.dist(T.SN==s,:);
            d2      = T.dist(T.SN~=s,:);
            r(i)    = mean(corr(d1',d2'));
            i       = i+1;
        end;
        
        % 4. Display intersubject reliabilities
        fprintf('Mean: %2.3f, SE: %2.3f\n',mean(r),stderr(r))
        varargout = {r};
    case 'FIG_groupRelPassive'             % reliability differences between spatial and tuning metrics
        metric = 3;
        reg = [1 2];
        grp = [1];
        
        S = [];
        for r=1:length(reg)
            for g=1:length(grp)
                r1 = df1_imana('SPAT_REL_betweenSubj','group',g,'reg',r,'metric',metric,'cond',2);
                r2 = df1_imana('TUNING_REL_betweenSubj','group',g,'reg',r,'stimType',1);
                
                v           = ones(length(r1),1);
                Si.sn       = [1:length(r1)]';
                Si.reg      = r*v;
                Si.group    = g*v;
                Si.type     = 1*v;
                Si.r        = r1;
                S           = addstruct(S,Si);
                
                Si.type     = 2*v;
                Si.r        = r2;
                S           = addstruct(S,Si);
                clear Si
            end;
        end;
        
        s = style.custom({'blue','black'});
        
        plt.subplot(221);
        plt.bar(S.type,S.r,'subset',S.group==1 & S.reg==1,'split',S.type,'leg','none','style',s);
        plt.subplot(222);
        plt.bar([S.type],S.r,'subset',S.group==1 & S.reg==2,'split',S.type,'leg','none','style',s);
        plt.set(gcf,'xticklabel',{'spatial','RDM','spatial','RDM'},...
            'ytick',0:0.2:1,'match','ylim');
        
        plt.labels([],'Pearsons r','S1','D',221);
        plt.labels([],'Pearsons r','M1',[],222);
        
        plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'1x2');
    case 'FIG_groupRelActive'             % reliability differences between spatial and tuning metrics
        metric = 3;
        reg = [1 2];
        grp = [1];
        
        S = [];
        for r=1:length(reg)
            for g=1:length(grp)
                r1 = df1_imana('SPAT_REL_betweenSubj','group',g,'reg',r,'metric',metric,'cond',1);
                r2 = df1_imana('TUNING_REL_betweenSubj','group',g,'reg',r,'stimType',0);
                
                v           = ones(length(r1),1);
                Si.sn       = [1:length(r1)]';
                Si.reg      = r*v;
                Si.group    = g*v;
                Si.type     = 1*v;
                Si.r        = r1;
                S           = addstruct(S,Si);
                
                Si.type     = 2*v;
                Si.r        = r2;
                S           = addstruct(S,Si);
                clear Si
            end;
        end;
        
        s = style.custom({'blue','black'});
        
        plt.subplot(221);
        plt.bar(S.type,S.r,'subset',S.group==1 & S.reg==1,'split',S.type,'leg','none','style',s);
        plt.subplot(222);
        plt.bar([S.type],S.r,'subset',S.group==1 & S.reg==2,'split',S.type,'leg','none','style',s);
        plt.set(gcf,'xticklabel',{'spatial','RDM','spatial','RDM'},...
            'ytick',0:0.2:1,'match','ylim');
        
        plt.labels([],'Pearsons r','S1','D',221);
        plt.labels([],'Pearsons r','M1',[],222);
        
        plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'1x2');
    case 'REL_betweenSubj_ttest'
        for i = 1:5 % motor condition (combining both left and right hemispheres)
            
            r=unique(regionNum);
            roi = [r(i), r(i+5)];
            g1 = roi_dystonia('REL_betweenSubj', 'region', roi, 'group', 1);
            g2 = roi_dystonia('REL_betweenSubj', 'region', roi, 'group', 2);
            fprintf(regionName{i})
            ttest(g1,g2,2,'independent')
        end
        for i = 1:5 % repeat for sensory condition (left hemisphere only)
            r=unique(regionNum);
            roi = [r(i)];
            g1 = roi_dystonia('REL_betweenSubj', 'region', roi, 'group', 1, 'stimType', 1);
            g2 = roi_dystonia('REL_betweenSubj', 'region', roi, 'group', 2, 'stimType', 1);
            fprintf(regionName{i})
            ttest(g1,g2,2,'independent')
        end
    case 'REL_intraHem' %
        dim = 5;
        roi      = [2 12];
        stimType    = 0;
        group       = 1;
        vararginoptions(varargin,{'region','stimType','group'})
        
        % 1. Get subset of data
        D = load(fullfile(regDir,'reg_distance_raw.mat'));
        D = getrow(D,ismember(D.region,roi) & D.stimtype==stimType & ...
            D.group==group & D.hand~=D.regSide);
        
        % 2. MANOVA fixed factor hemisphere (rmANOVA not enough data)
        [u,s,v] = svd(D.dist);
        D.distRed = u(:,1:dim)*s(1:dim,1:dim);
        % MANOVArp(D.regSide, D.SN, D.distRed);
        MANOVA1(D.regSide, D.distRed);
    case 'plot_bar' % plots psc and mean distance, fprint t-test
        % df2_fig('barplots', 1)
        D = load('reg_distance_raw.mat');
        regType=varargin{1};
        HAND=[0 1 1];
        STIMTYPE=[0 0 1];
        plotname={'Left Motor', 'Right Motor', 'Right Sensory'};
        facecolor={[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1]};
        edgecolor={[0 0 1],[1 0 0],[0 0 1],[1 0 0],[0 0 1],[1 0 0]};
        
        % 1. get region data (LH & RH)
        r = getrow(D,D.region==regType | D.region==(regType+10));
        
        % 2. plot ipsilateral barplots % print stats on screen
        ipsi = getrow(r, r.hand==r.regSide);
        figure(1), subplot(2,1,1), barplot([ipsi.hand,ipsi.stimtype], ipsi.psc, 'split', ipsi.group, 'facecolor', facecolor, 'edgecolor', edgecolor, 'gapwidth', [0.2 0.1 0.05]);
        set(gca, 'XTickLabel', []), ylabel([])% ylabel({(regionName{regType}); 'Ipsilateral'}), set(gca, 'XTickLabel', []) %{'LM' , 'LM', 'RM' , 'RM', 'RS' , 'RS'}, 'fontSize', 8), title('percentage signal change')
        for i = 1:3;
            tmp = getrow(ipsi, ipsi.hand==HAND(i) & ipsi.stimtype==STIMTYPE(i));
            fprintf('hand %d stimtype %d \n', HAND(i), STIMTYPE (i))
            [t,p] = ttest(tmp.psc(tmp.group==1), tmp.psc(tmp.group==2),2,'independent')
            clear tmp
        end
        figure(2), subplot(2,1,1), barplot([ipsi.hand,ipsi.stimtype], (mean(ipsi.dist,2)), 'split', ipsi.group,'facecolor', facecolor, 'edgecolor', edgecolor, 'gapwidth', [0.2 0.1 0.05]);
        set(gca, 'XTickLabel', []), ylabel([])% ylabel({(regionName{regType}); 'Ipsilateral'}), set(gca, 'XTickLabel', {'Left Motor' , 'Right Motor' , 'Right Sensory'}, 'fontSize', 8), title('mean distance')
        % mean distance t-tests
        for i = 1:3;
            tmp = getrow(ipsi, ipsi.hand==HAND(i) & ipsi.stimtype==STIMTYPE(i));
            fprintf('hand %d stimtype %d \n', HAND(i), STIMTYPE (i))
            [t,p] = ttest(mean(tmp.dist(tmp.group==1),2), mean(tmp.dist(tmp.group==2),2), 2,'independent')
            clear tmp
        end
        clear ipsi
        
        % 3. plot contralateral barplot % print stats on screen
        contra = getrow(r, r.hand~=r.regSide);
        hold on, figure(1), subplot(2,1,2),barplot([contra.hand,contra.stimtype], contra.psc, 'split', contra.group, 'facecolor', facecolor, 'edgecolor', edgecolor, 'gapwidth', [0.2 0.1 0.05]);
        set(gca, 'XTickLabel', []), ylabel([])%ylabel({(regionName{regType}); 'Contralateral'}), set(gca, 'XTickLabel', {'LM' , 'LM', 'RM' , 'RM', 'RS' , 'RS'}, 'fontSize', 8), title('percentage signal change')
        for i = 1:3; % ttest
            tmp = getrow(contra, contra.hand==HAND(i) & contra.stimtype==STIMTYPE(i));
            fprintf('hand %d stimtype %d \n', HAND(i), STIMTYPE (i))
            [t,p] = ttest(tmp.psc(tmp.group==1), tmp.psc(tmp.group==2),2,'independent')
            clear tmp
        end
        hold on, figure(2), subplot(2,1,2),barplot([contra.hand,contra.stimtype], mean(contra.dist,2), 'split', contra.group, 'facecolor', facecolor, 'edgecolor', edgecolor, 'gapwidth', [0.2 0.1 0.05]);
        set(gca, 'XTickLabel', []), ylabel([])%ylabel({(regionName{regType}); 'Contralateral'}), set(gca, 'XTickLabel', {'LM' , 'LM', 'RM' , 'RM', 'RS' , 'RS'}, 'fontSize', 8), title('mean distance')
        for i = 1:3;
            tmp = getrow(contra, contra.hand==HAND(i) & contra.stimtype==STIMTYPE(i));
            fprintf('hand %d stimtype %d \n', HAND(i), STIMTYPE (i))
            [t,p] = ttest(mean(tmp.dist(tmp.group==1),2), mean(tmp.dist(tmp.group==2),2),2,'independent')
            clear tmp
        end
        clear contra
    case 'plot_line' % lineplots, fprint MANOVA
        % df1_imana_roi('plot_line',1,'dist');
        
        color={[0 0 1],[1 0 0],[0 0 1],[1 0 0]};
        style_s={'-','-',':',':'};
        D=load(fullfile(regDir,'reg_distance_raw.mat'));
        regType=varargin{1};
        type = varargin{2};
        HAND=[0 1 1];
        STIMTYPE=[0 0 1];
        ipsi_x=[0.2, 0.3, 0.5, 0.1, 0.1];
        contra_x=[1.4, 0.8, 0.5,  0.1, 0.1];
        plotname={'Left Motor','Right Motor','Right Sensory'};
        xticklabels={'1/2', '1/3', '1/4', '1/5', '2/3', '2/4', '2/5', '3/4', '3/5', '4/5'};
        % D.ndist=bsxfun(@rdivide,D.dist,sqrt(sum(D.dist.^2,2)));
        
        %Ipsi
        for i=1:3
            subplot(2,3,i);
            traceplot(1:10,D.(type),'split',[D.group],...
                'subset',   D.regType==regType & D.hand==HAND(i) & D.stimtype==STIMTYPE(i) & D.regSide==D.hand, ...
                'linecolor',color,'patchcolor',color,'linestyle',style_s,'errorfcn','stderr');
            set(gca,'YLim',[-0.02 ipsi_x(find(unique(D.regType)==regType))], 'fontsize', 8, 'XTick', [1:10], 'XTickLabel', xticklabels);
            ylabel({(regname{regType}); 'Ipsilateral'});
            d = getrow(D,   D.regType==regType & D.hand==HAND(i) & D.stimtype==STIMTYPE(i) & D.regSide==D.hand);
            title(plotname{i})
            MANOVA1(d.group,d.dist)
        end;
        %Contra
        for i=1:3
            subplot(2,3,i+3);
            traceplot(1:10,D.(type),'split',[D.group],...
                'subset',   D.regType==regType & D.hand==HAND(i) & D.stimtype==STIMTYPE(i) & D.regSide~=D.hand, ...
                'linecolor',color,'patchcolor',color,'linestyle',style_s,'errorfcn','stderr' );
            d = getrow(D,   D.regType==regType & D.hand==HAND(i) & D.stimtype==STIMTYPE(i) & D.regSide~=D.hand);
            set(gca,'YLim',[-0.02 contra_x(find(unique(D.regType)==regType))], 'fontsize', 8, 'XTick', [1:10], 'XTickLabel', xticklabels);
            ylabel({(regname{regType}); 'Contralateral'});
            title(plotname{i})
            MANOVA1(d.group,d.dist)
        end;
    case 'laterality_index'
        D=load(fullfile(regDir,'reg_distance_raw.mat'));
        regType = varargin{1};
        % type = varargin{2};
        for i = 1:length(subj_name);
            % get subset of data
            d = getrow(D, D.SN==i & D.stimtype==0 & D.regType==regType);
            % build variables into data stuct 'l'
            l.mdist_ipsi(i)    = mean(mean(d.dist(d.regSide==d.hand,:))); % ipsi mean
            ipsi               = d.dist(d.regSide==d.hand,:);
            l.corr_ipsi(i)     = corr(ipsi(1,:)', ipsi(2,:)');
            l.mdist_contra(i)  = mean(mean(d.dist(d.regSide~=d.hand,:))); % contra mean
            l.group(i)         = d.group(1);
            l.SN(i)            = d.SN(1);
            contra             = d.dist(d.regSide~=d.hand,:);
            l.corr_contra(i)   = corr(contra(1,:)', contra(2,:)');
            l.lat_index(i)     = (l.mdist_contra - l.mdist_ipsi)/(l.mdist_contra+l.mdist_ipsi);
        end
        
        % save('lat_index', '-struct', 'l')
        
        
        % % compare motor task in S1, controls then dystonia
        % d = getrow(D, D.regSide~=D.hand & D.regType==regType & D.stimtype==0 & D.group==1);
        % MANOVA1(d.hand, d.dist)
        % d = getrow(D, D.regSide~=D.hand & D.regType==regType & D.stimtype==0 & D.group==2);
        % MANOVA1(d.hand, d.dist)
        % end
        
    case 'Fig_repStructurePassive'          % representational structure across groups
        D = load(fullfile(regDir,'reg_distance_raw.mat'));
        D.normdist = bsxfun(@rdivide,D.dist,sqrt(mean(D.dist.^2,2)));
        
        % D = rmfield(D,'subj');
        D1 = getrow(D,D.stimtype==1 & D.hand~=D.regSide & D.regType==1);
        D2 = getrow(D,D.stimtype==1 & D.hand~=D.regSide & D.regType==2);
        
        style.use('group_smallmarker');
        
        plt.subplot(1,3,2);
        plt.trace([],D1.normdist,'split',D1.group,'leg','none');
        plt.subplot(1,3,3);
        plt.trace([],D2.normdist,'split',D2.group,'leg','none');
        
        plt.set(gcf,'xtick',1:10,'ytick',0:0.5:1.5,'ylim',[-0.1 2],'yprecision','%1.1f');
        plt.labels([],{'dissimilarity','(normalized)'},'S1 (RH)','B',132);
        plt.labels([],[],'M1 (RH)','C',133);
        plt.set(gcf,'ax','square','yprecision','%1.1f');
        anot.hide_axis('y');
        plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'1x2');
    case 'Fig_repStructureActive'          % representational structure across groups
        D = load(fullfile(regDir,'reg_distance_raw.mat'));
        D.normdist = bsxfun(@rdivide,D.dist,sqrt(mean(D.dist.^2,2)));
        
        % D = rmfield(D,'subj');
        D1 = getrow(D,D.stimtype==0 & D.hand~=D.regSide & D.regType==1);
        D2 = getrow(D,D.stimtype==0 & D.hand~=D.regSide & D.regType==2);
        
        style.use('group_smallmarker');
        
        plt.subplot(141);
        plt.trace([],D1.normdist,'split',D1.group,'subset',D1.hand==0,'leg','off');
        
        plt.subplot(142);
        plt.trace([],D2.normdist,'split',D2.group,'subset',D2.hand==0,'leg','off');
        anot.hide_axis('y');
        plt.subplot(143);
        plt.trace([],D1.normdist,'split',D1.group,'subset',D1.hand==1,'leg','off');
        anot.hide_axis('y');
        plt.subplot(144);
        plt.trace([],D2.normdist,'split',D2.group,'subset',D2.hand==1,'leg','off');
        anot.hide_axis('y');
        
        plt.set(gcf,'xtick',1:10,'ytick',0:0.5:1.5,'ylim',[-0.1 2],'yprecision','%1.1f');
        plt.labels([],{'dissimilarity','normalized'},'S1 (LH)','A',141);
        plt.labels([],[],'M1 (LH)','B',142);
        plt.labels([],[],'S1 (RH)','C',143);
        plt.labels([],[],'M1 (RH)','D',144);
        
        plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'1x2');
    case 'FIG_somatotopicDistancesPassive'
        D = load(fullfile(regDir,'spatial_distances.mat'));
        D = getrow(D,ismember(D.region,[1 2]) & D.metric==3 & D.condition==2);
        
        style.use('group_smallmarker');
        
        plt.subplot(121);
        plt.trace(1:10,D.dist,'split',D.group,'subset',D.region==1,'leg','none');
        plt.labels([],'distance (cm)','CoG distances','A');
        plt.subplot(122);
        plt.trace(1:10,D.dist,'split',D.group,'subset',D.region==2,...
            'leg',{'non-dystonic','dystonic'});
        plt.set(gcf,'match','ylim','xtick',1:10);
        plt.labels([],[],'CoG distances','B');
        
        plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'0.75x2');
        
    case 'FIG_somatotopicDistancesActive'
        D = load(fullfile(regDir,'spatial_distances.mat'));
        D = getrow(D,ismember(D.region,[1 2]) & D.metric==3 & D.condition==1);
        
        style.use('group_smallmarker');
        
        plt.subplot(141);
        plt.trace(1:10,D.dist,'split',D.group,'subset',D.region==1 & D.hand==1,'leg','off');
        plt.labels([],'distance (cm)','CoG distances','A');
        
        plt.subplot(142);
        plt.trace(1:10,D.dist,'split',D.group,'subset',D.region==2 & D.hand==1,'leg','off');
        
        plt.subplot(143);
        plt.trace(1:10,D.dist,'split',D.group,'subset',D.region==1 & D.hand==2,'leg','off');
        
        plt.subplot(144);
        plt.trace(1:10,D.dist,'split',D.group,'subset',D.region==2 & D.hand==2,'leg','off');
        
        
        plt.set(gcf,'match','ylim','xtick',1:10);
        
        plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'1x2');
    case 'FIG_meanActivation'
        D = load(fullfile(regDir,'mean_betas.mat'));
        D1 = getrow(D,D.stimType==0 & D.regType==1 & D.hand~=D.regSide);
        D2 = getrow(D,D.stimType==0 & D.regType==2 & D.hand~=D.regSide);
        
        style.use('group_smallmarker');
        
        plt.subplot(3,2,[1 3]);
        plt.line([D1.hand D1.digit],D1.meanAct,'split',D1.group,...
            'leg',{'non-dystonic','dystonic'},'leglocation','south');
        plt.drawline(0,'style',style.custom('black','linewidth',1));
        plt.labels('finger',{'parameter','estimates (a.u.)'},'S1','A');
        plt.subplot(3,2,[2 4]);
        plt.line([D2.hand D2.digit],D2.meanAct,'split',D2.group,'leg','none');
        plt.drawline(0,'style',style.custom('black','linewidth',1));
        plt.set(gcf,'ylim',[-0.2 2.9]);
        plt.labels('finger',{'parameter','estimates (a.u.)'},'M1','B');
        
        plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'style','1x2');
    case 'Fig_MDS_Active'
        D = load(fullfile(regDir,'reg_distance_raw.mat'));
        % D = rmfield(D,'subj');
        Dc = getrow(D,D.stimtype==0 & D.hand~=D.regSide & D.regType==1 & D.group==1 & D.hand==1);
        D1 = getrow(D,D.stimtype==0 & D.hand~=D.regSide & D.regType==1 & D.group==2 & D.hand==0);
        D2 = getrow(D,D.stimtype==0 & D.hand~=D.regSide & D.regType==1 & D.group==2 & D.hand==1);
        
        % normalize distances
        Dc.dist = bsxfun(@rdivide,Dc.dist,sqrt(mean(Dc.dist.^2,2)));
        D1.dist = bsxfun(@rdivide,D1.dist,sqrt(mean(D1.dist.^2,2)));
        D2.dist = bsxfun(@rdivide,D2.dist,sqrt(mean(D2.dist.^2,2)));
        
        [yc,~] = cmdscale(mean(Dc.dist,1));
        [y1,~] = cmdscale(mean(D1.dist,1));
        [y2,~] = cmdscale(mean(D2.dist,1));
        
        yc(end+1,:) = yc(1,:);
        y1(end+1,:) = y1(1,:);
        y2(end+1,:) = y2(1,:);
        
        T = [];
        S.label = [1:5 1]';
        S.cat   = [1:6]';
        S.y     = yc;
        S.group = zeros(6,1);
        T       = addstruct(T,S);
        S.y     = y1;
        S.group = 1+zeros(6,1);
        T       = addstruct(T,S);
        S.y     = y2;
        S.group = 2+zeros(6,1);
        T       = addstruct(T,S);
        
        style.use('groupx3');
        plt.subplot(121); hold on;
        plt.scatter(T.y(:,1),-T.y(:,2),'label',T.label,'regression','none');
        plt.xy(T.y(:,1),-T.y(:,2),T.cat,'split',T.group,'markersize',5);
        
        plt.set(gcf,'ax','equal');
        plt.labels('dim 1','dim 2','MDS','E',121);
        % plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'style','1x2');
    case 'Fig_MDS_Passive'
        D = load(fullfile(regDir,'reg_distance_raw.mat'));
        
        Dc = getrow(D,D.stimtype==1 & D.hand~=D.regSide & D.regType==1 & D.group==1 & D.hand==1);
        D2 = getrow(D,D.stimtype==1 & D.hand~=D.regSide & D.regType==1 & D.group==2 & D.hand==1);
        
        % normalize distances
        Dc.dist = bsxfun(@rdivide,Dc.dist,sqrt(mean(Dc.dist.^2,2)));
        D2.dist = bsxfun(@rdivide,D2.dist,sqrt(mean(D2.dist.^2,2)));
        
        [yc,~] = cmdscale(mean(Dc.dist,1));
        [y2,~] = cmdscale(mean(D2.dist,1));
        
        yc(end+1,:) = yc(1,:);
        y2(end+1,:) = y2(1,:);
        
        T = [];
        S.label = [1:5 1]';
        S.cat   = [1:6]';
        S.y     = yc;
        S.group = zeros(6,1);
        T       = addstruct(T,S);
        S.y     = y2;
        S.group = 1+zeros(6,1);
        T       = addstruct(T,S);
        
        style.use('group');
        plt.subplot(121); hold on;
        plt.scatter(T.y(:,1),-T.y(:,2),'label',T.label,'regression','none','facealpha',1);
        plt.xy(T.y(:,1),-T.y(:,2),T.cat,'split',T.group,'markersize',5);
        
        plt.set(gcf,'ax','equal');
        plt.labels('dim 1','dim 2','MDS','E',121);
        % plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'1x2');
    case 'Fig_repStructureCerebellumAnterior'          % representational structure across groups
        D = load(fullfile(regDir,'reg_distance_raw.mat'));
        D = rmfield(D,'subj');
        D1 = getrow(D,D.stimtype==0 & D.hand==D.regSide & D.regType==9);
        D2 = getrow(D,D.stimtype==0 & D.hand~=D.regSide & D.regType==9);
        
        style.use('group_smallmarker');
        
        plt.subplot(221);
        plt.trace([],D1.dist,'split',D1.group,'subset',D1.hand==0,'leg','none');
        plt.subplot(222);
        plt.trace([],D2.dist,'split',D2.group,'subset',D2.hand==0,'leg','none');
        
        plt.subplot(223);
        plt.trace([],D1.dist,'split',D1.group,'subset',D1.hand==1,'leg','none');
        plt.subplot(224);
        plt.trace([],D2.dist,'split',D2.group,'subset',D2.hand==1,'leg','none');
        
        plt.set(gcf,'match','ylim','xtick',1:10);
        plt.labels([],'dissimilarity (a.u.)','Ipsi Anterior','A',221);
        plt.labels([],[],'Contra Anterior','B',222);
        plt.labels([],'dissimilarity (a.u.)',[],[],223);
        
        plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'style','2x2');
    case 'Fig_repStructureCerebellumPosterior'          % representational structure across groups
        D = load(fullfile(regDir,'reg_distance_raw.mat'));
        D = rmfield(D,'subj');
        D1 = getrow(D,D.stimtype==0 & D.hand==D.regSide & D.regType==10);
        D2 = getrow(D,D.stimtype==0 & D.hand~=D.regSide & D.regType==10);
        
        style.use('group_smallmarker');
        
        plt.subplot(221);
        plt.trace([],D1.dist,'split',D1.group,'subset',D1.hand==0,'leg','none');
        plt.subplot(222);
        plt.trace([],D2.dist,'split',D2.group,'subset',D2.hand==0,'leg','none');
        
        plt.subplot(223);
        plt.trace([],D1.dist,'split',D1.group,'subset',D1.hand==1,'leg','none');
        plt.subplot(224);
        plt.trace([],D2.dist,'split',D2.group,'subset',D2.hand==1,'leg','none');
        
        plt.set(gcf,'match','ylim','xtick',1:10);
        plt.labels([],'dissimilarity (a.u.)','Ipsi Posterior','A',221);
        plt.labels([],[],'Contra Posterior','B',222);
        plt.labels([],'dissimilarity (a.u.)',[],[],223);
        
        plt.save(fullfile(figureDir,sprintf('%s.pdf',what)),'style','2x2');
        
    case 'Fig_cog_reliability'
        D = load(fullfile(regDir,'spatial_distances_splithalf.mat'));
        
        % Looping over subjects/regions/conditions/hand
        S = [];
        for ty=1:2
            Dt = getrow(D,D.type==ty);
            
            for sn = unique(D.sn)'
                Ds = getrow(Dt,Dt.sn==sn);
                
                for reg = unique(Ds.region)'
                    Dr = getrow(Ds,Ds.region==reg);
                    
                    for cond = unique(Dr.condition)'
                        Dc = getrow(Dr,Dr.condition==cond);
                        
                        for h = unique(Dc.hand)'
                            Dh  = getrow(Dc,Dc.hand==h);
                            
                            for m = unique(Dh.metric)'
                                Dm     	= getrow(Dh,Dh.metric==m);
                                [x,dig] = pivottablerow(Dm.digit,[Dm.x Dm.y],'nanmean(x,1)');
                                d       = pdist(x,'euclidean');
                                
                                % Saving results
                                Si.sn           = sn;
                                Si.metric       = m;
                                Si.region       = reg;
                                Si.condition    = cond;
                                Si.hand         = h;
                                Si.group        = mean(Dm.group);
                                Si.type         = ty;
                                Si.dist         = d;
                                Si.x            = x(:,1)';
                                Si.y            = x(:,2)';
                                Si.digit        = dig';
                                S               = addstruct(S,Si);
                                
                                fprintf('Type: %d, SN: %d, Reg: %d, Cond: %d, Hand: %d, Metric: %d\n',ty,sn,reg,cond,h,m);
                            end;
                        end;
                    end;
                end;
            end;
        end;
        
        varargout = {S};
        save(fullfile(regDir,'spatial_distances_splithalf.mat'),'-struct','S');
        
    case 'Fig_spatial_distances_groups'
        D = load(fullfile(regDir,'spatial_distances.mat'));
        D = getrow(D,D.region==1 & D.condition==2);
        
        style.reset;
        plt.subplot(121);
        plt.trace([],D.dist,'split',D.k,'subset',D.group==1,'leglocation','northeast','leg','off');
        plt.set('xtick',1:10);
        plt.subplot(122);
        plt.trace([],D.dist,'split',D.k,'subset',D.group==2,'leglocation','northeast');
        plt.set('xtick',1:10);
        
        plt.labels([],{'euclidean distance','between CoGs'},'controls',[],121);
        plt.labels([],[],'dystonic',[],122);
        plt.legend('northeast',{'weighted mean','softmax, k=0.2','softmax, k=0.4','softmax, k=0.6','softmax, k=0.8','peak voxel'});
        
        anot.hide_axis('y');
        plt.match('y');
        
    case 'Fig_spatial_distances_reliability'
        % S1
        D = load(fullfile(regDir,'spatial_distances_splithalf.mat'));
        D = getrow(D,D.region==1 & D.condition==2);
        
        style.reset;
        D1a = getrow(D,D.type==1 & D.group==1);
        D1b = getrow(D,D.type==2 & D.group==1);
        D2a = getrow(D,D.type==1 & D.group==2);
        D2b = getrow(D,D.type==2 & D.group==2);
        
        % reliability measured through correlation
        D1a.r = diag(corr(D1a.dist',D1b.dist'));
        D2a.r = diag(corr(D2a.dist',D2b.dist'));
        T = addstruct(D1a,D2a);
        
        plt.subplot(121);
        plt.line(T.k,T.r,'plotfcn','nanmean','split',T.group,'leg',{'controls','dystonic'},'leglocation','southeast');
        plt.labels([],'Pearsons r','Split-half Rel (corr)','S1');
        plt.set('xtick',0:0.2:1);
        plt.set('xticklabel',{'w-mean','k=0.2','k=0.4','k=0.6','k=0.8','peak'});
        
        % M1
        D = load(fullfile(regDir,'spatial_distances_splithalf.mat'));
        D = getrow(D,D.region==2 & D.condition==2);
        
        style.reset;
        D1a = getrow(D,D.type==1 & D.group==1);
        D1b = getrow(D,D.type==2 & D.group==1);
        D2a = getrow(D,D.type==1 & D.group==2);
        D2b = getrow(D,D.type==2 & D.group==2);
        
        % reliability measured through correlation
        D1a.r = diag(corr(D1a.dist',D1b.dist'));
        D2a.r = diag(corr(D2a.dist',D2b.dist'));
        T = addstruct(D1a,D2a);
        
        plt.subplot(122);
        plt.line(T.k,T.r,'plotfcn','nanmean','split',T.group,'leg',{'controls','dystonic'},'leglocation','southeast');
        plt.labels([],'Pearsons r','Split-half Rel (corr)','M1');
        plt.set('xtick',0:0.2:1);
        plt.set('xticklabel',{'w-mean','k=0.2','k=0.4','k=0.6','k=0.8','peak'});
        
        plt.set('ylim',[-0.2 1.0]);
    case 'Fig_reliability_rsa'
        D = load(fullfile(regDir,'reg_distance_raw_splithalf.mat'));
        
        % stimtime = 1
        D = getrow(D,D.regSide ~= D.hand);
        D = getrow(D,D.stimtype == 1);
        D.r = diag(corr(D.dist1',D.dist2'));
        
        plt.box(D.region,D.r,'split',D.group);
        ylim([-0.2, 1]);
        plt.labels([],'Pearsons r','RSA Reliability');
        plt.legend('southeast',{'controls','dystonic'});
        plt.set('xticklabel',{'S1','S1','M1','M1'});
        
    case 'Calc_reliability' % Calculate split-half reliability
        D = load(fullfile(regDir,'spatial_distances_splithalf.mat'));
        R = load(fullfile(regDir,'reg_distance_raw_splithalf.mat'));
        T=[];
        for s=1:17
            for r=1:2 
                for h=1:2
                    for c=1:2
                        for m=1:max(D.metric)
                            i1=find(D.sn==s & D.region==r & D.hand==h & D.condition==c & D.metric==m & D.type==1);
                            i2=find(D.sn==s & D.region==r & D.hand==h & D.condition==c & D.metric==m & D.type==2);
                            if (~isempty(i1))
                                if length(i1)>1
                                    keyboard; % Must be an error
                                end
                                TT.sn=s;
                                TT.region = r;
                                TT.hand = h;
                                TT.condition = c;
                                TT.metric=m;
                                TT.k = D.k(i1);
                                TT.corr=corr(D.dist(i1,:)',D.dist(i2,:)');
                                T=addstruct(T,TT);
                            end
                        end
                        i=find(R.SN==s & R.regType==r & R.hand==h-1 & R.stimtype==c-1 & R.regSide~=R.hand);
                        if (~isempty(i1))
                            if length(i1)>1
                                keyboard; % Must be an error
                            end
                            TT.sn=s;
                            TT.region = r;
                            TT.hand = h;
                            TT.condition = c;
                            TT.metric=11;
                            TT.k = 0;
                            TT.corr=corr(R.dist1(i,:)',R.dist2(i,:)');
                            T=addstruct(T,TT);
                        end
                    end
                end
            end
        end
        save(fullfile(regDir,'reliability.mat'),'-struct','T');
    case 'Fig_reliability'  % Plot split-half relaibilities
        c=2;
        T=load(fullfile(regDir,'reliability.mat'));
        T.metricType=T.metric;
        T.metricType(T.metric>1 & T.metric<10)=3;
        T.metricType(T.metric==10)=2;
        T.metricType(T.metric==11)=4;
        for r=1:2
            subplot(1,2,r);
            lineplot([T.metricType T.k],T.corr,'subset',T.region==r & ~isnan(T.corr) & T.condition==c,'style_thickline');
            set(gca,'YLim',[0 1]);
        end;
         set(gcf,'PaperPosition',[2 2 6 2.3])
        wysiwyg

    case 'make_distance'
        df1_imana('SPAT_metrics2d','isplot',0);
        df1_imana('SPAT_metrics2d_splithalf','isplot',0);
        df1_imana('SPAT_estimateDistances');
        df1_imana('SPAT_estimateDistances_splithalf');
        df1_imana('Calc_reliability');
        df1_imana('Fig_reliability');
    case 'power_analysis'
        D=load('reg_distance_raw.mat'); 
        T=getrow(D,D.stimtype==1 & D.regSide==0 & D.regType==1); 
        mDist = sqrt(mean(T.dist,2));
        aDist = mean(mDist)*(1-0.29); 
        
        ttest(mDist(T.group==1),mDist(T.group==2),2,'independent'); 
        SE=sqrt(2*std(mDist).^2/8); 
        expectedT = (mean(mDist)-aDist)/SE;
        power = tcdf(expectedT - tinv(0.95,15),15); 
        fprintf('Power for 29 percent reduction: %f  ',power);
    case 'Bayesian_analysis' 
        D=load('reg_distance_raw.mat'); 
        T=getrow(D,D.stimtype==1 & D.regSide==0 & D.regType==1); 
        mDist = sqrt(mean(T.dist,2));
        SE=sqrt(2*std(mDist).^2/8); 
        aDist = mean(mDist(T.group==1))*(1-0.29);
        expEff = (mean(mDist(T.group==1)-aDist))/SE; 
        trueEff = (mean(mDist(T.group==1))-mean(mDist(T.group==2)))/SE; 
        BayesFact = tpdf(-trueEff,15)/tpdf(expEff-trueEff,15); 
        fprintf('bayes factor against a 29 percent reduction: %f  ',BayesFact);        
        
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
