function varargout=df1_imana(what,varargin)
% baseDir=        fullfile('/Users','tob', 'Projects','FingerPattern_dystonia');
%  baseDir=        fullfile('/media','DATA', 'Projects','FingerPattern_dystonia');
% baseDir=        fullfile('/Volumes/MacintoshHD2/fingerPattern_dystonia');
% baseDir=        '/Users/naveed/Documents/data/FingerPattern_dystonia';
% baseDir=        fullfile('~/Projects/fingerPattern_dystonia');
% % baseDir=        fullfile('/Volumes/MotorControl/project/FingerPattern_dystonia);
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

% All these variables should go into the patient_list.txt document 
subj_MTname={   'MT02554', 'MT02613', 'MT02614', 'MT02689', 'MT02815', ...
    'MT02816', 'MT02852', 'MT02851', 'MT02855', 'MT02853', ...
    'MT02850', 'MT02859', 'MT02860'};%'MT02481','MT02481',
subj_name={ 'd01', 's01', 'd02', 's02', 'd04', ...
            's04', 'd06', 's03', 'd07', 'd08', ...
            'd09', 'd10', 'd11'};  
subj_group=[2 1 2 1 2 ...
    1 2 1 2 2 ...
    2 2 2]';  % 2 for dystonic 1 for pianists
subj_NumVol= [  150 144 144 144 144 ...
    144 144 144 144 144 ...
    144 144 144]; 

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
    case 'make_nii'%________________STEP 1____________________________
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
                spmj_tar2nii ([subj_MTname{sn} '.r', num2str(i) '.tar'], [subj_name{sn} '_run' run{i} '.nii'], 'startTR', dummyScans+1)
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
    case 'setAC'%___________________STEP 1.0__________________________
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
            -95 -140 -164];     %d11
        %----set the AC
        V= spm_vol(ana_name);
        R= spm_read_vols(V);
        V.mat(1:3,4)= AC_coord(sn,:);
        V.mat
        spm_write_vol(V,R);
    case 'segmentation'%____________STEP 1.1__________________________
        % df1_imana('segmentation',1)
        sn=varargin{1};
        for s=sn
            spmj_segmentation(fullfile(anatomicalDir, subj_name{s}, [subj_name{s}, '_anatomical.nii']))
        end;
    case 'slice_timing'%____________STEP 2____________________________
        %df1_imana('slice_timing',1)
        % [2:2:32 1:2:32] interleave
        prefix='';
        sn=varargin{1};
        for r= 1:8
            for i=1:(subj_NumVol(sn))
                N{i} = [fullfile(baseDir, 'imaging_data_raw',subj_name{sn}, [prefix subj_name{sn},'_run',run{r},'.nii,',num2str(i)])];
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
    case 'makefieldmap'%____________STEP 3____________________________
        %df1_imana('makefieldmap', 1)
        prefix='a';
        %prefix='';
        sn=varargin{1};
        spmj_makefieldmap(baseDir, subj_name{sn}, run,'prefix',prefix, 'subfolderRawdata', '');
    case 'make_realign_unwarp'%_____STEP 4____________________________
        % df1_imana('make_realign_unwarp',2)
        %prefix='ab'; %with bias corred epi's
        prefix='a'; %with slice time corrected images
        %prefix=''; % on raw images
        sn=varargin{1}
        spmj_realign_unwarp(baseDir, subj_name{sn}, run, 1, subj_NumVol(sn),'prefix',prefix, 'subfolderRawdata','') %76
        % maybe produce mean.nii
        % spm_imcalc_ui({},'mean.nii','mean(X)',{1,[],[],[]})
    case 'move_images'%_____________STEP 5____________________________
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
    case 'meanimage_bias_correction'%STEP 5.1_________________________
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
    case 'plot_movementparameters'%___________________________________
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
    case 'coreg'%___________________STEP 6____________________________
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
    case 'coreg_2'%_________________STEP 6.1__________________________
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
    case 'coreg_3'%_________________STEP 6.2__________________________
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
    case 'make_samealign'%__________STEP 7____________________________
        % df1_imana('make_samealign',1)
        %prefix='ra';
        prefix='ua';
        %prefix='u';
        sn=varargin{1};
        Q={};
        cd(fullfile(baseDir, 'imaging_data',subj_name{sn}));
        for r= 1:numel(run)
            for i=1:(subj_NumVol(sn))
                Q{end+1} = [fullfile(baseDir, 'imaging_data',subj_name{sn}, [prefix, subj_name{sn},'_run',run{r},'.nii,',num2str(i)])];
            end
        end
        P{1} = fullfile(baseDir, 'imaging_data', subj_name{sn}, ['meanepi_' subj_name{sn} '.nii']);
        spmj_makesamealign_nifti(char(P),char(Q))
        %delete(char(P))
    case 'check_samealign'%_________STEP 8____________________________
        % df1_imana('check_samealign',1)
        prefix='ua';
        %prefix='u';
        sn=varargin{1};
        Q={};
        cd(fullfile(baseDir, 'imaging_data',subj_name{sn}));
        for r= 1:numel(run)
            for i=1:(subj_NumVol(sn))
                Q{end+1} = [fullfile(baseDir, 'imaging_data',subj_name{sn}, [prefix, subj_name{sn},'_run',run{r},'.nii,',num2str(i)])];
            end
        end
        P{1}= fullfile(baseDir, 'imaging_data',subj_name{sn}, ['meanepi_' subj_name{sn} '.nii']);
        spmj_checksamealign(char(P),char(Q))
    case 'make_maskImage'%__________STEP 8.1__________________________
        % df1_imana('make_maskImage',1)
        % Do this on the parcelation maps
        sn=varargin{1};
        for s=sn
            cd(fullfile(baseDir,'imaging_data', subj_name{s}));
            %----load mean image
            nam{1}=  fullfile(baseDir, 'imaging_data',subj_name{s}, ['meanepi_' subj_name{s} '.nii']);
            %---- mask with c1,c2, & c3 from segmentation
            nam{2}= fullfile(anatomicalDir, subj_name{s},  ['c1' subj_name{s}, '_anatomical.nii'])
            nam{3}= fullfile(anatomicalDir, subj_name{s},  ['c2' subj_name{s}, '_anatomical.nii'])
            nam{4}= fullfile(anatomicalDir, subj_name{s},  ['c3' subj_name{s}, '_anatomical.nii'])
            spm_imcalc_ui(nam, 'rmask_noskull.nii', 'i1>1 & (i2+i3+i4)>0.2')
            %----reslice volume
            % V= spm_vol(nam{1});
            % spmj_reslice_vol('mask_noskull.nii', V.dim, V.mat, 'rmask_noskull.nii')
            
        end
    case 'make_glm_1'%________________STEP 9____________________________
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
    case 'make_glm_2'%__________One regressor per run, but with extra regressors for error trials
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
    case 'estimate_glm'%____________STEP 10___________________________
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
    case 'hrf_getdata'%_______________________________________________
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
    case 'hrf_plot'%__________________________________________________
        %df1_imana('hrf_plot',1)
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
    case 'contrast'%________________STEP 11
        %df1_imana('contrast',1, 1)
        sn=varargin{1};
        glmType=varargin{2};
        for s=sn
            cd(fullfile(baseDir,glmName{glmType}, subj_name{s}));
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
                            conName= sprintf('T_hem_%i_stimType_%i_hand_%i', h,s,d);
                            SPM.xCon(end+1)=spm_FcUtil('Set',conName, 'T', 'c',con',SPM.xX.xKXs);
                        end
                    end
                end
            end
            %
            %             % OPTIONAL: t contrast error against rest left
            %             for h=1:2
            %                 if (any(T.regType==2 & T.hand==h))
            %                     con=zeros(1,size(SPM.xX.X,2));
            %                     con(T.hand==h & T.regType==2)=1;
            %                     con=con/sum(con);
            %                     SPM.xCon(end+1)=spm_FcUtil('Set','t_err_left', 'T', 'c',con',SPM.xX.xKXs);
            %                 end;
            %             end;
            
            
            %____do the constrasts
            SPM=spm_contrasts(SPM,[1:length(SPM.xCon)]);
            save SPM SPM;
        end
    case 'SUIT_isolate'%____________DF_______________________________________
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
    case 'SUIT_normalize'%__________DF_______________________________________
        %df1_imana('SUIT_normalize',1)
        sn=varargin{1};
        %----got to suit DIR
        cd(fullfile(anatomicalDir,subj_name{sn},'suit'))
        %----run suit normalize
        suit_normalize(['c_',subj_name{sn},'_anatomical.nii'], 'mask',['c_',subj_name{sn},'_anatomical_pcereb_corr.nii'])
    case 'make_mask_suit'%__________DF_______________________________________
        %df1_imana('make_mask_suit',1)
        sn=varargin{1};
        mask_all= fullfile(baseDir, glmName{1},subj_name{sn}, 'mask.img');
        mask_cerebllum= fullfile(anatomicalDir,subj_name{sn},'suit',['c_',subj_name{sn},'_anatomical_pcereb_corr.nii'])
        mask_suit= fullfile(baseDir, glmName{1},subj_name{sn}, 'mask_suit.nii');
        %----calc suit mask
        spm_imcalc_ui({mask_all,mask_cerebllum}, mask_suit, 'i1>0 & i2>0');
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
            [LI,voxmin,voxmax,voxel,node,surfindx,depth,radvox]= lmva_voxelselection_surf(surf, [6 100],volDef);
            LI=LI';
            save(fullfile(refDir,subj_name{s}, 'vol_roi_100vox.mat'), 'LI','voxmin','voxmax','voxel','radvox');
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
            searchlightname=fullfile(refDir,subj_name{s},'vol_roi_100vox.mat');
            S=load(searchlightname);

            % Kill the voxels that belong to cortical searchlights
            volDef.mask(S.voxel)=0; 
            
            % Compute the search lights for the cerebellum and add them to
            % the file 
            [LI,voxmin,voxmax,voxel,radvox]= lmva_voxelselection_volume([20 100],volDef);
            S.LI=[S.LI;LI];
            S.voxel=[S.voxel;voxel];
            S.voxmin=[S.voxmin;voxmin];
            S.voxmax=[S.voxmax;voxmax];
            S.radvox=[S.radvox;radvox];
            save(fullfile(refDir,subj_name{s}, 'vol_roi_100vox.mat'), '-struct','S');
        end;
    case 'MVA_do'%__________________DF____________________________________
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
            
            %----maybe do it in the futur with this code
            %             idx=find(T.regType==1);
            %                 for i=idx'
            %                     beta_images{i}=  fullfile(glmDirSubj,SPM.Vbeta(i).fname);
            %                 end
            
            lda_names= {fullfile(workDir, subj_name{s}, 'acc_L_motor.nii'),...
                fullfile(workDir, subj_name{s}, 'acc_R_motor.nii'),...
                fullfile(workDir, subj_name{s}, 'acc_R_sens.nii')};
            %----do LDA
            lmva_spm('vol_roi_100vox.mat',beta_images,lda_names,ldafunction,'params',{T.run',T.hand',T.stimType',T.digit'})%
        end
    case 'MVA_do_dist'
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
            
            surfaces=load(fullfile(glmDirSubj,['vol_roi_100vox','.mat']));
            lmva_spm(surfaces,SPM.xY.P,outfiles,ldafunction,'params',{SPM,D});
        end;

    case 'MVA_zValue'%______________DF____________________________________
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
    case 'MNI_normalization_write'%___________________________________
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
                [dir,name,ext]=spm_fileparts(images{j});
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
    case 'MNI_mask_anatomical'%__________________________________________________
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
    case 'MNI_smooth'%________________________________________________
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
        % df1_imana('surf_freesurfer', sn);
        df1_imana('surf_xhemireg',sn);
        df1_imana('surf_map_ico',sn,2);
    case 'surf_freesurfer'%_________________1_________________________
        %df1_imana('surf_freesurfer', 5)
        sn=varargin{1};
        freesurfer_reconall(freesurferDir,subj_name{sn},fullfile(anatomicalDir,subj_name{sn},[subj_name{sn} '_anatomical.nii']));
    case 'surf_xhemireg'
        % df1_imana('surf_xhemireg',1)
        sn=varargin{1};
        for i=sn
            freesurfer_registerXhem({subj_name{i}},freesurferDir,'hemisphere',[1:2]);
        end;
    case 'surf_map_ico'%____________DF________________________________________
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
    case 'surf_make_caret'%_________DF________________________________________
        % df1_imana('surf_make_caret',3:5,2)
        sn=varargin{1};
        atlas=varargin{2};
        for i=sn
            caret_importfreesurfer([atlasA{atlas} subj_name{i}],freesurferDir,caretDir);
        end;
        %=========================================================================
        % Map functional volumes on surface over caret
        %=========================================================================
    
    case 'surf_map_con'%_________DF______________New mappint algorithm_____
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
    case 'surf_map_finger'%______DF______________New mappint algorithm_____
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
                
                C1=caret_load(white);
                C2=caret_load(pial);
                
                for f=1:length(fileList)
                    images{f}=fullfile(baseDir, glmName{glm},subj_name{s},fileList{f});
                end;
                M=caret_vol2surf_own(C1.data,C2.data,images,'ignore_zeros',1);
                caret_save(fullfile(caretSDir,[subj_name{s} '_finger.metric']),M);
            end;
        end;
    case 'surf_avrgsurfshape'%___DF______________Generates group metric____
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
    case 'surf_fingerpatterns'%__DF______________Plot finger patterns______
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
    case 'surf_spatial_cog'%_____ToDo______________________________________
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
    case 'surf_make_ROIpaint'%___DF______________Generates surface ROI paint file
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
    case 'surf_makeGroup'%_______
        atlas=2;
        
        INname={'func','func','func',...
            'accuracy_L_motor','accuracy_R_motor','accuracy_lda_R_sens'};
        OUTname={'func_L_motor','func_R_motor',...
            'func_R_sens',...
            'acc_L_motor','acc_R_motor','acc_R_sens'};
        
        inputcol= [1 2 3 1 1 1];
        replaceNaN=[1 1 1 0 0 0];
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
    case 'surf_zacc'% Make the zacc files (on group level) and remove the NaN
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
        Smooth_iterations='';
        vararginoptions(varargin,{'excludeSubj', 'Smooth_iterations'});
        if isempty(Smooth_iterations)
            Smooth_iterations=5; smoothdir= '';
        else
            smoothdir= [num2str(Smooth_iterations), 'iterSmooth'];
        end
        
        subj=1:11;
        
        SPMname={'func_left','func_right',...
            'zacc_left', 'zacc_right'};
        for h=1:2
            surfaceGroupDir=[caretDir filesep 'fsaverage'  filesep hemName{h} filesep smoothdir];
            cd(surfaceGroupDir)
            %----get the full directory name of the metric files and the smoothed metric files that we create below
            for i=1:length(SPMname);
                filenames{i}=[surfaceGroupDir filesep hem{h} '.' SPMname{i} '.metric']; % unsmoothed
            end;
            %----define name of coord and topology
            coordfile=[caretDir filesep 'fsaverage'  filesep hemName{h} filesep hem{h} '.FIDUCIAL.coord'];
            topofile=[caretDir filesep 'fsaverage'  filesep hemName{h} filesep hem{h} '.CLOSED.topo'];
            %----smooth the metric files and save them with the prefix 's'
            sfilenames=caret_smooth(filenames,'coord',coordfile,'topo',topofile,'iterations', Smooth_iterations);
        end;
    case 'surf_groupCon'              % Different group-GLM models  ....
        %sl1_pattern_behavior('surf_mapCorr',1:3)
        
        type=varargin{1};
        SPMname={'func_left','func_right',...
            'zacc_left', 'zacc_right'};
        vararginoptions({varargin{2:end}},{'SPMname'});
        numGroups=3;
        groupName={'ltSham','ltAnod','ltCath'};
        numSubj=16;
        numSPM=length(SPMname);
        
        for h=1:2
            surfaceGroupDir=[caretDir filesep 'fsaverage'  filesep hemName{h} ];
            cd(surfaceGroupDir);
            data=[];column_name={};
            outDir=[surfaceGroupDir filesep ];
            switch (type)
                case 'all'      % Left trained all
                    X=ones(numSubj,1);  % Averaged over all
                case 'sep'      % Left trained seperate
                    X=[];
                    for i=1:numGroups
                        X=[X double((tDCS_group==i)')];
                    end;
                case 'pre'
                    D=dload(fullfile(baseDir,'SequenceLearning','analyze','MT','df1_preposttest.dat'));
                    [X,R]=pivottable([D.SN D.tDCS],[],D.MT,'mean', 'subset',D.day==1); % despite hand differences, the left and right performance are correlated with 0.98
                    X=log(X);
                    X=X-mean(X);
                    % Add group difference and intercept
                    X=[ones(size(X,1),1) R(:,2)-mean(R(:,2)) X];
                otherwise
                    error('unknown type');
            end;
            
            for s=1:length(SPMname)
                %----get the full directory name of the metric files and the smoothed metric files that we create below
                filename=['s' hem{h} '.' SPMname{s} '.metric']; % smoothed
                outname=[type '_' SPMname{s} ];
                cSPM=caret_getcSPM('regression',X,'data',filename,[1:length(subj_name)],'no_intercept');
                save([outDir, 'cSPM.' outname, '.mat'],'cSPM');
                caret_savecSPM([outDir, 'stats.',outname,'.metric'],cSPM);
                switch(type)
                    case 'all'
                        data(:,s)=cSPM.con(1).con; % mean
                        data(:,s+length(SPMname))=cSPM.con(1).Z; % T
                        column_name{s}=['mean_' SPMname{s}];
                        column_name{s+length(SPMname)}=['T_' SPMname{s}];
                    case 'sep'
                        for g=1:numGroups
                            data(:,(s-1)*numGroups+g)=cSPM.con(g).con; % mean of group
                            data(:,(s-1)*numGroups+g+numSPM*numGroups)=cSPM.con(g).Z; % T
                            column_name{(s-1)*numGroups+g}=['mean_' SPMname{s} '_' groupName{g}];     %----intercept
                            column_name{(s-1)*numGroups+g+numSPM*numGroups}=['T_' SPMname{s} '_' groupName{g}];   %----slope
                        end;
                    case 'pre'
                        data(:,s)=cSPM.con(1).con; % mean
                        data(:,s+numSPM)=cSPM.con(2).con; % group diff
                        data(:,s+numSPM*2)=cSPM.con(3).con*1000; % Slope with MT
                        data(:,s+numSPM*3)=cSPM.con(1).Z; % T
                        data(:,s+numSPM*4)=cSPM.con(2).Z;
                        data(:,s+numSPM*5)=cSPM.con(3).Z;
                        
                        column_name{s}=['mean_' SPMname{s}];     %----intercept
                        column_name{s+numSPM}=['diff_' SPMname{s}];   %----slope
                        column_name{s+numSPM*2}=['slope_' SPMname{s}];   %----slope
                        column_name{s+numSPM*3}=['Tmean_' SPMname{s}];
                        column_name{s+numSPM*4}=['Tdiff_' SPMname{s}];
                        column_name{s+numSPM*5}=['Tslope_' SPMname{s}];
                end;
            end;
            C=caret_struct('metric','data',data,'column_name',column_name);
            caret_save([surfaceGroupDir '/' 'summary.' type '.metric'],C);
        end;
    case 'ROI_suit'%_____________DF______________Generates suit ROI nii file
        %df1_imana('suit_ROI', 1)
        sn=varargin{1};
        
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
        refImage= fullfile(baseDir, 'imaging_data',subj_name{sn}, ['meanepi_',subj_name{sn},'.nii']);
        defMat= fullfile(anatomicalDir,subj_name{sn},'suit',['mc_',subj_name{sn},'_anatomical_snc.mat']);
        source= fullfile(anatomicalDir,subj_name{sn},'suit', ['ROI_cerebellum.nii']);
        suit_reslice_inv(source, defMat,'reference',refImage,'prefix', '');
        movefile(fullfile(anatomicalDir,subj_name{sn},'suit', ['ROI_cerebellum.nii']), fullfile(anatomicalDir,subj_name{sn}, ['ROI_cerebellum.nii']))
    case 'ROI_define'%___________DF______________Reads the ROI definition from sequence learning paper
        %df1_imana('ROI_define', 1:4)
        sn= varargin{1};
        
        for s=sn
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
                R=region_calcregions(R);
                cd(regDir);
%                 for r=1:length(R);
%                     if (s==1 && ~isempty(R{r}))
%                         region_saveasimg(R{r},spm_vol(file));
%                     end;
%                 end;
                save([subj_name{s} '_regions.mat'],'R');
            end;
        end;
    case 'ROI_data'%_____________DF______________This gives an ROI structure: relies on the fact that the regions
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
            
            load(fullfile(regDir,[subj_name{s} '_regions.mat']));
            for i=1:length(j)
                P{i}=sprintf('beta_%4.4d.img',j(i));
            end;
            
            % Add a few extra images
            %----task against rest
            P{end+1}='con_0005.img'; %L_motor
            P{end+1}='con_0006.img'; %R_motor
            P{end+1}='con_0007.img'; %R_sens
            P{end+1}='ResMs.img';
            %----F-value figners
            P{end+1}='spmF_0002.img'; %L_motor
            P{end+1}='spmF_0003.img'; %R_motor
            P{end+1}='spmF_0004.img'; %R_sens
            
            E=getrow(SI,j); % Retain only regressors of interest
            
            V=spm_vol(char(P));
            data = region_getdata(V,R);
            
            for i=1:length(R)
                if (~isempty(R{i}))
                    D=[];
                    vec=ones(size(data{i},2),1);
                    D.beta=data{i}(1:length(j),:)';
                    D.SN=vec*s;
                    D.regNum=vec*i;
                    D.xyz=R{i}.data;
                    for z=1:3 %----task against rest
                        D.handAct(:,z)=data{i}(length(j)+z,:)';
                    end;
                    D.ResMs=data{i}(length(j)+4,:)';
                    for z=1:3%----digit F-test
                        D.Ftest(:,z)=data{i}(length(j)+z+4,:)';
                    end;
                    T=addstruct(T,D);
                end;
            end;
        end;
        regDir=fullfile(baseDir,'RegionOfInterest');
        cd(regDir);
        save(fullfile(regDir,['reg_data_' dataext '.mat']),'-struct','T');
        varargout={T};
    case 'ROI_stats'%____________DF______________caluculates some statistic on the region...
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
    case 'ROI_decomposition'%____DF______________
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
    case 'ROI_raw'                  % Do extraction of time series to LDA-t values
        % T=df1_imana('ROI_raw',[1:11]);
        % save('reg_distance_ldat.mat','-struct','T');
        
        selection='none';
        fcn='stats';
        prct=0;
        regions=[1 2 3 9 10    11 12 13 19 20];
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
                [dir,name,ext,num]=spm_fileparts(raw_data(i,:));
                Raw{i}=fullfile(baseDir,'imaging_data',subj_name{s},[name ext num]);
            end;
            V=spm_vol(char(Raw));
            
            % Loop over the possile regions 
            for r=regions
                indx=(Act.SN==s & Act.regNum==r);
                M=mean(Act.handAct(indx,1:2),2)./sqrt(Act.ResMs(indx,:));
                [~,i]=sort(M,1,'descend');
                R{r}.data=R{r}.data(i(1:min(100,length(i))),:);
                Ya = region_getdata(V,R{r});  % Data is N x P
                P=size(Ya,2);
                    
                % Make the two contrast matrices:
                C1=indicatorMatrix('allpairs_p',(double(D.hand==0 & D.stimType==0) .* D.digit));
                C2=indicatorMatrix('allpairs_p',(double(D.hand==1 & D.stimType==0) .* D.digit));
                C3=indicatorMatrix('allpairs_p',(double(D.hand==1 & D.stimType==1) .* D.digit));
                    
                [dist,~,~,conw]=distance_ldt_spm(Ya,SPM,[C1 C2 C3]',D.run','C0',indicatorMatrix('identity',D.run));
                c=mean(conw)/P; 
                
                % end;
                vec=[1;1;1];
                S.SN=s*vec;
                S.region=r*vec;
                S.hand=[1;2;2];
                S.stimtype=[1;1;2];
                S.group=subj_group(s)*vec;
                
                S.ldt=[dist(1:10);dist(11:20);dist(21:30)];
                S.ldc=[c(1:10);c(11:20);c(21:30)];
                T=addstruct(T,S);
                fprintf('%d %d\n',s,r);
            end;
        end;
                
        Side=[ones(1,10) ones(1,10)*2]; 
        Type=[1:10 1:10];
        T.regSide=Side(T.region)';
        T.regType=Type(T.region)';
        
        varargout={T};
    case 'accstats_subspace'
        side={'L','R'};
        
        y=varargin{1};
        D=varargin{2};
        T=varargin{3};
        sizeSubspace=varargin{4};
        c=D.digit';
        h=D.hand';
        [V,X]=lmva_randomsubspace(y,@df1_acc,'params',{D.digit',D.hand',D.stimType',D.run'},'num_iter',200,'size_subspace',sizeSubspace);
        
        if (~isempty(X))
            accS=mean(X.acc);
            E.accLm=accS(1)';
            E.accRm=accS(2)';
            E.accRs=accS(3)';
            E.sizeSubspace=sizeSubspace;
        else
            E.accLm=[];
            
        end;
        varargout={E};
    case 'distance+'%____________DF______________Compute difference between the patterns + distance from zero
        % T=df1_imana('ROI_stats','selection','none','fcn','distance+','regions',[6 12]);        %
        regular=0.001;
        y=varargin{1};
        [P,N]=size(y);
        D=varargin{2};
        R=varargin{3};
        F=[];
        for s=0:1
            for h=0:1
                j = find(D.hand==h & D.stimType==s);  % select datapoints in this class
                if (~isempty(j))
                    vec=[1:10]';
                    [U,ss,V]=svd(y(:,j),0);
                    E.dist_mahal_all=distance_mahalanobis(ss*V',D.digit(j)')';
                    E.dist_all=distance_crossval_pairs(ss*V',D.digit(j)',D.run(j)')';
                    [X,M]=lmva_randomsubspace(y(:,j),@distance_crossval_pairs,...
                        'params',{D.digit(j)',D.run(j)'},...
                        'num_iter',400,...
                        'size_subspace',100);
                    E.dist_sub=mean(M.acc);
                    E.hand=h;
                    E.stimtype=s;
                    F=addstruct(F,E);
                end;
            end;
        end;
        varargout={F};
    case 'distance_compare'%_____DF______________
        
        color={[0 0 1],[1 0 0],[0 1 0]};
        normalize=0;
        var='dist_sub';
        
        vararginoptions(varargin,{'normalize','var'});
        
        field_dist={'dist_mahal_all','dist_all','dist_sub','dist_ldaT'};
        field_oth={'SN','hand','region','regSide','regType'};
        
        
        TD=load(fullfile(regDir,'reg_distance_comp.mat'));
        TDr=load(fullfile(regDir,'reg_distance_ldaT.mat'));
        TD.dist_ldaT=TDr.ldaT;
        
        T1=load(fullfile('/Users/jdiedrichsen/Projects/FingerPattern/tendigit1/RegionOfInterest','reg_distance_comp.mat'));
        T1r=load(fullfile('/Users/jdiedrichsen/Projects/FingerPattern/tendigit1/RegionOfInterest','reg_distance_ldaT.mat'));
        T1.dist_ldaT=T1r.ldaT;
        
        % Select contralateral Motor activity
        TD.hand=TD.hand+1; % Make hands 1,2: DO THIS IN GENERAL??
        TD=getrow(TD,TD.hand~=TD.regSide & TD.stimtype==0 & TDr.regType==6);
        T1=getrow(T1,T1.hand~=T1.regSide & T1.regType==6);
        
        % wrangle the data into a new format: In general it may be useful
        % to have this type of structure
        for i=1:length(field_dist)
            x=TD.(field_dist{i}); % Unnormalized
            if (normalize)
                x=bsxfun(@rdivide,x,mean(x,2));
            end;
            TDa.(field_dist{i})=x (:);
            
            x=T1.(field_dist{i}); % Unnormalized
            if (normalize)
                x=bsxfun(@rdivide,x,mean(x,2));
            end;
            T1a.(field_dist{i})=x(:);
        end;
        for i=1:length(field_oth)
            TDa.(field_oth{i})=repmat(TD.(field_oth{i}),10,1);
            T1a.(field_oth{i})=repmat(T1.(field_oth{i}),10,1);
        end;
        TDa.pair=kron([1:10]',ones(length(TD.SN),1));
        T1a.pair=kron([1:10]',ones(length(T1.SN),1));
        
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
        varargout={T};
    case 'mva_decomp_figure' % This is Figure X
        % T= getrow(D, D.region== 1 |D.region== 2 | D.region== 6 | D.region== 16 |D.region== 17 |D.region== 18)
        % df1_imana('mva_decomp_figure', T )
        % df1_imana('mva_decomp_figure', getrow(T, T.SN==1))
        T=varargin{1};
        %T=getrow(T,T.reg~=2);
        %xlabel= regname;
        %xlabel1= regname;
        %figure(1);
        % Sepate estimate for sensory and motor
        subplot(4,2,1);barplot(T.region,[T.Dvar_cs T.Dvar_cm]);ylabel('var_c');
        subplot(4,2,2);barplot(T.region,T.Dcov_c./(sqrt(T.Dvar_cs).*sqrt(T.Dvar_cm)));ylabel('cov_c/var_c');
        
        subplot(4,2,3);barplot(T.region,[T.Dvar_fs T.Dvar_fm]);ylabel('var_f');
        subplot(4,2,4);barplot(T.region,T.Dcov_f./(sqrt(T.Dvar_fs).*sqrt(T.Dvar_fm)));ylabel('cov_f/var_f');
        subplot(4,2,5);barplot(T.region,[T.Dvar_bs T.Dvar_bm]);ylabel('var_b');
        subplot(4,2,6);barplot(T.region,T.Dcov_b./(sqrt(T.Dvar_bs).*sqrt(T.Dvar_bm)));ylabel('cov_b/var_b');
        subplot(4,2,7);barplot(T.region,T.Dvar_e);ylabel('var_e');
        %         for i=1:7
        %             subplot(4,2,i);
        %             if (mod(i,2)==0 )
        %                 set(gca,'YLim',[-0.2 0.7]);
        %                 %set(gca,'XTickLabel',xlabel);
        %             else
        %                 %set(gca,'XTickLabel',xlabel1);
        %             end;
        %         end;
        %         subplot(4,2,8);
        %         set(gcf,'PaperPosition',[2 2 11 16]/2.54);
        %         wysiwyg;
        %         figure(2);
        %         subplot(2,2,2);
        %         barplot(T.region,[T.Tcorr_fs T.Tcorr_fd],'leg',{'same','diff'});
        %         set(gca,'XTickLabel',xlabel1);
        %         ylabel('correlation');
        %         set(gca,'YLim',[-0.1 0.4]);
        %
        %         subplot(2,2,3);
        %         barplot(T.region,[T.Tcorr_sc_s T.Tcorr_sc_m]);
        %         ylabel('correlation');
        %         set(gca,'XTickLabel',xlabel1);
        %
        %         subplot(2,2,4);
        %         barplot(T.region,[T.Tcorr_sc_sb T.Tcorr_sc_db]);
        %         ylabel('correlation');
        %         set(gca,'XTickLabel',xlabel1);
        %         set(gcf,'PaperPosition',[2 2 11 8]/2.54);
        %         wysiwyg;
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
    case 'data_ROI_compare1' % Compare different ways of getting the betas - just make sure we are replicating results from sl1
        S=[];
        for s   = [1:13];        % Subject
            A1=load(fullfile(baseDir,'glm_firstlevel_1',subj_name{s},'SPM.mat'));
            A2=load(fullfile(baseDir,'glm_firstlevel_2',subj_name{s},'SPM.mat'));
            B1=load(fullfile(baseDir,'glm_firstlevel_1',subj_name{s},'SPM_info.mat'));
            B2=load(fullfile(baseDir,'glm_firstlevel_2',subj_name{s},'SPM_info.mat'));
            % B.SPM=spmj_move_rawdata(B.SPM,fullfile(baseDir, 'imaging_data',subj_name{s}, condition{sc}));
            load(fullfile(regDir,[subj_name{s} '_regions.mat']));
            
            for reg = [1 2 3 9 10 11];       % Regions
                % Redo the analysis directly from the raw data
                D=region_getdata(A2.SPM.xY.VY,{R{reg}});
                y_raw=D{1};
                P=size(y_raw,2); % Number of voxels
                
                % Filter the data with W from the first GLM. Also filter the
                % Designmatrix for the second GLM using exactly the same W
                y_filt = spm_filter(A1.SPM.xX.K,y_raw);
                A1f.xKXs   = spm_sp('Set',spm_filter(A1.SPM.xX.K,A1.SPM.xX.X));       % KWX
                A1f.xKXs.X = full(A1f.xKXs.X);
                A1f.pKX    = spm_sp('x-',A1f.xKXs);                        % projector
                
                % Do the raw version
                A1n.xKXs   = spm_sp('Set',A1.SPM.xX.X);       % KWX
                A1n.xKXs.X = full(A1n.xKXs.X);
                A1n.pKX    = spm_sp('x-',A1n.xKXs);                        % projector
                
                A2f.xKXs   = spm_sp('Set',spm_filter(A1.SPM.xX.K,A2.SPM.xX.X));       % KWX
                A2f.xKXs.X = full(A2f.xKXs.X);
                A2f.pKX    = spm_sp('x-',A2f.xKXs);                        % projector
                
                A2n.xKXs   = spm_sp('Set',A2.SPM.xX.X);       % KWX
                A2n.xKXs.X = full(A2n.xKXs.X);
                A2n.pKX    = spm_sp('x-',A2n.xKXs);                        % projector
                
                b1f = A1f.pKX*y_filt;                              %-Parameter estimates of the 4 reg GLM
                b1n = A1n.pKX*y_raw;
                b2f = A2f.pKX*y_filt;
                b2n = A2n.pKX*y_raw;                              %-Parameter estimates of the 4 reg GLM
                
                % Set up the design matrix on the second level (Z) for
                % 1 regressor version
                B1.type=B1.hand+B1.stimType+1;                      %
                Z1a=indicatorMatrix('hierarchicalI',[B1.type B1.digit]);
                Z1a=bsxfun(@minus,Z1a,mean(Z1a)); 
                Z1i=indicatorMatrix('identity',B1.run);
                partition1=B1.run;
                
                % Set up the design matrix on the second level (Z) for
                % 2 regressor version
                B2.type=B2.hand+B2.stimType+1;             % Right movement
                Z2a=indicatorMatrix('hierarchicalI',[B2.type B2.digit]);
                Z2a=bsxfun(@minus,Z2a,mean(Z2a)); 
                Z2i=indicatorMatrix('identity',B2.run);
                partition2=B2.run;
                
                partitionT = [];
                for i=1:length(A2.SPM.nscan)
                    partitionT = [partitionT; ones(A2.SPM.nscan(i),1)*i];
                end;
                type = [ones(1,5) ones(1,5)*2 ones(1,5)*3];
                
                % Do random subspace approach do get more reliable
                % measures
                for n=1:100
                    voxind=sample_wor([1:P],80,1);
                    T.SN=s;
                    T.regNum= reg;
                    T.acc1  = cnfn_crossval_regress4(b1f(1:120,voxind),Z1a,[],Z1i,partition1,'split',type);            % firstlevel_4
                    T.acc2  = cnfn_crossval_regress4(b1n(1:120,voxind),Z1a,[],Z1i,partition1,'split',type);            % firstlevel_4
                    T.acc3  = cnfn_crossval_regress4(b2f(1:240,voxind),Z2a,[],Z2i,partition2,'split',type);            % firstlevel_4
                    T.acc4  = cnfn_crossval_regress4(b2n(1:240,voxind),Z2a,[],Z2i,partition2,'split',type);            % firstlevel_4
                    T.acc5  = cnfn_crossval_regress6(y_filt(:,voxind),Z2a,Z2i,partitionT,A2f.xKXs.X,'split',type);             % Straight from time series
                    T.acc6  = cnfn_crossval_regress6(y_raw(:,voxind),Z2a,Z2i,partitionT,A2n.xKXs.X,'split',type);             % Straight from time series
                    S=addstruct(S,T);
                end;
                fprintf('%d %d\n',s,reg);
            end;
        end;
        
        T=S; 
        for i=1:3
            subplot(3,1,i);
            barplot([],[T.acc1(:,i) T.acc2(:,i) T.acc3(:,i) T.acc4(:,i) T.acc5(:,i) T.acc6(:,i) ]);
            set(gca,'XTickLabel',{'1f','1n','2f','2n','rf','rn'});
        end;
        
        varargout={S};
    case 'data_ROI_compare2' % Try different amounts of high-pass filtering
        S=[];
        for s   = [1:5];        % Subject
            A1=load(fullfile(baseDir,'glm_firstlevel_1',subj_name{s},'SPM.mat'));
            A2=load(fullfile(baseDir,'glm_firstlevel_2',subj_name{s},'SPM.mat'));
            B1=load(fullfile(baseDir,'glm_firstlevel_1',subj_name{s},'SPM_info.mat'));
            B2=load(fullfile(baseDir,'glm_firstlevel_2',subj_name{s},'SPM_info.mat'));
            % B.SPM=spmj_move_rawdata(B.SPM,fullfile(baseDir, 'imaging_data',subj_name{s}, condition{sc}));
            load(fullfile(regDir,[subj_name{s} '_regions.mat']));
            
            for reg = [2];       % Regions
                % Redo the analysis directly from the raw data
                D=region_getdata(A2.SPM.xY.VY,{R{reg}});
                y_raw=D{1};
                P=size(y_raw,2); % Number of voxels
                
                % Filter the data with W from the first GLM. Also filter the
                % Designmatrix for the second GLM using exactly the same W
                
                X{1} = A2.SPM.xX.X;       % KWX
                for i=1:4
                    X{i+1}=X{i};
                    for sess=1:8
                        X{i+1}(A1.SPM.xX.K(sess).row,end+1)=A1.SPM.xX.K(sess).X0(:,i);
                    end;
                end;
                
                
                % Set up the design matrix on the second level (Z) for
                % 1 regressor version
                B1.type=B1.hand+B1.stimType+1;                      %
                Z1a=indicatorMatrix('hierarchicalI',[B1.type B1.digit]);
                Z1i=indicatorMatrix('identity',B1.run);
                partition1=B1.run;
                
                % Set up the design matrix on the second level (Z) for
                % 2 regressor version
                B2.type=B2.hand+B2.stimType+1;             % Right movement
                Z2a=indicatorMatrix('hierarchicalI',[B2.type B2.digit]);
                Z2i=indicatorMatrix('identity',B2.run);
                partition2=B2.run;
                
                partitionT = [];
                for i=1:length(A2.SPM.nscan)
                    partitionT = [partitionT; ones(A2.SPM.nscan(i),1)*i];
                end;
                type = [ones(1,5) ones(1,5)*2 ones(1,5)*3];
                
                % Do random subspace approach do get more reliable
                % measures
                for n=1:100
                    voxind=sample_wor([1:P],80,1);
                    T.SN=s;
                    T.regNum= reg;
                    for i=1:5
                        T.numTrend = i-1;
                        %                             b2f=pinv(X{i})*y_raw(:,voxind);
                        %                             T.acc1  = cnfn_crossval_regress4(b2f(1:240,voxind),Z2a,[],Z2i,partition2,'split',type);            % firstlevel_4
                        T.acc2  = cnfn_crossval_regress6(y_raw(:,voxind),Z2a,Z2i,partitionT,X{i},'split',type);             % Straight from time series
                        S=addstruct(S,T);
                    end;
                    S=addstruct(S,T);
                end;
                fprintf('%d %d\n',s,reg);
            end;
        end;
        
        
        set(gca,'XTickLabel',{'1f','1n','4f','4n','rf','rn'});
        varargout={S};
    case 'data_ROI_compare3' % Weighting factor of OLS - weighting seems to give small advatange
        S=[];
        for s   = [1:5];        % Subject
            A1=load(fullfile(baseDir,'glm_firstlevel_1',subj_name{s},'SPM.mat'));
            A2=load(fullfile(baseDir,'glm_firstlevel_2',subj_name{s},'SPM.mat'));
            B2=load(fullfile(baseDir,'glm_firstlevel_2',subj_name{s},'SPM_info.mat'));
            % B.SPM=spmj_move_rawdata(B.SPM,fullfile(baseDir, 'imaging_data',subj_name{s}, condition{sc}));
            load(fullfile(regDir,[subj_name{s} '_regions.mat']));
            
            for reg = [2];       % Regions
                % Redo the analysis directly from the raw data
                D=region_getdata(A2.SPM.xY.VY,{R{reg}});
                y_raw=D{1};
                P=size(y_raw,2); % Number of voxels
                
                % Filter the data with W from the first GLM. Also filter the
                % Designmatrix for the second GLM using exactly the same W
                y_w = A1.SPM.xX.W*y_raw;
                
                A2w.xKXs   = spm_sp('Set',A1.SPM.xX.W*A2.SPM.xX.X);       % KWX
                A2w.xKXs.X = full(A2w.xKXs.X);
                A2w.pKX    = spm_sp('x-',A2w.xKXs);                        % projector
                
                A2n.xKXs   = spm_sp('Set',A2.SPM.xX.X);       % KWX
                A2n.xKXs.X = full(A2n.xKXs.X);
                A2n.pKX    = spm_sp('x-',A2n.xKXs);                        % projector
                
                b2w = A2w.pKX*y_w;
                b2n = A2n.pKX*y_raw;                              %-Parameter estimates of the 4 reg GLM
                
                % Set up the design matrix on the second level (Z) for
                % 2 regressor version
                B2.type=B2.hand+B2.stimType+1;             % Right movement
                Z2a=indicatorMatrix('hierarchicalI',[B2.type B2.digit]);
                Z2a=bsxfun(@minus,Z2a,mean(Z2a)); 
                Z2i=indicatorMatrix('identity',B2.run);
                partition2=B2.run;
                 type = [ones(1,5) ones(1,5)*2 ones(1,5)*3];               
                % Do random subspace approach do get more reliable
                % measures
                for n=1:100
                    voxind=sample_wor([1:P],80,1);
                    T.SN=s;
                    T.regNum= reg;
                    T.acc1  = cnfn_crossval_regress4(b2w(1:240,voxind),Z2a,[],Z2i,partition2,'split',type);            % firstlevel_4
                    T.acc2  = cnfn_crossval_regress4(b2n(1:240,voxind),Z2a,[],Z2i,partition2,'split',type);            % firstlevel_4
                    S=addstruct(S,T);
                end;
                fprintf('%d %d\n',s,reg);
            end;
        end;
        for i=1:3
            subplot(3,1,i);
            barplot([],[S.acc1(:,i) S.acc2(:,i) ]);
            set(gca,'XTickLabel',{'1w','1n'});
        end;
        
        varargout={S};
    case 'force_scanner'
        sn=varargin{1};
        cd(fullfile(behaviourDir,subj_name{sn}))
        D=dload(['DF1_', subj_name{sn}, '.dat']);
        barplot(D.hand, D.maxForce, 'split', [D.stimType D.digit], 'subset', [D.announce==0& D.digit~=0], 'facecolor',[0.5 0.5 0.5]);
        set(gcf,'PaperPosition',[1 1 6 2]);
        wysiwyg;
        saveas(gcf, 'maxForce', 'jpg')
        
    case 'make_alldat_IN2b'     % case for Anna behavioural analysis
        S=[];
        for i=1:length(subj_name)
            Si = load(fullfile(behaviourEMGDir,'data',subj_name{i},sprintf('IN2b_%s.mat',subj_name{i})));
            Si.SN = zeros(length(Si.BN),1)+i;
            S=addstruct(S,Si);
        end;
        save(fullfile(behaviourEMGDir,'analysis','IN2b_alldat.mat'),'-struct','S');
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

