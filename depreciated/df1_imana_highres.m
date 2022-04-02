function varargout=df1_imana(what,varargin)
baseDir=        fullfile('/Users','tob', 'Projects','FingerPattern_dystonia');
% baseDir=        fullfile('/media','DATA', 'Projects',''SequenceLearning','FingerPattern_dystonia', 'FingerPattern_dystonia');
% baseDir=        fullfile('/Users/jdiedrichsen/Projects/FingerPattern_dystonia');
% baseDir=        fullfile('/Volumes/MotorControl/project/FingerPattern_dystonia);

behaviourDir=   fullfile(baseDir, 'Behavioural_data');%, 'analyze');
groupData=      fullfile(baseDir, 'group_data');
groupDir=       fullfile(baseDir, 'group_analysis');
anatomicalDir=  fullfile(baseDir, 'anatomicals');
freesurferDir=  fullfile(baseDir, 'surfaceFreesurfer');
caretDir=       fullfile(baseDir, 'surfaceCaret');
regDir=         fullfile(baseDir, 'RegionOfInterest');
glmName= {'glm_firstlevel_1'}; 

% atlasA={'i','x'};
% atlasname={'fsaverage','fsaverage_sym'}; 

subj_MTname={'MT02481','MT02481'};
subj_name={'s01-1','s01-2',};  % s01-1 -> old seq s01-2 -> new seq 

radius= 12;
numVox= 60;
run={'01','02','03','04','05','06','07'};
hem={'lh','rh'};
hemName={'LeftHem','RightHem'};
condition= {'', ''};


ROI_coords={[18 -59 -34],[-32 -33 62]; % right: cerebellum left: S1 SN 1
    [18 -59 -34],[-32 -33 62];}; %SN2

switch(what)
    case 'preprocess'
        %df1_imana('preprocess', 2)
        sn=varargin{1};
        for s=sn
            %             df1_imana('make_nii', s)
            %             %NOT DONE df1_imana('slice_timing',s);
            %             df1_imana('makefieldmap',s);
            %             df1_imana('make_realign_unwarp',s);
            %             df1_imana('move_images',s);
            %             df1_imana('meanimage_bias_correction',s);
            %             df1_imana('segmentation',s)
            %             figure(s); df1_imana('plot_movementparameters',s);
            %%need anatomical sxx_2_MT.nii at this point
            %             df1_imana('coreg', s)
            %             df1_imana('make_samealign',s);
            %             df1_imana('check_samealign',s);
            %             df1_imana('make_maskImage',s)
            
            %             df1_imana('make_glm_expl_mask',s)
            %             df1_imana('estimate_glm',s, 'glm_firstlevel_1')
            %             df1_imana('contrast',s)
            %             df1_imana('MVA_search',s)
            %             df1_imana('MVA_do',s)
            %             df1_imana('MVA_zValue',s,1)
            %       df1_imana('MNI_normalization_write',s, 1)
            %             df1_imana('surf_MVA_search',s, 1, 0); %LEFT
            %             df1_imana('surf_MVA_search',s, 2, 0); %RIGHT
            %             df1_imana('surf_MVA_do',s,1);
            %             df1_imana('surf_MVA_do',s,2);
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
        for i=1:7
            if strcmp(subj_name{sn}, 's01-2')
                spmj_tar2nii ([subj_MTname{sn} '.r', num2str(i) '.tar'], [subj_name{sn} '_run' run{i} '.nii']);
            else
                spmj_tar2nii ([subj_MTname{sn} '.r', num2str(i) '.tar'], [subj_name{sn} '_run' run{i} '.nii'], 'startTR', 4)
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
    case 'segmentation'%____________STEP 1.1___________________________
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
        for r= 1:7
            for i=1:79
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
        spmj_realign_unwarp(baseDir, subj_name{sn}, run, 1, 79,'prefix',prefix, 'subfolderRawdata','') %76
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

        copyfile(source,dest);
    case 'meanimage_bias_correction'%_________________________________
        % df1_imana('meanimage_bias_correction',1)
        sn=varargin{1};
        %prefix='ua';
        prefix='u';
        %full epi
        P{1}={fullfile(baseDir, 'imaging_data_raw',subj_name{sn}, [subj_name{sn},'_epi.nii'])};
        spm_biascorrect(P{1});
        source = fullfile(baseDir, 'imaging_data_raw',subj_name{sn},['b',subj_name{sn},'_epi.nii']);
        dest = fullfile(baseDir, 'imaging_data',subj_name{sn}, [subj_name{sn}, '_epi.nii']);
        copyfile(source,dest);
        %functional mean image
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
            x= dlmread (fullfile(baseDir, 'imaging_data',subj_name{sn}, ['rp_' subj_name{sn},'_run',run{j},'.txt']));
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
    case 'make_samealign'%__________STEP 7____________________________
        % df1_imana('make_samealign',1)
        %prefix='ra';
        %prefix='ua';
        prefix='u';
        sn=varargin{1};
        Q={};
        cd(fullfile(baseDir, 'imaging_data',subj_name{sn}));
        for r= 1:numel(run)
            Q{end+1} = [fullfile(baseDir, 'imaging_data',subj_name{sn}, [prefix, subj_name{sn},'_run',run{r},'.nii'])];
        end
        P{1} = fullfile(baseDir, 'imaging_data', subj_name{sn}, ['meanepi_' subj_name{sn} '.nii']);
        spmj_makesamealign_nifti(char(P),char(Q))
        %delete(char(P))
    case 'check_samealign'%_________STEP 8____________________________
        % df1_imana('check_samealign',1)
        %prefix='ua';
        prefix='u';
        sn=varargin{1};
        Q={};
        cd(fullfile(baseDir, 'imaging_data',subj_name{sn}));
        for r= 1:numel(run)
            for i=1:76
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
    case 'make_glm_1'%______________STEP 9_________________________
        %df1_imana('make_glm_1',1)
        % set threshold for the masking
        spm fmri
        subj=varargin{1};
        for sn=subj
            %----
            delay=  1*32;
            dur=    6%3.8;
            prefix= 'ua'; %'a' for new sequence
            numVolumes= 79; %76 for new seq
            num_dummys=3;
            timePerSlice= 85/1000; % time pro slice in sec
            D=dload( fullfile(behaviourDir,subj_name{sn}, ['DF1_',subj_name{sn},'.dat']));
            %D=getrow(D,D.day==7 & D.lastTrial==0); % Take scan data
            %sequences= unique(D.seqType);
            T=[];
            %correct for the number of dummy scans (not in nii file)
            %TR= D.startTR-num_dummys;
            Slice= D.startSlice-num_dummys*32;
            
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
                %for h=1:2 % hand
                for s= 1%:2 % stimType
                    for d=1:5 % digit
                        %idx_digit=find(D.BN==r & D.hand==h & D.stim==s & D.digit==c & D.announce==1);
                        idx_digit=find(D.BN==r  & D.stimType==s & D.digit==d & D.announce==1);
                        %J.sess(r).cond(end+1).name = ['B%dH%dC%dD%d',r,h,s,d];
                        J.sess(r).cond(end+1).name = sprintf('B%dS%dD%d',r,s,d);
                        J.sess(r).cond(end).onset = [Slice(idx_digit)- delay]*timePerSlice;
                        %----print information on screen
                        fprintf('run: %d \t stimType: %d \t digit: %d\n',r,s,d);
                        fprintf('onset time in seconds: %2.2fs\n',([Slice(idx_digit)- delay]*timePerSlice)/60);
                        %----
                        J.sess(r).cond(end).duration =  dur;
                        J.sess(r).cond(end).tmod = 0;
                        J.sess(r).cond(end).pmod = struct('name', {}, 'param', {}, 'poly', {});
                        S.sn=sn;
                        S.run=r;
                        S.numEvents=length(idx_digit);
                        %S.hand=h;
                        S.regType=1;  % Digit regressor
                        S.digit=d; % Sequences laled after its type
                        S.stimType=s;         % Sequences labled 1-4
                        T=addstruct(T,S);
                    end
                end;
                %end
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
            spm_jobman('run',matlabbatch);
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
        end;
        % varargout={J};
    case 'make_glm_1_1'%__high res sequence____STEP 9_________________________
        %df1_imana('make_glm_1_1',2)
        % set threshold for the masking
        spm fmri
        subj=varargin{1};
        for sn=subj
            global defaults;
            defaults.mask.thresh= 0.0008;
            %----
            delay=  2.5*32;
            dur=    3.6;
            prefix= 'u'; %'a' for new sequence
            numVolumes= 76; 
            num_dummys=0;
            timePerSlice= 70/1000; % time pro slice in sec
            D=dload( fullfile(behaviourDir,subj_name{sn}, ['DF1_',subj_name{sn},'.dat']));
            %D=getrow(D,D.day==7 & D.lastTrial==0); % Take scan data
            %sequences= unique(D.seqType);
            T=[];
            %correct for the number of dummy scans (not in nii file)
            %TR= D.startTR-num_dummys;
            Slice= D.startSlice-num_dummys*32;
            
            J.dir = {fullfile(baseDir,'glm_firstlevel_1', subj_name{sn})};
            if (~exist(J.dir{1},'dir'))
                mkdir(J.dir{1});
            end;
            
            J.timing.units = 'secs';% 'scans';%
            J.timing.RT = 2.72; %2.8;%
            J.timing.fmri_t = 16;
            J.timing.fmri_t0 = 1;
            
            for r=1:numel(run) %run
                for i=1:numVolumes
                    N{i} = fullfile(baseDir, 'imaging_data',subj_name{sn},[prefix,subj_name{sn},'_run',run{r},'.nii,',num2str(i)]); % ua uab 'ra'
                end;
                J.sess(r).scans= N;
                J.sess(r).cond=[];
                %for h=1:2 % hand
                for s= 1%:2 % stimType
                    for d=1:5 % digit
                        %idx_digit=find(D.BN==r & D.hand==h & D.stim==s & D.digit==c & D.announce==1);
                        idx_digit=find(D.BN==r  & D.stimType==s & D.digit==d & D.announce==1);
                        %J.sess(r).cond(end+1).name = ['B%dH%dC%dD%d',r,h,s,d];
                        J.sess(r).cond(end+1).name = sprintf('B%dS%dD%d',r,s,d);
                        J.sess(r).cond(end).onset = [Slice(idx_digit)- delay ]*timePerSlice;
                        %----print information on screen
                        fprintf('run: %d \t stimType: %d \t digit: %d\n',r,s,d);
                        fprintf('onset time in seconds: %2.2fs\n',([Slice(idx_digit)- delay]*timePerSlice)/60);
                        %----
                        J.sess(r).cond(end).duration =  dur;
                        J.sess(r).cond(end).tmod = 0;
                        J.sess(r).cond(end).pmod = struct('name', {}, 'param', {}, 'poly', {});
                        S.sn=sn;
                        S.run=r;
                        S.numEvents=length(idx_digit);
                        %S.hand=h;
                        S.regType=1;  % Digit regressor
                        S.digit=d; % Sequences laled after its type
                        S.stimType=s;         % Sequences labled 1-4
                        T=addstruct(T,S);
                    end
                end;
                %end
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
            spm_jobman('run',matlabbatch);
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
        end;
        % varargout={J};    
    case 'make_glm_2'%__________One regressor per run, but with extra regressors for error trials
        global defaults;
        defaults.mask.thresh= 0.7;
        prefix= 'ua';
        sn=varargin{1};
        num_dummys=3; %%% duration of mini-block
        D=load( fullfile(behaviourDir,'D.mat'));
        D=getrow(D,D.sn==sn & D.lastTrial==0);
        T=[];
        
        
        sequences= unique(D.seqType(D.seqType~=0));
        TR= D.startTR-num_dummys;
        
        J.dir = {fullfile(baseDir,'glm_firstlevel_2', subj_name{sn})};
        if (~exist(J.dir{1},'dir'))
            mkdir(J.dir{1});
        end;
        J.timing.units = 'scans';%'secs';
        J.timing.RT = 2.72;
        J.timing.fmri_t = 16;
        J.timing.fmri_t0 = 1;
        
        for r=1:numel(run)
            for i=1:156
                N{i} = fullfile(baseDir, 'imaging_data',subj_name{sn},[prefix,subj_name{sn},'_run',run{r},'.nii,',num2str(i)]); % ua uab 'ra'
            end;
            J.sess(r).scans= N;
            J.sess(r).cond=[];
            for h=1:2
                % Define Sequence regressors
                for c=1:4
                    idx=find(D.seqType==sequences(c)& D.BN==r & D.announce==0 & D.hand==h & D.myError==0);
                    if (length(idx)>0)
                        J.sess(r).cond(end+1).name = sprintf('C%dT%d', h,c);
                        J.sess(r).cond(end).onset = [TR(idx)-1];
                        J.sess(r).cond(end).duration =  0;
                        J.sess(r).cond(end).tmod = 0;
                        J.sess(r).cond(end).pmod = struct('name', {}, 'param', {}, 'poly', {});
                        % Make structure for contrast, and later analysis
                        S.sn=sn;
                        S.run=r;
                        S.numEvents=length(idx);
                        S.hand=h;
                        S.regType=1;  % Sequence regressor
                        S.seqType=sequences(c); % Sequences laled after its type
                        S.seqNum=c;         % Sequences labled 1-4
                        T=addstruct(T,S);
                    end;
                end;
                % Define error regressors
                idx=find( D.BN==r & D.announce==0 & D.hand==h & D.myError==1);
                if (length(idx)>0)
                    J.sess(r).cond(end+1).name = sprintf('C3T%d', h);
                    J.sess(r).cond(end).onset = [TR(idx)-1];
                    J.sess(r).cond(end).duration =  0;
                    J.sess(r).cond(end).tmod = 0;
                    J.sess(r).cond(end).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    % Make structure for contrast, and later analysis
                    S.sn=sn;
                    S.run=r;
                    S.numEvents=length(idx);
                    S.hand=h;
                    S.regType=2;  % Error regressor
                    S.seqType=NaN; % Sequences laled after its type
                    S.seqNum=NaN;         % Sequences labled 1-4
                    T=addstruct(T,S);
                end;
            end;
            
            
            
            J.sess(r).multi = {''};
            J.sess(r).regress = struct('name', {}, 'val', {});
            J.sess(r).multi_reg = {''};
            J.sess(r).hpf = 128;
        end;
        J.fact = struct('name', {}, 'levels', {});
        J.bases.hrf.derivs = [0 0];
        J.volt = 1;
        J.global = 'None';
        J.mask = {''};
        J.cvi =  'wls';
        matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec=J;
        spm_jobman('run',matlabbatch);
        save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
        
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
        % load(SPM)
        % spm_rwls_resstats(SPM)
    case 'hrf_getdata'%_______________________________________________
        % to check the model quality of the glm
        %df1_imana('hrf_getdata',1)
        sn=varargin{1};
        pre=4; post=12;
        for s=sn
            T=[];
            cd(fullfile(baseDir,'glm_firstlevel_1', subj_name{s}));
            load SPM;
            for r=1:size(ROI_coords,2)
                R{r}=region('sphere',ROI_coords{s,r},4);
            end
            [y_raw, y_adj, y_hat, y_res,B] = region_getts(SPM,R);
            
            D=spmj_get_ons_struct(SPM);
            ons= D.ons/2.72;
            for r=1:size(y_raw,2)
                for i=1:size(D.block,1);
                    D.y_adj(i,:)=cut(y_adj(:,r),pre,round(ons(i))-1,post,'padding','nan')';
                    D.y_hat(i,:)=cut(y_hat(:,r),pre,round(ons(i))-1,post,'padding','nan')';
                    D.y_res(i,:)=cut(y_res(:,r),pre,round(ons(i))-1,post,'padding','nan')';
                end;
                D.region=ones(size(D.event,1),1)*r;
                T=addstruct(T,D);
            end;
            save ROI_M1_IRF.mat T
        end;
    case 'hrf_plot'%__________________________________________________
        %df1_imana('hrf_plot',1)
        sn=varargin{1}; 
        cd(fullfile(baseDir,'glm_firstlevel_1', subj_name{sn}));
        load ROI_M1_IRF;
        figure
        for r=1:size(ROI_coords,2)
            subplot(size(ROI_coords,2),1,r);
            traceplot([-4:12],T.y_adj,'errorfcn','stderr','subset',T.region==r);
            hold on;
            traceplot([-4:12],T.y_hat,'subset',T.region==r,'linestyle',':');
            %             hold on;
            %             traceplot([-4:12],T.y_res,'subset',T.region==r,'linestyle','-');
            hold off;
            xlabel('TR');
            ylabel('activation');
            drawline(0);
        end;
    case 'contrast'%________________STREP 11 
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
            SPM.xCon(1)=spm_FcUtil('Set','t_seq_left', 'T', 'c',con',SPM.xX.xKXs);

%             %____F contrast between seq left
%             for h=1:2
%                 con=zeros(3,size(SPM.xX.X,2));
%                 con(1,T.seqNum==1 & T.hand==h)=1;con(1,T.seqNum==2 & T.hand==h)=-1;
%                 con(2,T.seqNum==2 & T.hand==h)=1;con(2,T.seqNum==3 & T.hand==h)=-1;
%                 con(3,T.seqNum==3 & T.hand==h)=1;con(3,T.seqNum==4 & T.hand==h)=-1;
%                 if (any(sum(con,2)~=0))
%                     keyboard;
%                 end;
%                 SPM.xCon(h)=spm_FcUtil('Set','F_seq_left', 'F', 'c',con',SPM.xX.xKXs);
%             end;
%             
%             %_____t contrast movement against rest left
%             for h=1:2
%                 con=zeros(1,size(SPM.xX.X,2));
%                 con(T.hand==h & T.regType==1)=1;
%                 con=con/sum(con);
%                 SPM.xCon(end+1)=spm_FcUtil('Set','t_seq_left', 'T', 'c',con',SPM.xX.xKXs);
%             end;
%             
%             
%             % Optional: t-contrast for each sequence / hand against rest 
%             for h=1:2
%                 for c=1:4 
%                     con=zeros(1,size(SPM.xX.X,2));
%                     con(T.hand==h & T.seqNum==c)=1;
%                     con=con/sum(con);
%                     SPM.xCon(end+1)=spm_FcUtil('Set',sprintf('seq %h %d',con), 'T', 'c',con',SPM.xX.xKXs);
%                 end;
%             end;
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
    case 'MVA_search'%______________STEP  12___________________________
        % df1_imana('MVA_search',1, 60)
        sn=varargin{1};
        numVox= varargin{2};
        for s=sn
            cd(fullfile(baseDir,'glm_firstlevel_1', subj_name{s}));
            V=spm_vol('rmask_noSkull.nii');   %V=spm_vol(['mask.img']);
            X=spm_read_vols(V);
            [vox_i,vox_j,vox_k]=ind2sub(size(X),find(X~=0));
            vox=[vox_i,vox_j,vox_k];
            [LI,voxmin,voxmax,n]=lmva_voxelselection(vox',vox',[radius, numVox],V.mat,V.dim,[],'mva_numvox.nii');
            save (sprintf('volsearch_%i.mat',numVox), 'vox', 'LI', 'voxmin', 'voxmax', 'n')
        end
        %___STEP 10______________________________________
    case 'MVA_do'%____________________________________________________
        % df1_imana('MVA_do',1, 60)
        sn=varargin{1};
        numVox= varargin{2};
        workDir= fullfile(baseDir, ['glm_firstlevel_1']);
        ldafunction= @df1_calcAcc_pilot;
        for s=sn
            cd(fullfile(workDir, subj_name{s}));
            T= load(fullfile(baseDir,'glm_firstlevel_1', subj_name{sn},'SPM_info.mat'));
            load('SPM');
            %define beta images
            beta_images=  {SPM.Vbeta(SPM.xX.iC).fname}';
            %left hand
            lda_names= {fullfile(workDir, subj_name{s}, ['lda_',num2str(numVox),'.nii']),...
                fullfile(workDir, subj_name{s},['lda_numberVoxel',num2str(numVox),'.nii'])};
            %----do LDA
            lmva_spm(sprintf('volsearch_%i.mat',numVox),beta_images,lda_names,ldafunction,'params',{T.run,0,T.stimType,T.digit})%T.hand
        end
    case 'MVA_search_regress'%________________________________________
        % df1_imana('MVA_search_regress',1)
        sn=varargin{1};
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            V=spm_vol('mask_noSkull.nii');   %V=spm_vol(['mask.img']);
            X=spm_read_vols(V);
            [vox_i,vox_j,vox_k]=ind2sub(size(X),find(X~=0));
            vox=[vox_i,vox_j,vox_k];
            [LI,voxmin,voxmax,n]=lmva_voxelselection(vox',vox',[radius, 160],V.mat,V.dim,[],'mva_numvox.nii');
            save volsearch.mat vox LI voxmin voxmax n
        end
    case 'CHECK_MVA_do_regress'%______________________________________
        % df1_imana('MVA_do_regress',1, 1, 'MT')
        sn=varargin{1}; seqCondition=varargin{2}; parameter=varargin{3}
        ldafunction= @sl1_calculate_accuracy_4reg_regress;
        %----load behavioural data
        load (fullfile(behaviourDir,'D.mat'));
        for s=sn
            T_behav= tapply(D,{'bn', 'bn_seq','seqType'},{parameter,'median','name','params'},'subset',D.phase==3 & D.subj==s & D.seqType>0 & D.trainSeq==seqCondition);
            %----define the images
            for i=1:128 %----get the beta images
                beta_resliced{i}= fullfile(resliceDir,subj_name{s}, [conditions{seqCondition+1}(1),'_beta_',num2str(i),'.nii']);
            end
            %----name for lda nifti
            lda_names= {fullfile(resliceDir, subj_name{s},  [conditions{seqCondition+1}(1),'_lda_vox_',num2str(radius),'_reg',parameter, '.nii']),...
                fullfile(resliceDir, subj_name{s}, [conditions{seqCondition+1}(1),'_lda_vox_',num2str(radius),'_reg',parameter, '_numberVoxel','.nii'])};
            %----design matrix
            X= [ones(size(T_behav.params)),T_behav.params];
            %----do LDA
            cd(fullfile(resliceDir,subj_name{s}));
            lmva_spm('volsearch.mat',beta_resliced,lda_names,ldafunction, 'params', {X})
        end
    case 'MVA_zValue'%________________________________________________
        %df1_imana('MVA_zValue',1:3,1)
        sn=varargin{1}; 
        glmType=varargin{2}
        if glmType==1
            workDir= [glmDir, '_1'];
            numRegressors= 32;
        elseif glmType==4
            workDir= fullfile(baseDir, 'glm_firstlevel_3');
            numRegressors= 32*3;
        end
        mu= 1/4*100;

        sigma= sqrt(1/4 * 3/4 * 1/numRegressors)*100;
        images= {['lda_L_',num2str(numVox),'.nii'], ['lda_R_',num2str(numVox),'.nii']}
        outimages= {['zlda_L_',num2str(numVox),'.nii'], ['zlda_R_',num2str(numVox),'.nii']}
        
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
    case 'regres_resliceImages'
        % df1_imana('regres_resliceImages',1)
        %define imput
        sn=varargin{1}; glmNum=varargin{2};
        for s= sn
            if ~exist(fullfile([resliceDir, '_', num2str(glmNum)],subj_name{s}))
                mkdir(fullfile([resliceDir, '_', num2str(glmNum)],subj_name{s}));
            end
            for c= 1:size(condition,2)
                %----do the reslice for mask
                spmj_reslice_vol(fullfile([glmDir, '_', num2str(glmNum)] , subj_name{s}, condition{c},   'mask_noskull.nii'), [100 100 65], spm_matrix([-100 -120 -40 0 0 0 2 2 2]) , fullfile([resliceDir, '_', num2str(glmNum)], subj_name{s}, [condition{c}(1),'_mask_noskull.nii']));
                %----do the reslice for con0006; movement against rest
                spmj_reslice_vol(fullfile([glmDir, '_', num2str(glmNum)], subj_name{s}, condition{c},   'con_0006.img'), [100 100 65], spm_matrix([-100 -120 -40 0 0 0 2 2 2]) , fullfile([resliceDir, '_', num2str(glmNum)], subj_name{s}, [condition{c}(1),'_con_0006.nii']));
                %----do the reslice for spmT_0006; movement against rest
                spmj_reslice_vol(fullfile([glmDir, '_', num2str(glmNum)], subj_name{s}, condition{c},   'spmT_0006.img'), [100 100 65], spm_matrix([-100 -120 -40 0 0 0 2 2 2]) , fullfile([resliceDir, '_', num2str(glmNum)], subj_name{s}, [condition{c}(1),'_spmT_0006.nii']));
                %----do the reslice for spmF_0001; movement against rest
                spmj_reslice_vol(fullfile([glmDir, '_', num2str(glmNum)], subj_name{s}, condition{c},   'spmF_0001.img'), [100 100 65], spm_matrix([-100 -120 -40 0 0 0 2 2 2]) , fullfile([resliceDir, '_', num2str(glmNum)], subj_name{s}, [condition{c}(1),'_spmF_0001.nii']));
                %----load SPM structure
                load(fullfile([glmDir, '_', num2str(glmNum)], subj_name{s}, condition{c},'SPM'));
                %----load beta images
                V= spm_vol(SPM.Vbeta(SPM.xX.iC));
                %----loop over beta images
                for i=SPM.xX.iC
                    %----do the reslice
                    spmj_reslice_vol(fullfile([glmDir, '_', num2str(glmNum)], subj_name{s}, condition{c},V(i).fname), [100 100 65], spm_matrix([-100 -120 -40 0 0 0 2 2 2]) , fullfile([resliceDir, '_', num2str(glmNum)],subj_name{s}, [condition{c}(1),'_beta_',num2str(i),'.nii']));
                end
            end
            %----calculate the AND image from the mask images over the conditions
            spm_imcalc_ui({fullfile([resliceDir, '_', num2str(glmNum)],subj_name{s},[condition{1}(1),'_mask_noskull.nii']), fullfile([resliceDir, '_', num2str(glmNum)],subj_name{s},[condition{2}(1),'_mask_noskull.nii'])}, ...
                fullfile([resliceDir, '_', num2str(glmNum)],subj_name{s},['mask_noskull.nii']), 'and(i1, i2)');
        end
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
        %df1_imana('MNI_normalization_write',1:3, 1)
        sn=varargin{1}; glmType=varargin{2}
        if glmType==1
            workDir= fullfile(baseDir, 'glm_firstlevel_1');
            outPrefix=[];
            images={'con_0004.img','con_0005.img','con_0006.img','rmask_noSkull.nii', 'zlda_L_160.nii', 'zlda_R_160.nii'} %
        elseif glmType==2
            workDir= fullfile(baseDir, 'glm_firstlevel_3');
            outPrefix='_3reg';
            images={''} %
        end
        for s=sn
            if s==8
                defor= fullfile(anatomicalDir, subj_name{s},'ana2', [subj_name{s}, '_2_R1_seg_sn.mat']);
            else
                defor= fullfile(anatomicalDir, subj_name{s},'ana2', [subj_name{s}, '_2_MT_seg_sn.mat']);
            end
            %defor= fullfile(anatomicalDir, subj_name{s},'ana2', [subj_name{s}, '_2_T1w_seg_sn.mat']);
            for j=1:numel(images)
                [dir,name,ext]=spm_fileparts(images{j});
                sn_images{j}= fullfile(workDir,subj_name{s},images{j});
                
                out_images{j}= fullfile(groupData,[name, '_' subj_name{s}, outPrefix,  '.nii']);
            end
            spmj_normalization_write(defor, sn_images,'outimages',out_images);
            
            % Do the anatomical
            sn_images={}; out_images={};
            sn_images{1}=fullfile(baseDir, 'anatomicals', subj_name{s},'ana2', [subj_name{s}, '_2_MT.nii']);
            out_images{1}=fullfile(groupData,[subj_name{s}, '_2_MT.nii']);
            spmj_normalization_write(defor, sn_images,'outimages',out_images);
        end;
    case 'MNI_mask_anatomical'%__________________________________________________
        %df1_imana('MNI_mask_anatomical', 1:3)
        sn=varargin{1};
        cd(groupData);
        for s=sn
            images{s}= fullfile(groupData,['rmask_noSkull_',subj_name{s},  '.nii']);
            ana_images{s}= fullfile(groupData,[subj_name{s}, '_2_MT.nii']);
        end
        spmj_imcalc_mtx(images,'mask_avrg.nii','nanmean(X)');
        spmj_imcalc_mtx(ana_images,'avrg_2_MT.nii','nanmean(X)');
        %spmj_imcalc_mtx([],'mask_avrg.nii','nanmean(X)');
        spmj_imcalc_mtx('mask_avrg.nii','mask_thres.nii','X>0.4');
        %     case 'MNI_train-untrain'
        %         %df1_imana('MNI_train-untrain',[1:8],'_4reg')
        %         sn=varargin{1}; prefix= varargin{2};
        %         images= {'zlda_radius12_', 'con_0006_'}
        %         nam=[];
        %         for i= 1:length(images)
        %             for s=sn
        %                 nam{1} = fullfile(groupData, [images{i}, condition{1},'_',subj_name{s},prefix, '.nii,1']);
        %                 nam{2} = fullfile(groupData, [images{i}, condition{2},'_',subj_name{s},prefix, '.nii,1']);
        %                 out= fullfile(groupData, [images{i}, 'trained-untrained','_',subj_name{s},prefix, '.nii,1']);
        %                 spm_imcalc_ui(nam,out,'i1-i2')
        %             end
        %         end
    case 'MNI_smooth'%________________________________________________
        %df1_imana('MNI_smooth',1:3)
        sn=varargin{1}
        images= {'zlda_L_160', 'zlda_R_160','con_0004','con_0005','con_0006'}
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
    case 'surf_ROI_coord'%__________________________________________
        % coord= df1_imana('coord_roiData_surface','subj', 's05', 'centerNode', 26804, 'Vdim', V.dim, 'Vin', V, 'LI', n2v)
        vararginoptions(varargin,{'subj', 'centerNode', 'Vdim', 'Vin', 'LI'});
        cd(fullfile(dataDir,'glm_firstlevel', subj));
        for i=1:numel(centerNode)
            [I,J,K]=ind2sub(Vdim, LI{centerNode});
            R= zeros(Vdim);
            for i=1:length(Vin)
                R(LI{centerNode})=1;%spm_sample_vol(Vin(i),double(I),double(J),double(K),0);
            end
        end;
        Vin.fname= 'roiData_surface.nii'
        spm_write_vol(Vin, R)
        varargout={[I; J; K]'};
    case 'surf_xhemireg'
        sn=varargin{1};
        for i=sn 
            freesurfer_registerXhem({subj_name{i}},freesurferDir,'hemisphere',[1:2]);
        end; 
    case 'surf_map_ico'%______________________________________________
        % df1_imana('surf_map_ico',6:7)
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
    case 'surf_make_caret'%___________________________________________
        % df1_imana('surf_make_caret',6:7)
        sn=varargin{1};
        atlas=varargin{2}; 
        for i=sn
              caret_importfreesurfer([atlasA{atlas} subj_name{i}],freesurferDir,caretDir);
        end;
        %=========================================================================
        % Map functional volumes on surface over caret
        %=========================================================================
    case 'surf_MVA_search' %__________________________________________
        %df1_imana('surf_MVA_search',1, 1);
        
        sn=varargin{1};   
        hem=varargin{2};
        atlas=varargin{3};
        refDir= fullfile(baseDir,glmName{1});
        
        for s=sn
            for h=hem
                if h== 1
                    caret_subjDIR = fullfile(caretDir,[atlasA{atlas},subj_name{s}],'LeftHem');
                    coord_pial= caret_load(fullfile(caret_subjDIR, 'lh.PIAL.coord'));
                    coord_white= caret_load(fullfile(caret_subjDIR, 'lh.WHITE.coord'));
                    topo= caret_load(fullfile(caret_subjDIR, 'lh.CLOSED.topo'));
                elseif h==2
                    caret_subjDIR = fullfile(caretDir,[atlasA{atlas},subj_name{s}],'RightHem');
                    coord_pial= caret_load(fullfile(caret_subjDIR, 'rh.PIAL.coord'));
                    coord_white= caret_load(fullfile(caret_subjDIR, 'rh.WHITE.coord'));
                    topo= caret_load(fullfile(caret_subjDIR, 'rh.CLOSED.topo'));
                else
                    disp('wrong hemisphere! Should be lh or rh ');
                end
                %get surface roi
                refV= spm_vol(fullfile(refDir, subj_name{s},'mask.img'));
                epiInfo.mat= refV.mat;
                epiInfo.dim= refV.dim;
                epiInfo.mask=spm_read_vols(refV);
                node_range=1:size(coord_white.data,1);
                [LI,voxmin,voxmax,vORr]= surfing_voxelselection(coord_white.data',coord_pial.data',topo.data', [radius 160],epiInfo, node_range ,[5,0,1]);
                LI=LI';
                save(fullfile(caret_subjDIR, ['surface_roi_160vox', '.mat']), 'LI','voxmin','voxmax','vORr');
                M=caret_struct('metric','data',vORr);
                caret_save(fullfile(caret_subjDIR, ['surface_roi_160vox.metric',]),M);
            end;
        end;
    case 'surf_MVA_do'%_______________________________________________
        sn=varargin{1};  
        hem=[1:2];
        atlas=2;
        glm=1; 
        vararginoptions({varargin{2:end}},{'atlas','hem','glm'}); 
        
        ldafunction= @df1_calcAcc;
        for s=sn
            for h=hem
                glmDirSubj=fullfile(baseDir,glmName{glm},subj_name{s});
                caretSubjDir = fullfile(caretDir,[atlasA{atlas},subj_name{s}],hemName{h});
                cd(caretSubjDir);                %name for metric outputfile
                
                metric_left=fullfile(caretSubjDir,sprintf('%s_accuracy_L.metric',subj_name{s}));
                metric_right=fullfile(caretSubjDir,sprintf('%s_accuracy_R.metric',subj_name{s}));
                %metric=fullfile(caretSubjDir,sprintf('%s_%s_accuracy_4reg.metric',subj_name{s}, seqCondition));
                metric_out = {[metric_left, ',accuracy'];[metric_right, ',accuracy']};
                
                surfaces=load(fullfile(caretSubjDir,['surface_roi_160vox','.mat']));
                load(fullfile(glmDirSubj,'SPM.mat'));
                T=load(fullfile(glmDirSubj,'SPM_info.mat')); 

                idx=find(T.regType==1); 
                for i=idx'
                    beta_images{i}=  fullfile(glmDirSubj,SPM.Vbeta(i).fname);
                end
                
                lmva_spm(surfaces,beta_images,metric_out,ldafunction,'params',{T.hand(idx)',T.seqNum(idx)',T.run(idx)'})
            end;
        end;
    case 'surf_MVA_do_decomp'%_______________________________________________
        sn=varargin{1};  
        hem=varargin{2};
        atlas=2; 
        
        ldafunction= @df1_decomp;
        for s=sn
            for h=hem
                glmDirSubj=fullfile(baseDir,glmName{1}, subj_name{s});
                %glmDirSubj=fullfile(glm4regDir, subj_name{s});
                if h== 1
                    caretSubjDir = fullfile(caretDir,[atlasA{atlas},subj_name{s}],'LeftHem');
                elseif h==2
                    caretSubjDir = fullfile(caretDir,[atlasA{atlas},subj_name{s}],'RightHem');
                end;
                %
                cd(caretSubjDir);                %name for metric outputfile
                
                metric_left=fullfile(caretSubjDir,sprintf('%s_decomp.metric',subj_name{s}));
                
                lda_metric = {[metric_left,  ',var_left'];...
                    [metric_left, ',var_right'];...
                    [metric_left, ',cov_condition'];...
                    [metric_left, ',var_sequences_left'];...
                    [metric_left, ',var_sequences_right'];...
                    [metric_left, ',covE'];...
                    [metric_left, ',covI'];...
                    [metric_left, ',noise']};
                
                surfaces=load(fullfile(caretSubjDir,['surface_roi_160vox','.mat']));
                load(fullfile(glmDirSubj,'SPM.mat'));
                
                for i=1:64
                    beta_images{i}=  fullfile(glmDirSubj,SPM.Vbeta(i).fname);
                end
                
                % surfaces=getrow(surfaces,[2000:2050]);
                
                %init decomp 
                % LOAD SPM_INFO.MAT AND BASE IT ON THS 
                numblocks=8;
                cc=repmat([1:4 1:4],1,numblocks);
                hh=repmat([1 1 1 1  2 2 2 2 ],1,numblocks);
                b=[];
                for i=1:numblocks
                    b=[b ones(1,8)*i];
                end;
                lmva_spm(surfaces,beta_images,lda_metric,ldafunction, 'params', {cc,hh,b})
            end;
        end;
    case 'surf_MVA_do_corr'%_______________________________________________
        sn=varargin{1};  
        hem=varargin{2};
        atlas=2; 
        
        ldafunction= @df1_calcCorr;
        for s=sn
            for h=hem
                glmDirSubj=fullfile(baseDir,glmName{1}, subj_name{s});
                if h== 1
                    caretSubjDir = fullfile(caretDir,[atlasA{atlas},subj_name{s}],'LeftHem');
                elseif h==2
                    caretSubjDir = fullfile(caretDir,[atlasA{atlas},subj_name{s}],'RightHem');
                end;
                cd(caretSubjDir);                %name for metric outputfile
                
                metric_left=fullfile(caretSubjDir,sprintf('%s_corr.metric',subj_name{s}));
                
                lda_metric = {...
                    [metric_left,  ',corrO'];...
                    [metric_left, ',corrE'];...
                    [metric_left, ',corrI'];...
                    [metric_left, ',corrED'];...
                    [metric_left, ',corrID'];...
                    [metric_left,  ',corrOm'];...
                    [metric_left, ',corrEm'];...
                    [metric_left, ',corrIm'];...
                    [metric_left, ',corrEDm'];...
                    [metric_left, ',corrIDm']};
                
                surfaces=load(fullfile(caretSubjDir,['surface_roi_160vox','.mat']));
                load(fullfile(glmDirSubj,'SPM.mat'));
                
                for i=1:64
                    beta_images{i}=  fullfile(glmDirSubj,SPM.Vbeta(i).fname);
                end
                
                % surfaces=getrow(surfaces,[2000:2050]);
                
                %init decomp 
                % LOAD SPM_INFO.MAT AND BASE IT ON THS 
                numblocks=8;
                cc=repmat([1:4 1:4],1,numblocks);
                hh=repmat([1 1 1 1  2 2 2 2 ],1,numblocks);
                b=[];
                for i=1:numblocks
                    b=[b ones(1,8)*i];
                end;
                lmva_spm(surfaces,beta_images,lda_metric,ldafunction, 'params', {cc,hh,b})
            end;
        end;
    case 'surf_MVA_dim'         % Runs dimensionality analysis on surface
        sn=varargin{1}; seqCondition=varargin{2};  hem=varargin{3};
        %df1_imana('surf_MVA_do',1, 'trained', 1);
        %df1_imana('surf_MVA_do',1, 'trained', 2);
        %df1_imana('surf_MVA_do',1, 'untrained', 1);
        %df1_imana('surf_MVA_do',1, 'untrained', 2);
        ldafunction= @sl1_calcDim;
        for s=sn
            for h=hem
                glmDirSubj=fullfile([glmDir,'_1'], subj_name{s});
                %glmDirSubj=fullfile(glm4regDir, subj_name{s});
                if h== 1
                    caret_subjDIR = fullfile(caretDir,['i',subj_name{s}],'LeftHem');
                elseif h==2
                    caret_subjDIR = fullfile(caretDir,['i',subj_name{s}],'RightHem');
                end;
                %
                cd(caret_subjDIR);                %name for metric outputfile
                metric=fullfile(caret_subjDIR,sprintf('%s_%s_accuracy_dim.metric',subj_name{s}, seqCondition));
                %metric=fullfile(caret_subjDIR,sprintf('%s_%s_accuracy_4reg.metric',subj_name{s}, seqCondition));
                lda_metric = {[metric, ',acc1'];[metric, ',acc2'];[metric, ',acc3'];[metric, ',dim'];[metric, ',numberVoxel']};
                
                surfaces=load(fullfile(caret_subjDIR,['surface_roi_160vox_',seqCondition,'.mat']));
                load(fullfile(glmDirSubj,seqCondition,'SPM.mat'));
                k=1;
                for i=SPM.xX.iC
                    beta_images{k}= fullfile(glmDirSubj,seqCondition,SPM.Vbeta(i).fname);
                    k=k+1;
                end
                if size(surfaces.LI,1)==1
                    surfaces.LI= surfaces.LI';
                end
                lmva_spm(surfaces,beta_images,lda_metric,ldafunction)%sl1_mva
            end;
        end;
    case 'surf_map_con'         % New mappint algorithm
        %df1_imana('surf_map_con', 1:11, 1:2)
        % map volume images to metric file and save them in individual surface folder
        sn=varargin{1};
        hemisphere=[1:2];
        atlas=2; 
        
        vararginoptions({varargin{2:end}},{'atlas','hemisphere'});
        glm=1; 
        fileList={'con_0003.img','con_0004.img', 'spmT_0003.img', 'spmT_0004.img','spmF_0001.img','spmF_0002.img'};
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
    case 'surf_map_seq'         % New mappint algorithm
        %df1_imana('surf_map_con', 1:11, 1:2)
        % map volume images to metric file and save them in individual surface folder
        sn=varargin{1};
        hemisphere=[1:2];
        atlas=2; 
        
        vararginoptions({varargin{2:end}},{'atlas','hemisphere'});
        glm=1; 
        fileList={'spmT_0005.img','spmT_0006.img','spmT_0007.img','spmT_0008.img',...
                  'spmT_0009.img','spmT_0010.img','spmT_0011.img','spmT_0012.img'};
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
                caret_save(fullfile(caretSDir,[subj_name{s} '_seq.metric']),M);
            end;
        end;
    case 'surf_avrgcoord'       % Makes an average FIDICUAL (WHITE) surface 
        for h=1:2
            surfaceGroupDir=[caretDir filesep 'fsaverage'  filesep hemName{h} ];
            cd(surfaceGroupDir);
            for i=1:length(subj_name);
                coord_name{i}=[caretDir filesep 'i' subj_name{i} filesep hemName{h} filesep hem{h} '.WHITE.coord'];
            end;
            topo_name= [surfaceGroupDir filesep hem{h} '.CLOSED.topo'];
            caret_avrgsurface('dir', caretDir,'topo_name', topo_name,'coord_name', coord_name,'fsaverage', h);
        end
    case 'surf_avrgsurfshape'   % generates 
        % Make group metric for surface shape into symmetric template (left
        % template) 
        INname={'surface_shape','surface_shape','surface_shape'};
        OUTname={'curv.surface_shape','sulc.surface_shape','area.surface_shape'}; 
        inputcol=[1 2 3]; 
        atlas=2;
        sn=[1:length(subj_name)];
        
        vararginoptions(varargin,{'atlas','sn'}); 
        for h=1:2
            surfaceGroupDir=[caretDir filesep atlasname{atlas}  filesep hemName{h} ];
            cd(surfaceGroupDir);
            for j=1:length(INname); %----loop over each input metric file and make a group metric file
                for i=sn; %----define names of subj metric files
                    infilenames{j}{i}=[caretDir filesep atlasA{atlas} subj_name{i} filesep hemName{h} filesep hem{h} '.' INname{j}];
                end;
                %----name for group metric file in average surface folder
                outfilenames{j}=[surfaceGroupDir filesep hem{h} '.' OUTname{j}];
                %----make the group metric
                fprintf('hem: %i  image: %i \n', h,j);
                caret_metricpermute(infilenames{j},'outfilenames',outfilenames(j),'inputcol',inputcol(j));
            end;
        end;
    case 'surf_correctAcc'
        name={'accuracy_L','accuracy_R'}; 
        for i=1
            for h=2 
                for s=1:2
                    infilename=[caretDir filesep 'i' subj_name{i} filesep hemName{h} filesep subj_name{i} '_' name{s} '.metric'];
                    A=caret_load(infilename); 
                    A.data(:,1)=A.data(:,1)*10e19; 
                    fprintf('%s\n',infilename); 
                    caret_save(infilename,A);
                end;
            end;
        end;
    case 'surf_makeGroup' % Make group metric on func and accuracy data. Do contrast and zacc later on these files!
        %df1_imana('surf_makeGroup',OPTIONAL: 'atlas',x,...)
       atlas=2; 
       
        INname={'func','func',...
            'func','func',...
            'func','func',...
             'accuracy_L','accuracy_R'};
        OUTname={'func_left','func_right',...
            'ttest_left','ttest_right',...
            'Ftest_left','Ftest_right',...
            'acc_left','acc_right'};
        inputcol= [1 2 3 4 5 6 1 1];
       
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
                if j<=4 % Normal functional: replace NaNs;
                    caret_metricpermute(infilenames{j},'outfilenames',outfilenames(j),'inputcol',inputcol(j),'replaceNaNs',1);
                else
                    caret_metricpermute(infilenames{j},'outfilenames',outfilenames(j),'inputcol',inputcol(j));
                end;
            end;
        end;
    case 'surf_makeGroupDecomp' % Make group metric on correlation analysis 
        %df1_imana('surf_makeGroup',OPTIONAL: 'atlas',x,...)
       atlas=2; 
       
            
          % 'corr','corr','corr',...
           %  'corr','corr',...
            
        INname={'decomp','decomp','decomp',...            
            'decomp','decomp',...
             'decomp','decomp'};
        inputcol= [1 2 3 4 5 6 7];
        OUTname={'decomp_handL','decomp_handR','decomp_covHand','decomp_seqL','decomp_seqR',...
            'decomp_covE','decomp_covI'};
       
        vararginoptions(varargin,{'atlas'});
        
        for h=1:2
            surfaceGroupDir=[caretDir filesep atlasname{atlas} filesep hemName{h} ];
            cd(surfaceGroupDir);
            for j=1:length(INname); %----loop over each input metric file and make a group metric file
                infilenames{j}={}; 
                for i=[1:10] %----define names of subj metric files
                    infilenames{j}{end+1}=[caretDir filesep atlasA{atlas} subj_name{i} filesep hemName{h} filesep subj_name{i} '_' INname{j} '.metric'];
                end;
                %----name for group metric file in average surface folder
                outfilenames{j}=[surfaceGroupDir filesep hem{h} '.' OUTname{j} '.metric'];
                %----make the group metric
                fprintf('hem: %i  image: %i \n', h,j);
                caret_metricpermute(infilenames{j},'outfilenames',outfilenames(j),'inputcol',inputcol(j),'replaceNaNs',1);
                M=caret_load(outfilenames{j}); 
                m(:,j)=mean(M.data,2); 
            end;
            eps=0.01; 
            m(:,8)=m(:,3)./sqrt((m(:,1)+eps).*(m(:,2)+eps)); 
            m(:,9)=m(:,6)./sqrt((m(:,4)+eps).*(m(:,5)+eps)); 
            m(:,10)=m(:,7)./sqrt((m(:,4)+eps).*(m(:,5)+eps)); 
            OUTname{8}='decompCorrH';
            OUTname{9}='decompCorrE';
            OUTname{10}='decompCorrI';
            M=caret_struct('metric','data',m,'column_name',OUTname); 
            caret_save('summary_decomp.metric',M); 
            
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
    case 'surf_zfunc'% Make the zacc files (on group level) and remove the NaN
        %df1_imana('surf_Tvalue_2_Zvalue')
        atlas=2; 
        INname={'ttest_left','ttest_right'};
        OUTnameZ={'zfunc_left','zfunc_right'};
        
        vararginoptions(varargin,{'atlas'});  

        for h=1:2
            surfaceGroupDir=[caretDir filesep atlasname{atlas}  filesep hemName{h} ];
            cd(surfaceGroupDir);
            for i=1:length(INname)
                iname=[hem{h} '.' INname{i} '.metric'];
                onameZ=[hem{h} '.' OUTnameZ{i} '.metric'];
                fprintf('hem: %i  image: %i \n', h,i);
                %----calculate the p-value from the T value with the cumulative function tcdf
                %df=784
                x= caret_load(iname);
                x.data=tcdf(x.data,1176); 
                eps=0.00000001; 
                x.data(x.data>1-eps)=1-eps; 
                x.data(x.data<eps)=eps; 
                x.data=norminv(x.data); 
                x.data(isnan(x.data))= 0;
                caret_save(onameZ,x); 
            end;
        end;
    case 'surf_zFtest'% Make the zacc files (on group level) and remove the NaN
        %df1_imana('surf_zFtest')
        INname={'ftest_left','ftest_right'};
        OUTnameZ={'zFtest_left','zFtest_right'};
        atlas=2; 
        vararginoptions(varargin,{'atlas'}); 
        
        for h=1:2
            surfaceGroupDir=[caretDir filesep atlasname{atlas} filesep hemName{h} ];
            cd(surfaceGroupDir);
            for i=1:2
                iname=[hem{h} '.' INname{i} '.metric'];
                onameZ=[hem{h} '.' OUTnameZ{i} '.metric'];
                fprintf('hem: %i  image: %i \n', h,i);
                %----calculate the p-value from the F value with the cumulative function fcdf
                %df1= 3.0 df2=1176
                
                caret_metriccalc({iname},onameZ,'norminv(fcdf(i1,3, 1176))'); %1-fcdf(i1,3, 1176) for the real p-value
                %replace NaN with zeros
                x= caret_load(onameZ);
                idx=find(isinf(x.data) & x.data<0); 
                x.data(idx)= -8.3;
                idx=find(isinf(x.data) & x.data>0); 
                x.data(idx)= +8.3;
                idx=find(isnan(x.data)); 
                x.data(idx)= 0;
                caret_save(onameZ,x);
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
    case 'surf_RGBfile'               % Make the RGB-overlay file for sham / tDCS groups 
        %sl1_imana('surf_RGBfile')
        % Makes flipped RGB comparision 
        for h=1:2
            surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
            cd(surfaceGroupDir)
            C=caret_load('summary.sep.metric');
            RGBdata(:,1)=C.data(:,2); % Red: Anodal BOLD signal
            RGBdata(:,3)=C.data(:,1); % Blue: Sham BOLD signal
            RGBdata(RGBdata(:,1)<0.8,1)=0;
            RGBdata(RGBdata(:,3)<0.8,3)=0;
            RGBdata(:,4)=C.data(:,8); % Red: Anodal Z-Accuracy
            RGBdata(:,6)=C.data(:,7); % Blue: Sham Z-Accuracy
            RGBdata(RGBdata(:,4)<1.96,4)=0;
            RGBdata(RGBdata(:,6)<1.96,6)=0;
            
            sc1=[0.8 3;0 10;0.8 3];
            sc2=[1.96 3;0 10;1.96 3];
            C=caret_struct('RGBpaint','data',RGBdata,'scales',{sc1,sc2},'column_name',{'contrast','accuracy'});
            caret_save([surfaceGroupDir '/' hem{h} '.summary.RGB_paint'],C);
        end;
    case 'surf_make_ROIpaint'         % Generates ROI paint file 
        h=1; 
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
            Paint=caret_struct('paint','data',ROI,'paintnames',{'S1','M1','PM','SMA','SPL'},'column_name',{'ROI'});
            caret_save(['ROIs_BA.paint'],Paint);
    case 'surf_area_stat'             % Compute activated / classification area for SPL,M1,S1,SMA together 
        %df1_imana('surf_count_voxel')
        atlas=2; 
        thres=[1 1 1.96 1.96];
        infilename={'%s.func_left.metric',...
                    '%s.func_right.metric',...
                    '%s.zacc_left.metric',...
                    '%s.zacc_right.metric'};
        D=[];
        vararginoptions(varargin,{'thres','infilename','atlas'}); 
        
        
        for h=1:2
            surfaceGroupDir=[caretDir filesep atlasname{atlas}  filesep hemName{h} ];
            cd(surfaceGroupDir);
            ROI=caret_load('ROIs_BA.paint'); 
            Area=caret_load([hem{h} '.area.surface_shape']); 
            Mask=ROI.data>0; 
            for i=1:length(infilename)
                M= caret_load(sprintf(infilename{i},hem{h}));
                numSubj=size(M.data,2);
                above=M.data>thres(i); 
                C.numVert= sum(above,1)';  
                C.area=sum(Area.data.*above,1)'; 
                C.sn= [1:numSubj]';
                C.tDCS_group=tDCS_group';
                % C.hand_group=hand_group';
                C.hand=ones(numSubj,1)*mod(i-1,2)+1;
                C.type=ones(numSubj,1)*floor((i-1)/2)+1; 
                C.mean=nanmean(M.data(Mask,:))';
                C.hem= ones(numSubj,1)*h;
                D= addstruct(D,C);
            end
        end;
        
        % Add behavioral stats from within scan analysis 
        % Add additional information on performance in pre and post-test 
        BS=dload(fullfile(behaviourDir,'df1_scan.dat')); 
        BP=dload(fullfile(behaviourDir,'SequenceTrained_MT','df1_preposttest.dat')); 
        BT=dload(fullfile(behaviourDir,'SequenceTrained_MT','df1_training.dat')); 
        for i=1:length(D.sn)
            indx=find(BS.sn==D.sn(i) & BS.hand==D.hand(i));
            D.ScanError(i,1)=BS.Error(indx);
            D.ScanMT(i,1)=BS.MT(indx);
            
            indx=find(BP.SN==D.sn(i) & BP.hand==D.hand(i) & BP.day==1) ;
            D.PreMT(i,1)=mean(BP.MT(indx));
       
            indx=find(BP.SN==D.sn(i) & BP.hand==D.hand(i) & BP.day==6 & BP.seqCondition==1) ;
            D.PostMT(i,1)=mean(BP.MT(indx));
            
            indx=find(BT.SN==D.sn(i) & BT.hand==D.hand(i)) ;
            if (isempty(indx))
                D.TrainMT(i,1)=NaN; 
            else
                D.TrainMT(i,1)=mean(BT.MT(indx)); 
            end; 
        end; 
        
        
        
        varargout={D}; 
    case 'surf_area_stat_plot'        % Plot the results of the simple area statistic
        D=df1_imana('surf_area_stat');
        subplot(2,1,1); 
        barplot([D.tDCS_group D.hand],D.area,'subset',D.type==1 & D.hand==1 & D.tDCS_group<3,'split',D.hem)
        subplot(2,1,2); 
        barplot([D.tDCS_group D.hand],D.area,'subset',D.type==2 & D.hand==1 & D.tDCS_group<3,'split',D.hem)
        [X,g]=pivottable([D.tDCS_group D.sn],D.hand,D.mean,'mean','subset',D.type==2); 
        ttest(X(g(:,1)==1,1),X(g(:,1)==2,1),2,'independent'); 
        varargout={D};      
    case 'surf_plotPattern'
        sn=varargin{1};
        h=varargin{2};
        groupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
        cd(groupDir);
        border=fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.CS_flat.border']);
        switch(h)
            case 2
                coord='rh.FLAT.coord';
                topo='rh.CLOSED.topo';
                data='rh.surface_shape';
                xlims=[4 33];
                ylims=[14 45];
            case 1
                coord='lh.FLATr.coord';
                topo='lh.CUT.topo';
                data='lh.surface_shape';
                xlims=[-20 10];
                ylims=[-10 20];
        end;
        
        B=caret_load(border);
        data=fullfile(caretDir,['x' subj_name{sn}],hemName{h},[subj_name{sn} '_seq.metric']);
        sshape=fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.surface_shape']);
        
        subplot(2,3,1);
        M=caret_plotflatmap('data',sshape,'border',B.Border,'topo',topo,'coord',coord,...
            'xlims',xlims,'ylims',ylims,'cscale',[-0.4 0.4]);

        colormap(gray);
        for i=1:4
            subplot(2,3,i+1);
            [M,d]=caret_plotflatmap('M',M,'col',i,'data',data,'cscale',[-6 12],...
                'border',B.Border,'topo',topo,'coord',coord);
            maxT(i)=max(d(:));
        end;
        
        
        mm=max(maxT);
        for i=1:4
            subplot(2,3,i+1);
            caxis([-2*mm/3 mm]);
        end;
        set(gcf,'PaperPosition',[1 1 8 5]);
        wysiwyg;
        % See also: tendigit1_imana('plot_fingerpatterns',5,1);
    case 'surf_meanRadius'
        %df1_imana('surf_meanRadius')
        D=[];
        for s=1:11
            for h=1:2
                if h== 1
                    caret_subjDIR = fullfile(caretDir,['i',subj_name{s}],'LeftHem');
                elseif h==2
                    caret_subjDIR = fullfile(caretDir,['i',subj_name{s}],'RightHem');
                end;
                %
                cd(caret_subjDIR);                %name for metric outputfile
                load(fullfile(caret_subjDIR,['surface_roi_160vox_trained.mat']));
                D(end+1, 1)= nanmean(vORr);
                load(fullfile(caret_subjDIR,['surface_roi_160vox_untrained.mat']));
                D(end, 2)= nanmean(vORr);
            end
        end
        D
        fprintf('mean radius for surface searchlight, %2.2f  \n', nanmean(D(:)));
    case 'surf_ROI_new'
        %df1_imana('surf_ROI_new')
        %df1_imana('ROI_surfaceROI2region', 1:11, 3)
        %R= df1_imana('ROI_data', 1:11, '_masked')
        %save('reg_data_masked.mat','-struct', 'R')
        % surface files
        T=[];
        Nodes.nam_lh={'???','1_IPC', '2_IPC', '3_S1', '4_S1',  '5_M1',  '6_vPM', '7_dPM', '8_SMA', '9_V1'};
        %Nodes.lh= [150852 72461 108510 46596 138083 145540 16296 109706 52293];
        Nodes.lh= [6838 151153 72168 105626 2564  75967 11741 23868 107669];
        %Nodes.nam_rh ={'???','1_IPC', '2_dPM', '3_M1', '4_SMA','5_vPM', '6_V1'};
        %Nodes.rh=[82436 128773 17873 50014 19972 112279];
        Nodes.nam_rh={'???','1_IPC', '2_IPC', '3_S1', '4_S1',  '5_M1',  '6_vPM', '7_dPM', '8_SMA', '9_V1'};
        Nodes.rh= [142190 59735 14739 19695 76678  152810 70857 133011 69960];
        
        hem
        
        for h=1:2
            surfaceGroupDir=[caretDir filesep 'fsaverage'  filesep hemName{h} ];
            cd(surfaceGroupDir);
            topo=caret_load([hem{h} '.CLOSED.topo']);
            coord=caret_load([hem{h} '.FIDUCIAL.coord']);
            P=zeros(size(coord.data,1),1);
            
            for r=1:size(Nodes.(char(hem(h))),2);
                coords=surfing_circleROI(coord.data,topo.data,Nodes.(char(hem(h)))(r),15); %8
                P(coords,1)=r;
            end;
            
            PF=caret_struct('paint','data',P, 'paintnames', Nodes.(['nam_',char(hem(h))]));
            caret_save([hem{h} '.circleROI.paint'],PF);
        end;
    case 'surf_ROI'
        % T= df1_imana('surf_ROI',1:7)
        nam_metric= {'.zacc_left','.zacc_right','.zTtest_left','.zTtest_right', '.ttest_untrained','.ttest_trained','.func_untrained','.func_trained',};
        %perctile=0; %---- do you like to look at 10:10:100 percent of the data
        threshold=1.5; %----data are thresholded on the mean of the ttest trained and untrained images
        paintROI= varargin{1};
        T=[];
        for h=1:2
            %----directories
            surfaceGroupDir=[caretDir filesep 'fsaverage'  filesep hemName{h} ];
            %----load the paint file
            C=caret_load([surfaceGroupDir filesep hem{h} '.ROI_26_05_2011_2.paint']);
            for i=1:numel(nam_metric)
                M_file{i}=caret_load([surfaceGroupDir filesep hem{h} nam_metric{i} '.metric']);
            end;
            d=[];
            %calculate mean of trained and untrained in the t-map
            d_mean=nanmean([M_file{5}.data(:,1:11) M_file{6}.data(:,1:11)],2);
            d_mean(d_mean< threshold)=NaN;
            for r=1:numel(paintROI)
                %define the nodes of the ROI
                nodes=C.data==paintROI(r);
                %----loop over the images
                for i=1:length(nam_metric)
                    for j=1:length(subj_name)
                        D=[];
                        %----subj data
                        d=M_file{i}.data(nodes,j);
                        %----apply threshold
                        d= d(~isnan(d_mean(nodes)));
                        %----look only at the x percent highest
                        %cat= split_data(d_mean,'numquant',100);
                        %d= d(cat>=perctile);
                        %----do struct
                        D.meanAcc= nanmean(d);
                        D.maxAcc= max(d);
                        D.minAcc= min(d);
                        D.numNodes= numel(d);
                        D.SN=j;
                        D.region= r;
                        D.image= i;
                        D.hem= h;
                        D.condition= mod(i,2)+1; %1 for trained 2 for untrained
                        T=addstruct(T,D);
                    end
                end
            end
        end
        %left hem
        L= tapply(T,{'SN' 'region' 'image'},  {T.meanAcc,'mean','name','acc'}, ...
            'subset', T.hem==1);
        %right hem
        R= tapply(T,{'SN' 'region' 'image'},  {T.meanAcc,'mean','name','acc'}, ...
            'subset', T.hem==2);
        %struct are in equal shape!
        %isequal([R.SN R.image R.region], [L.SN L.image L.region])
        %do plot lh - rh
        % L.difAcc= L.acc-R.acc
        % barplot([L.region   L.image ], L.difAcc, 'subset', (L.image==1|L.image==3) )
        varargout={T};
    case 'plot_hemAsymetrie'
        % df1_imana('plot_hemAsymetrie', T)
        %nam_metric= {'.zacc_untrained','.zacc_trained','.zTtest_untrained','.zTtest_trained', '.ttest_untrained','.ttest_trained','.func_untrained','.func_trained',};
        T= varargin{1};
        %left hem
        L= tapply(T,{'SN' 'region' 'image'},  {T.meanAcc,'mean','name','acc'}, ...
            'subset', T.hem==1);
        %right hem
        R= tapply(T,{'SN' 'region' 'image'},  {T.meanAcc,'mean','name','acc'}, ...
            'subset', T.hem==2);
        %struct are in equal shape!
        isequal([R.SN R.image R.region], [L.SN L.image L.region]);
        %do plot lh - rh
        L.difAcc= L.acc-R.acc;
        
        regionName={'inferior IPS'  'M1'  'PMd'  'PMv'  'S1' 'superior IPS' 'SMA'}
        for i=1:7
            %figure(1);
            subplot(4, 7, i);
            barplot([L.image ], L.difAcc, 'subset', (L.image==1|L.image==5) & L.region==i );
            title(['diff untr ',regionName{i}]);
        end
        for i=1:7
            %figure(3);
            subplot(4, 7, i+7);
            barplot([T.image T.hem], T.meanAcc, 'subset', (T.image==1|T.image==5) & T.region==i );
            title(['untr ',regionName{i}]);
        end
        for i=1:7
            %figure(2);
            subplot(4, 7, i+14);
            barplot([L.image ], L.difAcc, 'subset', (L.image==2|L.image==6) & L.region==i );
            title(['diff tr ',regionName{i}]);
        end
        for i=1:7
            %figure(4);
            subplot(4, 7, i+21);
            barplot([T.image T.hem], T.meanAcc, 'subset', (T.image==2|T.image==6) & T.region==i );
            title(['tr ',regionName{i}]);
        end
        for i=1:7
            [tt(i, 1) p(i, 1)]= ttest(L.acc(L.region==i& L.image==1), R.acc(R.region==i& R.image==1),  2, 'paired');
            [tt(i, 2) p(i, 2)]= ttest(L.acc(L.region==i& L.image==2), R.acc(R.region==i& R.image==2), 2, 'paired');
            [tt(i, 3) p(i, 3)]= ttest(L.acc(L.region==i& L.image==5), R.acc(R.region==i& R.image==5),  2, 'paired');
            [tt(i, 4) p(i, 4)]= ttest(L.acc(L.region==i& L.image==6), R.acc(R.region==i& R.image==6),  2, 'paired');
        end
        
        [tt(:, 1) tt(:, 3) tt(:, 2) tt(:, 4)]
    case 'plot_primVSsec_motorAreas'
        % df1_imana('plot_primVSsec_motorAreas', T)
        %nam_metric= {'.zacc_untrained','.zacc_trained','.zTtest_untrained','.zTtest_trained', '.ttest_untrained','.ttest_trained','.func_untrained','.func_trained',};
        T= varargin{1};
        %left hem
        L= tapply(T,{'SN' 'region' 'image'},  {T.meanAcc,'mean','name','acc'}, ...
            'subset', T.hem==1);
        %right hem
        R= tapply(T,{'SN' 'region' 'image'},  {T.meanAcc,'mean','name','acc'}, ...
            'subset', T.hem==2);
        %struct are in equal shape!
        isequal([R.SN R.image R.region], [L.SN L.image L.region]);
        
        %define primary and secondary motor areas
        x= zeros(size(L.region));
        %x(L.region==1|L.region==3|L.region==5|L.region==6|L.region==7)= 2;
        x(L.region==1|L.region==3|L.region==4|L.region==6|L.region==7)= 2;
        x(L.region==2|L.region==5)= 1;
        L.primSec= x;
        R.primSec= x;
        
        regionName={'inferior IPS'  'M1'  'PMd'  'PMv'  'S1' 'superior IPS' 'SMA'}
        subplot(2, 2, 1);
        barplot([L.image L.primSec ], L.acc, 'subset', (L.image==1|L.image==5)); %L.region
        title(['lh untr LEFT: z-acc RIGHT: t-bold']);
        subplot(2, 2, 2);
        barplot([L.image L.primSec], L.acc, 'subset', (L.image==2|L.image==6));  %L.region
        title(['lh tr LEFT: z-acc RIGHT: t-bold']);
        subplot(2, 2, 3);
        barplot([R.image R.primSec ], R.acc, 'subset', (R.image==1|R.image==5)); %L.region
        title(['rh untr LEFT: z-acc RIGHT: t-bold']);
        subplot(2, 2, 4);
        barplot([R.image R.primSec ], R.acc, 'subset', (R.image==2|R.image==6)); %L.region
        title(['rh tr LEFT: z-acc RIGHT: t-bold']);
        
        %----calc t-values and compare
        RR= tapply(R,{'SN' 'primSec' 'image'},  {R.acc,'mean','name','acc'});
        LL= tapply(L,{'SN' 'primSec' 'image'},  {L.acc,'mean','name','acc'});
        
        [tt( 1) p( 1)]= ttest(LL.acc(LL.primSec==1& LL.image==1), LL.acc(LL.primSec==2& LL.image==1),  2, 'paired');
        [tt( 2) p( 2)]= ttest(LL.acc(LL.primSec==1& LL.image==2), LL.acc(LL.primSec==2& LL.image==2),  2, 'paired');
        
        [tt( 3) p( 3)]= ttest(LL.acc(LL.primSec==1& LL.image==5), LL.acc(LL.primSec==2& LL.image==5),  2, 'paired');
        [tt( 4) p( 4)]= ttest(LL.acc(LL.primSec==1& LL.image==6), LL.acc(LL.primSec==2& LL.image==6),  2, 'paired');
        
        LL=RR;
        [tt( 5) p( 5)]= ttest(LL.acc(LL.primSec==1& LL.image==1), LL.acc(LL.primSec==2& LL.image==1),  2, 'paired');
        [tt( 6) p( 6)]= ttest(LL.acc(LL.primSec==1& LL.image==2), LL.acc(LL.primSec==2& LL.image==2),  2, 'paired');
        
        [tt( 7) p( 7)]= ttest(LL.acc(LL.primSec==1& LL.image==5), LL.acc(LL.primSec==2& LL.image==5),  2, 'paired');
        [tt( 8) p( 8)]= ttest(LL.acc(LL.primSec==1& LL.image==6), LL.acc(LL.primSec==2& LL.image==6),  2, 'paired');
        fprintf('     untrained \t\t trained \n');
        fprintf('     acc \t bold \t acc \t bold \n')
        fprintf('left:  %2.3f  %2.3f\t%2.3f  %2.3f \n', [tt(1) tt(3) tt(2) tt(4)])
        fprintf('right: %2.3f  %2.3f\t%2.3f  %2.3f \n', [tt(5) tt(7) tt(6) tt(8)])
    case 'surf_spatial_cog'
        %df1_imana('surf_spatial_cog')
        h=1;
        reg=varargin{1};%----which region are you looking at
        subj= 1:11;%[1:4 6:11];%
        switch reg
            case 1 %preMotor
                lims{1}=[-70 -30 30 70 ];
                %lims{1}=[-100 0 0 100 ];
                lims{2}=[0 40 45 85];
                %lims{2}=[-100 0 0 100 ];
                paintROI=1; %----Which ROI should be used
                flat{1}= '.FLAT.coord';
                flat{2}= '.CUT.topo';
                metric_file{1}= '.zacc_trained_stats.metric';
                metric_file{2}= '.zacc_untrained_stats.metric';
                group_th=2.5; %---- threshold for the plotted group average activation
                subj_th=0.8 %0.4; % ---- additional threshold for the individual ROI region in which the COG is calculated
            case 2 %SMA
                lims{1}=[0 50 35 70 ];
                lims{2}=[-50 -10 25 65 ];
                %lims{1}=[-100 0 0 100 ];
                %lims{2}=[-100 0 0 100 ];
                paintROI=2;%----Which ROI should be used
                flat{1}= '.MEDIAL.coord';
                flat{2}= '.MEDIAL.topo';
                %flat{1}= '.FLAT.coord';
                %flat{2}= '.CUT.topo';
                metric_file{1}= '.zacc_trained_stats.metric';
                metric_file{2}= '.zacc_untrained_stats.metric';
                group_th=2%2;1; 1.5
                subj_th=0.8%0.9%0.8%0.7; 1
            case 3 %SMA diff
                lims{1}=[0 50 35 70 ];
                lims{2}=[-50 -10 25 65 ];
                paintROI=2;%----Which ROI should be used
                flat{1}= '.MEDIAL.coord';
                flat{2}= '.MEDIAL.topo';
                metric_file{1}= '.zacc_diff_stats.metric';
                metric_file{2}= '.zacc_diff_stats.metric';
                group_th=1.1;
                subj_th=1;
        end
        AllT=[];
        drawpoints=1;
        drawborder=0;
        
        set(gcf,'PaperPosition',[2 2 8 4]);
        wysiwyg;
        for h=1:2
            T=[];
            groupDir=[caretDir filesep 'fsaverage'  filesep hemName{h} ];
            cd(groupDir);
            %----load the surface information
            coord=([hem{h} flat{1}]);
            topo=([hem{h} flat{2}]);
            FLAT=caret_load(coord);
            %----load the functional data
            CT=caret_load([hem{h} metric_file{1}]);
            CU=caret_load([hem{h} metric_file{2}]);
            %----calculate functional average data for RGB map [red=trained green='' blue=untrained]
            DATA=[mean(CT.data(:,subj), 2) zeros(size(CU.data,1),1) mean(CU.data(:,subj), 2) ];
            %---set small z-values to zeros => they will not be presented in the RGB map
            DATA(DATA<group_th)=0;
            %----
            DSCALE=[0.97 5.7;0 0;0.97 1.5];%DSCALE=[0.97 5.7;0 0;0.97 5.7];DSCALE=[0.97 5.7;0 0;0.97 1.5];DSCALE=[2.45 5.7;0 0;0.97 1.5]; % 50% and 40%
            %----load the shape file for visualisation
            S=caret_load([hem{h} '.surface_shape']);
            %----load the ROI file in which SMA and dPM is defined
            P=caret_load([hem{h} '.ROI_cog.paint']);
            %----in case you like to draw boarders as well
            drawborder=1;
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
    case 'surf_stat'
        % df1_imana('surf_stat',1, 0.005, 0.08, 1,'cSPM_zacc_diff.mat',1)
        
        %df1_imana('surf_stat',1, 0.001, 0.001, 2,'cSPM_zacc_diff.mat',-1) cSPM_func_diff
        % df1_imana('surf_stat',1, 0.07, 0.05, 1,'cSPM_zacc_diff.mat',1)
        % df1_imana('surf_stat',1, 0.01, 0.05, 1,'cSPM_zacc_diff.mat',1)
        % df1_imana('surf_stat',1, 0.01, 0.05, 1,'/5interSmooth/cSPM_zacc_diff.mat',1)
        % df1_imana('surf_stat',1, 0.01, 0.05, 1,'10iterSmooth/cSPM_func_avrg.mat',1)
        %----
        %df1_imana('surf_stat',hem, uncorrectedP, clusterTh, ,maskTh,'cSPM_zacc_diff.mat',sign)
        %----
        %u: height-threshold (if u<1 it is interpreted as uncorrected p)
        %k: size-threshold(mm2) if k<1 it is interpreted as corrected p cluster)
        h= varargin{1}; u= varargin{2}; k= varargin{3}; maskTh= varargin{4}; cSPMname=varargin{5}; sign=varargin{6};
        maskImage=varargin{7};
        surfaceGroupDir=[caretDir filesep 'fsaverage'  filesep hemName{h} filesep '']; %filesep '20iterSmooth'
        cd(surfaceGroupDir);
        load ([ hem{h},'.avrgsurface.mat']);
        load(maskImage);
        %load ('cSPM_zacc_avrg.mat');
        %load ('cSPM_func_avrg.mat');
        %load ('cSPM_psc_avrg.mat');
        mask= cSPM.b>maskTh;
        load (cSPMname);
        [TabDat,threshold_con]=caret_list(S,cSPM,u,k, 'mask',mask,'sign', sign );
        if ~isempty(TabDat.dat)
            varargout={TabDat.dat(:,end)-1};
        else
            varargout={[]};
        end
    case 'surf_results'
        % Instructions on how to use the caret_toolbox for multiple comp.
        % correction
        h=1;
        conname='zacc_diff';
        maskname='zacc_avrg';
        u=3;
        k=50;
        u_mask=1.64; % p<0.05 uncorrected
        
        % FIRST MAKE SURFACE (DO THIS BEFORE)
        surfaceGroupDir=[caretDir filesep 'fsaverage'  filesep hemName{h} ];
        cd(surfaceGroupDir)
        load([hem{h} '.surface.mat']);
        % Fist load the masking image and get the mask
        load(['cSPM_' maskname '.mat']);
        mask=cSPM.con(1).Z>u_mask;
        % Now get the contrast and make the list
        load(['cSPM_' conname '.mat']);
        T=caret_list(S,cSPM,u,k,'sort_by','p_cluster','mask',mask,...
            'save_thresholded',[hem{h} '.summary.metric,17']);
    case 'surf_AccBold_plot'
        name={'.func_untrained','.func_trained','.zacc_untrained','.zacc_trained'};
        node=varargin{1};
        h=varargin{2};
        T=[];
        groupDir=[caretDir filesep 'fsaverage'  filesep hemName{h} ];
        cd(groupDir);
        % Left vPM, SMA, Post Parietal
        % T=df1_imana('surf_AccBold_plot',[6398 133890 157327],1)
        %
        
        for i=1:4
            M{i}=caret_load([hem{h} name{i} '.metric']);
        end;
        Nr=length(node);
        for r=1:Nr
            for i=1:4
                vec=ones(7,1);
                D.data=M{i}.data(node(r)+1,1:7)';
                D.SN=[1:7]';
                D.region=vec*r;
                D.seqCond=vec*mod(i-1,2)+1;
                D.measure=vec*floor((i-1)/2)+1;
                T=addstruct(T,D);
            end;
            subplot(Nr,2,1+(r-1)*2);
            barplot([],T.data,'split',T.seqCond,'subset',T.measure==1 & T.region==r);
            subplot(Nr,2,2+(r-1)*2);
            barplot([],T.data,'split',T.seqCond,'subset',T.measure==2 & T.region==r);
            set(gca,'YAxisLocation','right');
        end;
        set(gcf,'PaperPosition',[2 2 3 1.5*Nr]);
        wysiwyg;
        varargout={T};
    case 'surf_160ROI'
        %----get the ROI around the nodes that are defined in the average surfaces
        %----REMEMBER separate will have different voxel selected! only for
        %the resliced coordinates of voxels are equal between conditions
        %D= df1_imana('surf_160ROI', 1:11, 1, 'separate')
        %D= df1_imana('surf_160ROI', 1:11, 4, 'separate')
        %D= df1_imana('surf_160ROI', 1:11, 1, 'resliced')
        %D= df1_imana('surf_160ROI', 1:11, 4, 'resliced')
        sn=varargin{1}; glmType=varargin{2}; roiType=varargin{3};
        D=[];
        Node=[53771  34732; 67720  66799; 81360  74175; 121510   70611; 149400  152793]; %defined on fsaverage
        name={'M1','SMA','aIPS', 'dPM', 'vPM'};
        for s=sn %----loop over subj
            for c=1:2 %----loop over conditions
                if strcmp(roiType, 'separate')
                    %----define the images
                    for i=1:32*glmType%----get the beta images
                        b_images{i}= sprintf('%s%04.0f%s', fullfile([glmDir,'_', num2str(glmType)], subj_name{s}, condition{c}, 'beta_'),i, '.img');
                    end
                    %lda_image= fullfile([glmDir,'_', num2str(glmType)], subj_name{s}, condition{c},'_lda_radius12.img');
                    mv_image= fullfile([glmDir,'_', num2str(glmType)], subj_name{s}, condition{c},'con_0006.img');
                    tmv_image= fullfile([glmDir,'_', num2str(glmType)], subj_name{s}, condition{c},'spmT_0006.img');
                    fseq_image=fullfile([glmDir,'_', num2str(glmType)], subj_name{s}, condition{c},'spmF_0001.img');
                    refV= spm_vol(fullfile([glmDir,'_', num2str(glmType)], subj_name{s}, condition{c},'mask_noskull.nii'));
                elseif strcmp(roiType, 'resliced')
                    %----define the images
                    for i=1:32*glmType%----get the beta images
                        b_images{i}= fullfile([resliceDir,'_', num2str(glmType)],subj_name{s}, [condition{c}(1),'_beta_',num2str(i),'.nii']);
                    end
                    %lda_image= fullfile([resliceDir,'_', num2str(glmType)],subj_name{s}, [condition{c}(1),'_lda_radius12.nii']);
                    mv_image= fullfile([resliceDir,'_', num2str(glmType)],subj_name{s}, [condition{c}(1),'_con_0006.nii']);
                    tmv_image= fullfile([resliceDir,'_', num2str(glmType)],subj_name{s}, [condition{c}(1),'_spmT_0006.nii']);
                    fseq_image=fullfile([resliceDir,'_', num2str(glmType)],subj_name{s}, [condition{c}(1),'_spmF_0001.nii']);
                    refV= spm_vol(fullfile([resliceDir,'_', num2str(glmType)],subj_name{s}, 'mask_noskull.nii'));
                end
                beta_images= char(b_images);
                %----load the images
                V=spm_vol(char(beta_images, mv_image, tmv_image, fseq_image));
                %---------------------
                for h=1:length(hem) %----loop over hemispheres
                    if h== 1
                        caret_subjDIR = fullfile(caretDir,['i',subj_name{s}],'LeftHem');
                        coord_pial= caret_load(fullfile(caret_subjDIR, 'lh.PIAL.coord'));
                        coord_white= caret_load(fullfile(caret_subjDIR, 'lh.WHITE.coord'));
                        topo= caret_load(fullfile(caret_subjDIR, 'lh.CLOSED.topo'));
                    elseif h==2
                        caret_subjDIR = fullfile(caretDir,['i',subj_name{s}],'RightHem');
                        coord_pial= caret_load(fullfile(caret_subjDIR, 'rh.PIAL.coord'));
                        coord_white= caret_load(fullfile(caret_subjDIR, 'rh.WHITE.coord'));
                        topo= caret_load(fullfile(caret_subjDIR, 'rh.CLOSED.topo'));
                    end
                    %----load the surface_roi file
                    % load(fullfile(caret_subjDIR,['surface_roi_160vox_',condition{1},'.mat'])); %----this needs to be one condition because otherwise the ROI are unequale in size!!!
                    
                    epiInfo.mat= refV.mat;
                    epiInfo.dim= refV.dim;
                    epiInfo.mask=spm_read_vols(refV);
                    node_range=1:size(coord_white.data,1);
                    [LI,voxmin,voxmax,vORr]= surfing_voxelselection(coord_white.data',coord_pial.data',topo.data', [radius 160],epiInfo, Node(:,h) ,[5,0,1]);
                    
                    %---------------------
                    for i=1:size(Node,1) %----loop over nodes
                        [I,J,K]=ind2sub(V(1).dim,LI{i}); %----get voxel coords
                        data=[];
                        for j=1:length(V) %----loop over images and get data
                            data(j,:)=spm_sample_vol(V(j),double(I),double(J),double(K),0)';
                        end;
                        %-----------------------
                        %----make data struct---
                        %-----------------------
                        T=[];
                        if glmType==1
                            %----masking only include voxels an accuracy and beta value and are not zero!
                            %----NOT checked on acc!
                            todelete= isnan(data(33,:))| isnan(data(32,:)) | data(32,:)==0;
                            if sum(todelete)>0
                                %keyboard
                            end
                            %----build up beta struct
                            T.isNaN= todelete';
                            %----do not delete the NaNs yet, because they will make
                            %----the region unqual between the conditions
                            [x y z]= spmj_affine_transform(I,J,K, inv(V(1).mat));
                            T.coord= [x; y; z]';
                            T.beta= data(1:32,:)';
                            T.mean= [mean(T.beta(:,1:4:32), 2) mean(T.beta(:,[1:4:32]+1), 2) mean(T.beta(:,[1:4:32]+2), 2) mean(T.beta(:,[1:4:32]+3), 2)];
                            T.con_0006= data(33,:)';
                            T.tcon_0006= data(34,:)';
                            T.fseq= data(35,:)';
                            T.sn= ones(size(T.fseq))*s;
                            T.region= ones(size(T.fseq))*i;
                            T.hemisphere= ones(size(T.fseq))*h;
                            T.regionName= repmat(name(i), size(T.fseq, 1),1);
                            T.seqCondition= ones(size(T.fseq))*c;
                        elseif glmType==4
                            %----masking only include voxels an accuracy and beta value and are not zero!
                            %----NOT checked on acc!
                            todelete= isnan(data(129,:))| isnan(data(128,:))| data(128,:)==0;
                            if sum(todelete)>0
                                %keyboard
                            end
                            %----build up beta struct
                            T.isNaN= todelete';
                            %----do not delete the NaNs yet, because they will make
                            %----the region unqual between the conditions
                            [x y z]= spmj_affine_transform(I,J,K,inv(V(1).mat));
                            T.coord= [x; y; z]';
                            T.beta= data(1:128,:)';
                            T.mean= [mean(T.beta(:,1:4:128), 2) mean(T.beta(:,[1:4:128]+1), 2) mean(T.beta(:,[1:4:128]+2), 2) mean(T.beta(:,[1:4:128]+3), 2)];
                            T.con_0006= data(129,:)';
                            T.tcon_0006= data(130,:)';
                            T.fseq= data(131,:)';
                            T.sn= ones(size(T.fseq))*s;
                            T.region= ones(size(T.fseq))*i;
                            T.hemisphere= ones(size(T.fseq))*h;
                            T.regionName= repmat(name(i), size(T.fseq, 1),1);
                            T.seqCondition= ones(size(T.fseq))*c;
                        end
                        D= addstruct(D,T);
                    end %node
                end; %condition
            end %hemisphere
        end; %subj
        % at the moment the regions in the D structure have the same number
        % of voxels across conditions
        % pivottable([D.region D.sn], D.seqCondition, D.region, 'length' )
        % but they also have NaNs that are independent for conditions.
        %----now delete all voxels that have a NaNs for a condition
        %         U= getrow(D, D.seqCondition==1);
        %         T= getrow(D, D.seqCondition==2);
        %         toDelete= (U.isNaN | T.isNaN);
        %         D= addstruct(getrow(U, ~toDelete), getrow(T, ~toDelete));
        %         D= rmfield(D, 'isNaN');
        %----save ROI struct
        cd(regDir);
        if strcmp(roiType, 'separate')
            if glmType==1
                save reg_data_separate_run_160.mat -struct D
            elseif glmType==4
                save reg_data_separate_block_160.mat -struct D
            end
        elseif strcmp(roiType, 'resliced')
            if glmType==1
                save reg_data_resliced_run_160.mat -struct D
            elseif glmType==4
                save reg_data_resliced_block_160.mat -struct D
            end
        end
        varargout={D};
    case 'ROI_surfaceROI2region'%______reads out paint file____________________________
        %df1_imana('ROI_surfaceROI2region', 1:11, 3)
        %df1_imana('ROI_surfaceROI2region', 1:11, 2)
        %df1_imana('ROI_surfaceROI2region', 1:11, 1)
        sn= varargin{1};
        paintType= varargin{2};
        for s=sn
            R=[];
            for h=1:2
                if paintType==1
                    C=caret_load(fullfile(caretDir,'fsaverage',hemName{h},[hem{h} '.circleROI.paint']));
                elseif paintType==2
                    C=caret_load(fullfile(caretDir,'fsaverage',hemName{h},[hem{h} '.ROI_26_05_2011_2.paint']));
                elseif paintType==3
                    C=caret_load(fullfile(caretDir,'fsaverage',hemName{h},[hem{h} '.circleROI.paint']));
                    %'???'    '1_IPC'    '2_IPC'    '3_S1'    '4_S1'    '5_M1'    '6_vPM'    '7_dPM'    '8_SMA'    '9_V1'
                    C_mask=caret_load(fullfile(caretDir,'fsaverage',hemName{h},[hem{h} '.ROI_26_05_2011_2.paint']));
                    %???'    'inferior_IPS'    'M1'    'PMd'    'PMv'    'S1'    'superior_IPS'    'SMA'
                    %which anatomical area should be used for the threshold
                    mask_idx= [6 6 5 5 2 4 3 7 0];
                end
                
                caretSubjDir=fullfile(caretDir,['i' subj_name{s}]);
                file=fullfile([glmDir,'_1'],subj_name{s}, 'mask.img');
                regionCounter= 1;
                for i=1:numel(C.paintnames)
                    if ~strcmp(C.paintnames{i}, '???')
                        fprintf('Working on subj: %i hemisphere: %i region: %s \n', s, h, C.paintnames{i})
                        nR=[];
                        nR.type='surf_nodes';
                        if paintType==3
                            nR.location=find(C.data==i-1 &C_mask.data==mask_idx(regionCounter));
                        else
                            nR.location=find(C.data==i-1);
                        end
                        nR.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                        nR.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                        nR.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                        nR.linedef=[5,0,1];
                        nR.image=file;
                        nR.name=[C.paintnames{i}];
                        nR.hemisphere=h;
                        nR.region= regionCounter;
                        regionCounter= regionCounter+1;
                        R{end+1}= nR;
                    end
                    
                end;
            end;
            R=region_calcregions(R);
            cd(regDir);
            if paintType==1
                save([subj_name{s} '_regions.mat'],'R');
            elseif paintType==2
                save([subj_name{s} '_regions_cyto.mat'],'R');
            elseif paintType==3
                save([subj_name{s} '_regions_masked.mat'],'R');
            end
        end;
        varargout={R};
    case 'ROI_data'% This gives an ROI structure
        % R= df1_imana('ROI_data', 1:11)
        % R= df1_imana('ROI_data', 1:11, '_cyto')
        % R= df1_imana('ROI_data', 1:11, '_mask')
        
        sn=varargin{1};
        if numel(varargin)>=2
            regionType=varargin{2};
        else
            regionType='';
        end
        
        D=[]; %----outputvariable
        for s=sn
            fprintf('subj: %i \n', s)
            %----define the images
            load(fullfile(regDir,[subj_name{s} '_regions',regionType,'.mat']));
            
            %----define the images
            for i=1:64%----get the beta images
                b_images{i}= sprintf('%s%04.0f%s', fullfile([glmDir,'_1'], subj_name{s}, 'beta_'),i, '.img');
            end
            %lda_image= fullfile([glmDir,'_', num2str(glmType)], subj_name{s}, condition{c},'_lda_radius12.img');
            mvLeft_image= fullfile([glmDir,'_1'], subj_name{s},'con_0004.img');
            mvRight_image= fullfile([glmDir,'_1'], subj_name{s},'con_0005.img');
            mvLeft_Timage= fullfile([glmDir,'_1'], subj_name{s},'spmT_0004.img');
            mvRight_Timage= fullfile([glmDir,'_1'], subj_name{s},'spmT_0005.img');
            beta_images= char(b_images);
            %----load the images
            V=spm_vol(char(beta_images, mvLeft_image, mvRight_image, mvLeft_Timage, mvRight_Timage));
            %----get region data
            data = region_getdata(V,R);
            for i=1:numel(R)
                T=[];
                if R{i}.region~=0
                    %----masking only include voxels an accuracy and beta value and are not zero!
                    %----NOT checked on acc!
                    todelete= isnan(data{i}(1,:))| isnan(data{i}(32,:)) | data{i}(32,:)==0;
                    if sum(todelete)>0
                        %keyboard
                    end
                    %----build up beta struct
                    T.isNaN= todelete';
                    %----do not delete the NaNs yet, because they will make
                    %----the region unqual between the conditions
                    T.coord= R{i}.data;
                    T.beta= data{i}(1:64,:)';
                    for j=0:7
                        T.mean(:,j+1)= mean(T.beta(:,[1:8:64]+j), 2);
                    end
                    T.con_0004= data{i}(65,:)';
                    T.con_0005= data{i}(66,:)';
                    T.tcon_0004= data{i}(67,:)';
                    T.tcon_0005= data{i}(68,:)';
                    T.sn= ones(size(T.tcon_0005))*s;
                    T.region= ones(size(T.tcon_0005))*R{i}.region;
                    T.hemisphere= ones(size(T.tcon_0005))*R{i}.hemisphere;
                    T.regionName= repmat({R{i}.name}, size(T.tcon_0005, 1),1);
                    
                    D= addstruct(D,T);
                end
            end %region
        end %subject
        
        %----save ROI struct
        cd(regDir);
        varargout={D};
    case 'ROI_corr'%
        % df1_imana('ROI_corr')
        
        R= load(fullfile(regDir,sprintf('reg_data_masked.mat')));
        
        D=[];
        for h=1:2
            figure;
            for r=1:numel(unique(R.region(R.hemisphere==h)))
                for s=1:max(R.sn)
                    
                    C=[];
                    %----region subset
                    rindx=find(R.sn==s & R.region==r  & R.hemisphere==h & ~R.isNaN & R.accSubspace_cat>8);%& R.accSubspace_cat==10);R.accSubspace_cat>8
                    X=getrow(R,rindx);
                    C.acc_left= mean(X.accSubspace_left);
                    C.acc_right= mean(X.accSubspace_right);
                    
                    %----normal corr
                    [C.covO C.varO C.covI C.varI C.covE C.varE] = df1_imana('ROI_calc_covVar', X.mean);
                    %----calc correlations
                    C.corrO= C.covO/sqrt(C.varO);
                    C.corrE= C.covE/sqrt(C.varE+0.00000001);
                    C.corrI= C.covI/sqrt(C.varI+0.00000001);
                    
                    
                    %----decomp
                    numblocks=8;
                    cc=repmat([1:4 1:4],1,numblocks);
                    hh=repmat([1 1 1 1  2 2 2 2 ],1,numblocks);
                    b=[];
                    for i=1:numblocks
                        b=[b ones(1,8)*i];
                    end;
                    %
                    %[C.var_left, C.var_right, C.cov_condition, C.var_sequences_left, C.var_sequences_right, C.covE, C.covI, C.noise]
                    dC= df1_decomp_struct(X.beta,cc,hh,b);
                    C= addstruct(C, dC);
                    
                    %----
                    C.sn= s;
                    C.hemisphere= h;
                    C.region= r;
                    C.regionName{1}= char(unique(X.regionName));
                    %----calc correlations
                    C.d_corrE= C.d_covE/sqrt(C.d_var_sequences_left * C.d_var_sequences_right +0.0001); %+0.00000001
                    C.d_corrI= C.d_covI/sqrt(C.d_var_sequences_left * C.d_var_sequences_right +0.0001);
                    
                    D= addstruct(D, C);
                    
                    fprintf('hem: %i \t subj: %i \t region: %s \n ',h, s,char(unique(X.regionName)));
                end
                S= getrow(D, D.hemisphere==h& D.region==r);
                %myboxplot([ones(size(S.corrO)); ones(size(S.corrO))*2], [S.corrI; S.corrE]-[S.corrO; S.corrO])
                subplot(6, numel(unique(R.region(R.hemisphere==h))), r)
                myboxplot([ones(size(S.d_corrE)); ones(size(S.d_corrE))*2], [S.d_corrI; S.d_corrE])
                title(char(unique(X.regionName)));
                %
                subplot(6, numel(unique(R.region(R.hemisphere==h))), r+numel(unique(R.region(R.hemisphere==h))))
                myboxplot([ones(size(S.corrO)); ones(size(S.corrO))*2], [S.corrI; S.corrE]-[S.corrO; S.corrO])
                title('normal corr');
                
                subplot(6, numel(unique(R.region(R.hemisphere==h))), r+numel(unique(R.region(R.hemisphere==h)))*2)
                myboxplot([ones(size(S.corrO)); ones(size(S.corrO))*2], [S.acc_left; S.acc_right])
                title('acc');
                
                subplot(6, numel(unique(R.region(R.hemisphere==h))), r+numel(unique(R.region(R.hemisphere==h)))*3)
                myboxplot([ones(size(S.corrE)); ones(size(S.corrE))*2],...
                    [S.d_var_sequences_left; S.d_var_sequences_right])
                title('  seq L seq R ');
                
                subplot(6, numel(unique(R.region(R.hemisphere==h))), r+numel(unique(R.region(R.hemisphere==h)))*4)
                myboxplot([ones(size(S.corrE)); ones(size(S.corrE))*2; ones(size(S.corrE))*3],...
                    [S.d_var_left; S.d_var_right; S.d_noise])
                title('condition L condition R  noise');
                
                subplot(6, numel(unique(R.region(R.hemisphere==h))), r+numel(unique(R.region(R.hemisphere==h)))*4)
                plot(S.d_var_sequences_left .* S.d_var_sequences_right, S.d_covE, '*');
                hold on;
                plot(S.d_var_sequences_left .* S.d_var_sequences_right, S.d_covI, 'r*');
                hold off
                fprintf('intrinsic \n')
                ttest(S.d_corrI, [], 2, 'onesample')
                fprintf('extrinsic \n')
                ttest(S.d_corrE, [], 2, 'onesample')
                waitforbuttonpress
                %               keyboard
                %                 D=[];
            end;
        end;
        figure
        %plot_reg= [1 5 6 4 7 8 9; 1 5 6 4 7 8 9]
        plot_reg= [1:9; 1:9]
        for h=1:2
            a=1;
            for r=plot_reg(h, :)
                S= getrow(D, D.hemisphere==h& D.region==r);
                subplot(size(plot_reg,1), size(plot_reg,2), a+9*(h-1)); myboxplot([ones(size(S.d_corrE)); ones(size(S.d_corrE))*2], [S.d_corrI; S.d_corrE]);
                ylim([-1, 1]);
                title(char(unique(S.regionName)))
                a= a+1;
                fprintf('hem: %i  \t region: %s \n ',h,char(unique(S.regionName)));
                fprintf('intrinsic \n')
                ttest(S.d_corrI, [], 2, 'onesample')
                fprintf('extrinsic \n')
                ttest(S.d_corrE, [], 2, 'onesample')
            end
        end
    case 'ROI_calc_covVar'
        X=varargin{1};
        X_cov= cov(X);
        %----baseline
        C.covO= mean([X_cov(5, 3) X_cov(5, 4) X_cov(6, 3) X_cov(6, 4) X_cov(7, 1) X_cov(7, 2) X_cov(8, 1) X_cov(8, 2)]);
        C.varO=mean([X_cov(3, 3)* X_cov(5, 5), ...
            X_cov(5,5)* X_cov(4,4),...
            X_cov(6,6)* X_cov(3,3),...
            X_cov(6,6)* X_cov(4,4),...
            X_cov(1,1)* X_cov(7,7),...
            X_cov(2,2)* X_cov(7,7),...
            X_cov(1,1)* X_cov(8,8),...
            X_cov(2,2)* X_cov(8,8)]);
        
        %----extrinsic This was coded wrongly before as extrinic (c is on
        %the level extrinsic cueing 
        C.covE= mean([X_cov(5, 1) X_cov(6, 2) X_cov(7, 3) X_cov(8, 4)]);
        C.varE= mean([X_cov(5, 5)* X_cov(1, 1),...
            X_cov(6, 6)* X_cov(2, 2),...
            X_cov(7, 7)* X_cov(3, 3),...
            X_cov(8, 8)* X_cov(4, 4)]);
        
        %----intrinsic: This was coded wrongly before as extrinic (c is on
        %the level extrinsic cueing 
        C.covI= mean([X_cov(6, 1) X_cov(5, 2) X_cov(7, 4) X_cov(8, 3)]);
        C.varI= mean([X_cov(6, 6)* X_cov(1, 1),...
            X_cov(5, 5)* X_cov(2, 2),...
            X_cov(7, 7)* X_cov(4, 4),...
            X_cov(8, 8)* X_cov(3, 3)]);
        
        varargout={C.covO C.varO C.covE C.varE C.covI C.varI};
    case 'ROI_addBasalGangliaROI'
        %df1_imana('ROI_addBasalGangliaROI', 1:11)
        sn=varargin{1};
        D=[]; %----outputvariable
        for s=sn
            load(fullfile(regDir,[subj_name{s} '_regions.mat']));
            %----careful this number has to change if the previous region number changes!
            rmax= 7;
            for h=1:2
                for i=1:numel(bgNames)
                    fprintf('Working on subj: %i region: %s \n', s, [bgNames{i},'_',hem{h}])
                    V= spm_vol(fullfile(baseDir, 'basal_ganglia', [bgNames{i},'_',hem{h}, '.nii']));
                    C= spm_read_vols(V);
                    %----deform basal ganglia ROI to individual space
                    nam_def = fullfile(anatomicalDir,subj_name{s}, [subj_name{s},'_anatomical_seg_sn.mat']);
                    [def,def_mat] = spmdefs_get_sn2def(char(nam_def));
                    T=load(nam_def);
                    VSource=T.VF;
                    [i_def,i_defMat] = spmdefs_get_inv(def,def_mat,VSource);
                    %----do deformation
                    cd(fullfile(baseDir, 'basal_ganglia', 'subj_ROI'));
                    spmdefs_apply_def(i_def,i_defMat,fullfile(baseDir, 'basal_ganglia', [bgNames{i},'_',hem{h}, '.nii']),0,[subj_name{s},'_',bgNames{i},'_',hem{h}, '.nii']);
                    %xyz_tmp(:,1:3)= spmdefs_transformM(i_def,i_defMat,C.data>0);
                    %----define region
                    nR.file= fullfile(baseDir, 'basal_ganglia', 'subj_ROI',[subj_name{s},'_',bgNames{i},'_',hem{h}, '.nii']);
                    nR.threshold=0;
                    nR.type='image';
                    nR.hemisphere=h;
                    nR.name= bgNames{i};
                    nR.region= i+rmax;
                    R{end+1}= nR;
                end
            end
            R=region_calcregions(R);
            cd(regDir);
            save([subj_name{s} '_regions.mat'],'R');
        end
    case 'ROI_add_variable' % adds a new variable to the reg_data_... struct and overwrites the old one
        % df1_imana('ROI_add_variable', 'fcn', 'spatialsubspace', 'varName', {'accSubspace_left', 'accSubspace_right', 'radiusSubspace'});
        % df1_imana('ROI_add_variable', 'fcn', 'acc2z', 'varName', {'zaccSubspace_left', 'zaccSubspace_right',});
        % df1_imana('ROI_add_variable', 'fcn', 'quantiles_acc', 'varName', {'accSubspace_left_cat', 'accSubspace_right_cat','accSubspace_cat'});
        % df1_imana('ROI_add_variable', 'fcn', 'quantiles_func', 'varName', {'func_cat'});
        
        % df1_imana('ROI_add_variable', 'fcn', 'spatialsubspace', 'varName', {'accSubspace_left', 'accSubspace_right', 'radiusSubspace'}, 'dataext', {''});
        % df1_imana('ROI_add_variable', 'fcn', 'acc2z', 'varName', {'zaccSubspace_left', 'zaccSubspace_right'}, 'dataext', {''});
        % df1_imana('ROI_add_variable', 'fcn', 'quantiles_acc', 'varName', {'accSubspace_left_cat', 'accSubspace_right_cat','accSubspace_cat'}, 'dataext', {''});
        
        %l4_imana('ROI_add_variable', 'fcn', 'spatialsubspace', 'varName', {'accSubspace_left', 'accSubspace_right', 'radiusSubspace'}, 'dataext', {'_masked'});
        %df1_imana('ROI_add_variable', 'fcn', 'acc2z', 'varName', {'zaccSubspace_left', 'zaccSubspace_right'}, 'dataext', {'_masked'});
        %df1_imana('ROI_add_variable', 'fcn', 'quantiles_acc', 'varName', {'accSubspace_left_cat', 'accSubspace_right_cat','accSubspace_cat'}, 'dataext', {'_masked'});
        param={};
        fcn='';
        varName= {'newVar'};
        dataext={'_cyto'};
        crossCondition= 0;
        whichROIs= 1:length(regionName);
        vararginoptions(varargin,{'fcn','param','dataext','varName', 'crossCondition', 'whichROIs'});
        %mkdir(fullfile(regDir,dataext{1}, dataext{2}));
        %cd(fullfile(regDir,dataext{1}));
        %load region data
        R= load(fullfile(regDir,sprintf('reg_data%s.mat',dataext{1})));
        %load behavioural data
        % B= load(fullfile(behaviourDir,sprintf('D_%s.mat',dataext{2})));
        
        D=[];
        for s= unique(R.sn)'
            %----index for behaviour datastruct
            %bindx(:,1)= B.sn==s& B.seqCondition==1;
            %bindx(:,2)= B.sn==s& B.seqCondition==2;
            for h=1:2 %unique(R.hemisphere)'
                for r=1:numel(unique(R.region(R.hemisphere==h)))
                    fprintf('subj: %i hemisphere: %i region: %i \n', s, h, r);
                    %----index for imaging datastruct
                    rindx(:,1)=R.sn==s & R.hemisphere==h & R.region==r;
                    %---- if you like to calculate values across conditions
                    %----get region data
                    sR= getrow(R,rindx(:,1));
                    %----get behavioural data
                    %sB= getrow(B,bindx(:,1)| bindx(:,2));
                    
                    %----do calculations
                    %X= feval(@df1_imana,fcn, sR, sB, param{:});
                    X= feval(@df1_imana,fcn, sR, param{:});
                    %----write the results in the data struct with the defined variable names
                    for i= 1:numel(varName)
                        sR.(varName{i})= X.(varName{i});
                    end
                    D=addstruct(D,sR);
                    %---- if you like to calculate values separate for each condition
                    
                end % end region
            end %heimsphere
        end % end subjects
        cd(regDir);
        save (fullfile(regDir,sprintf('reg_data%s.mat',dataext{1})), '-struct', 'D');
    case 'spatialsubspace'%_______PLUG IN functions for ROI_add_variable
        T=varargin{1};
        left_idx=  repmat([ones(1,4), zeros(1,4)], 1, 8)==1;
        right_idx= repmat([ones(1,4), zeros(1,4)], 1, 8)==0;
        if size(T.beta,1)<=160
            %----
            a=sl1_calcAcc(T.beta(:,left_idx));
            Xleft.acc= repmat(a', size(T.beta,1), 1);
            %----
            a=sl1_calcAcc(T.beta(:,right_idx));
            Xright.acc= repmat(a', size(T.beta,1), 1);
        else
            Xleft= lmva_spatialsubspace(T.beta(:,left_idx), @sl1_calcAcc,T.coord,160);
            Xright= lmva_spatialsubspace(T.beta(:,right_idx), @sl1_calcAcc,T.coord,160);
        end
        D.accSubspace_left=  Xleft.acc(:,1);
        D.accSubspace_right= Xright.acc(:,1);
        D.radiusSubspace= Xleft.acc(:,2);
        varargout={D};
    case 'acc2z'
        T=varargin{1};
        D.zaccSubspace_left= (T.accSubspace_left-25)/76.5466;
        D.zaccSubspace_right= (T.accSubspace_right-25)/76.5466;
        varargout={D};
    case 'quantiles_acc'
        T= varargin{1};
        D.accSubspace_left_cat=split_data(T.accSubspace_left,'numquant',10);
        D.accSubspace_right_cat=split_data(T.accSubspace_right,'numquant',10);
        D.accSubspace_cat=split_data(mean([T.accSubspace_left T.accSubspace_right], 2),'numquant',10);
        varargout={D};
    case 'quantiles_func'
        T= varargin{1};
        D.func_cat=split_data(mean([T.con_0004 T.con_0005], 2),'numquant',10);
        varargout={D};
    case 'mni_coords'
        T=varargin{1};
        %anatomicalDir s07_anatomical_seg_sn.mat
        %----deforms region coords in template space
        %INIT deformation
        nam_def = fullfile(anatomicalDir, subj_name{unique(T.sn)},[subj_name{unique(T.sn)},'_anatomical_seg_sn.mat']);
        [def,def_mat] = spmdefs_get_sn2def(nam_def);
        C=load(nam_def);
        VSource=C.VF;
        [i_def,i_defMat] = spmdefs_get_inv(def,def_mat,VSource);
        %----do deformation
        D.mni_coord(:,1:3)= spmdefs_transformM(i_def,i_defMat,T.coord);
        varargout={D};
    case 'ROI_stats_new'% caluculates some statistic on the region...
        %T= df1_imana('ROI_stats_new','fcn', 'stats', 'param',{2.75, 0})
        param={};
        dataext={'cyto'};
        quantiles=1;
        threshold=''
        thVariable= 'func_cat'; %'accSubspace_cat';
        vararginoptions(varargin,{'fcn','param','dataext','quantiles', 'threshold', 'thVariable'});
        
        R= load(fullfile(regDir,sprintf('reg_data_%s.mat',dataext{1})));
        if isempty(threshold)
            rindx= R.(thVariable)>10-quantiles;
        else
            rindx= R.(thVariable)>threshold;
        end
        R= getrow(R,rindx)
        S=[];
        for s=1:max(R.sn)
            for r=1:numel(regionName)
                for h=1:2
                    %----behavioural subset
                    % bindx=find(B.sn==s & B.seqCondition==c);
                    
                    C=[];
                    %----region subset
                    rindx=find(R.sn==s & R.region==r & R.hemisphere==h);
                    if (isempty(rindx) & isempty(bindx))
                        fprinf('check struct');
                        break;
                    end
                    %C=feval(@df1_imana,fcn,getrow(R,rindx),getrow(B, bindx),param{:});
                    C=feval(@df1_imana,fcn,getrow(R,rindx),param{:});
                    fn=fieldnames(C);
                    vec=ones(size(C.(fn{1}),1),1);
                    C.sn=s*vec;
                    C.region=r*vec;
                    C.hemisphere=h*vec;
                    S=addstruct(S,C);
                    clc;
                    fprintf('subj: %i \t region: %i \n ', s,r);
                end
            end;
        end;
        varargout={S};
    case 'stats'%_______PLUG IN functions for ROI_stats
        %These are the basic examples to plug in to the ROI_stats function
        T=varargin{1};
        func_threshold=varargin{2};
        zAcc_threshold=varargin{3};
        %a=sl1_calcAcc(T.beta, {'standard'});
        %D.acc=a(1);
        D.funcVoxel= numel(T.con_0004);
        D.numfuncLeft= sum(T.tcon_0004>func_threshold);
        D.numfuncRight= sum(T.tcon_0005>func_threshold);
        D.meanfuncLeft= mean(T.tcon_0004(T.tcon_0004>func_threshold));
        D.meanfuncRight= mean(T.tcon_0005(T.tcon_0005>func_threshold));
        D.funcLeft= mean(T.con_0004);
        D.funcRight= mean(T.con_0005);
        
        D.numzAccLeft= sum(T.zaccSubspace_left>zAcc_threshold);
        D.numzAccRight= sum(T.zaccSubspace_right>zAcc_threshold);
        D.meanzAccLeft= mean(T.zaccSubspace_left(T.zaccSubspace_left>zAcc_threshold));
        D.meanzAccRight= mean(T.zaccSubspace_right(T.zaccSubspace_right>zAcc_threshold));
        
        %D.meanMT=  mean(B.MT_mean);   %----Mean movement time during the scan
        %D.meanF= mean (mean([B.f0_mean B.f1_mean B.f2_mean B.f3_mean B.f4_mean]));
        %X= feval(@df1_imana,'variance_components_MT',T,B);
        varargout={D};
        
        %         T=varargin{1};
        %         B=varargin{2};
        %         %----standard acc
        %         a=sl1_calcAcc(T.beta, {'standard'});
        %         D.acc=a(1);
        %
        %         %----do regression with mt
        %         %----substract the mean response of mt
        %         mt= B.MT_median- mean(B.MT_median);
        %         %----make design matrix
        %         Zgood= [kron(ones(4*8,1),eye(4))];
        %         Zfact= mt;
        %         a=sl1_calcAcc(T.beta, 'regress', Zgood(:,2:end), Zfact);
        %         D.acc_regress_mt=a(1);
        %
        %         %----do regression with force
        %         %----substract the mean response of force
        %         force= B.f_mean- mean(B.f_mean);
        %         Zfact= force;
        %         a=sl1_calcAcc(T.beta, 'regress', Zgood(:,2:end), Zfact);
        %         D.acc_regress_f=a(1);
        %
        %         %----calculate z-values
        %         D.zacc= (D.acc-25)/3.8273;
        %         D.zacc_regress_mt= (D.acc_regress_mt-25)/3.8273;
        %         D.zacc_regress_f= (D.acc_regress_f-25)/3.8273;
        %
        varargout={D};
    case 'ROI_asymmetri'
        %T= df1_imana('ROI_stats_new','fcn', 'stats', 'param',{2.75, 0})
        %df1_imana('ROI_asymmetri',T, 'numfunc')
        %figure
        %df1_imana('ROI_asymmetri',T, 'numzAcc')
        %figure
        %df1_imana('ROI_asymmetri',T, 'meanfunc')
        %figure
        %df1_imana('ROI_asymmetri',T, 'meanzAcc')
        T=varargin{1};
        nam=varargin{2};
        S=[];
        for r=1:numel(regionName)
            for s=1:max(T.sn)
                C=[];
                rindx=find(T.sn==s & T.region==r & T.hemisphere==2);
                lindx=find(T.sn==s & T.region==r & T.hemisphere==1);
                %----left hem
                C.lh_left=  T.([nam,'Left'])(lindx);
                C.lh_right= T.([nam,'Right'])(lindx);
                %----right hem
                C.rh_left=  T.([nam,'Left'])(rindx) ;
                C.rh_right= T.([nam,'Right'])(rindx);
                % movement I
                C.I_left=  C.lh_left/ ([C.lh_left + C.rh_left]);
                C.I_right= C.rh_right/ ([C.rh_right + C.lh_right]);
                %----hemispheric asymmetry (HA)
                C.HA= C.I_left-C.I_right;
                %----number Voxel in ROI
                C.rh_funcVoxel= T.funcVoxel(rindx);
                C.lh_funcVoxel= T.funcVoxel(lindx);
                
                C.sn=s;
                C.region=r;
                S=addstruct(S,C);
            end
            %----ttest
            fprintf('--------------\n')
            fprintf('asymmerie %s\n',char(regionName(r)))
            ttest(S.I_left(S.region==r), S.I_right(S.region==r), 2, 'paired');
        end
        subplot (2, 1, 1); myboxplot( S.region, [S.I_left S.I_right]);
        %ylim([0, 1])
        subplot (2, 1, 2); myboxplot( S.region, [S.lh_funcVoxel S.rh_funcVoxel]);
        fprintf('\n')
        for r=1:numel(regionName)
            %----ttest
            fprintf('--------------\n')
            fprintf('overall number of voxels %s\n',char(regionName(r)))
            ttest(S.lh_funcVoxel(S.region==r), S.rh_funcVoxel(S.region==r), 2, 'paired');
        end
    case 'ROI_stats'% caluculates some statistic on the region...
        %T= df1_imana('ROI_stats','fcn', 'stats')
        param={};
        dataext={'cyto'};
        quantiles=1:10;
        vararginoptions(varargin,{'fcn','param','dataext','quantiles'});
        
        R= load(fullfile(regDir,sprintf('reg_data_%s.mat',dataext{1})));
        %B= load(fullfile(behaviourDir,sprintf('D_%s.mat',dataext{2})));
        
        S=[];
        for s=1:max(R.sn)
            for r=1:numel(regionName)
                for h=1:2
                    for c=1:2
                        %----behavioural subset
                        % bindx=find(B.sn==s & B.seqCondition==c);
                        if ~isempty(quantiles)
                            for q=quantiles % loop over quantiles
                                C=[];
                                %----region subset
                                rindx=find(R.sn==s & R.region==r & R.seqCondition==c & R.hemisphere==h & R.accSubspace_cat>10-q);
                                if (isempty(rindx) & isempty(bindx))
                                    fprinf('check struct');
                                    break;
                                end
                                C=feval(@df1_imana,fcn,getrow(R,rindx),getrow(B, bindx),param{:});
                                fn=fieldnames(C);
                                vec=ones(size(C.(fn{1}),1),1);
                                C.sn=s*vec;
                                C.region=r*vec;
                                C.hemisphere=h*vec;
                                C.seqCondition=c*vec;
                                C.quantile=q*vec;
                                S=addstruct(S,C);
                                clc;
                                fprintf('subj: %i \t region: %i \t condition: %i ', s,r,c);
                            end;
                        else
                            C=[];
                            %----region subset
                            rindx=find(R.sn==s & R.region==r & R.seqCondition==c & R.hemisphere==h);
                            if (isempty(rindx) & isempty(bindx))
                                fprinf('check struct');
                                break;
                            end
                            C=feval(@df1_imana,fcn,getrow(R,rindx),getrow(B, bindx),param{:});
                            fn=fieldnames(C);
                            vec=ones(size(C.(fn{1}),1),1);
                            C.sn=s*vec;
                            C.region=r*vec;
                            C.hemisphere=h*vec;
                            C.seqCondition=c*vec;
                            S=addstruct(S,C);
                            clc;
                            fprintf('subj: %i \t region: %i \t condition: %i ', s,r,c);
                            
                        end
                    end
                end
            end;
        end;
        varargout={S};
    case 'ROI_plot_stats'   % Make a barplot across different ROIs
        %df1_imana('ROI_plot_stats',T, 'accM', [6  9 10 11 12 13 14])
        T=varargin{1};      % ROI-stats structure
        dep=varargin{2};    % This is the dependent variable to plot
        reg= varargin{3};   % regions that should be ploted
        for r=reg
            subplot(1,numel(reg),find(reg==r));
            barplot(T.side,T.(dep),'split',[T.hemisphere T.seqCondition],'subset',T.regType==r);
            if (strcmp(dep,'accM') | strcmp(dep,'accSubspaceM')| strcmp(dep,'accSubspaceM_prct' ) | strcmp(dep,'accM_decom') ...
                    |strcmp(dep,'acc_decom') | strcmp(dep,'acc_decomNoSVD')| strcmp(dep,'acc') | strcmp(dep,'acc_regress'))
                set(gca,'YLim',[0 45]);
            elseif (strcmp(dep,'zaccSubspace_prct'))
                set(gca,'YLim',[0 0.58]);
            else
                set(gca,'YLim',[-1.3 5.8]);
            end;
            title(regionName{r});
        end;
        %set(gcf,'PaperPosition',[2 2 100 10]);wysiwyg;
    case 'ROI_Ttest_stats_untr-tr'   % Make a Ttest across different ROIs
        %df1_imana('ROI_Ttest_stats_tr-untr',T, 'con_0006', 1:8)
        T=varargin{1};      % ROI-stats structure
        dep=varargin{2};    % This is the dependent variable to plot
        reg= varargin{3};     % regions that should be ploted
        %----calculate z-values for 1 regressor version
        % T.acc= (T.acc-25)/76.5466;
        % T.accM= T.accM/76.5466;
        % T.accSubspaceM_prct= T.accSubspaceM_prct/76.5466;
        for r=reg
            %----first left than right hemisphere
            for h=1:2
                tC= getrow(T,T.hemisphere==h & T.region==r & T.seqCondition==1);
                uC= getrow(T,T.hemisphere==h & T.region==r & T.seqCondition==2);
                
                %----do table and report results
                fprintf('\n \n');
                fprintf('REGION: \t %s \nHEMISPHERE: \t %s \n', regionName{r}, hemName{h});
                fprintf('subj \t untrained \t trianed \t difference \t untr: #Voxel \t tr: #Voxel\n');
                for s=1:max(T.sn)
                    fprintf('%d \t %6.4f \t %6.4f \t %6.4f \t %6.4f \t %6.4f  \n',s, uC.(dep)(s,1), tC.(dep)(s,1), uC.(dep)(s,1)-tC.(dep)(s,1), 1, 1);%uC.numVoxel(s,1), tC.numVoxel(s,1))
                end
                fprintf('----------------------------------------------------------------\n')
                fprintf('mean \t %6.4f \t %6.4f \t %6.4f \n', mean(uC.(dep)), mean(tC.(dep)), mean(uC.(dep))-mean(tC.(dep)))
                fprintf('----------------------------------------------------------------\n')
                [t,p]=ttest(uC.(dep),tC.(dep),1,'paired');
                fprintf('================================================================\n')
            end
        end;
    case 'ROI_Ttest_stats_right-left'   % Make a Ttest across different ROIs
        %df1_imana('ROI_Ttest_stats_right-left',T, 'con_0006', 1:8)
        T=varargin{1};      % ROI-stats structure
        dep=varargin{2};    % This is the dependent variable to plot
        reg= varargin{3};     % regions that should be ploted
        %----calculate z-values for 1 regressor version
        % T.acc= (T.acc-25)/76.5466;
        % T.accM= T.accM/76.5466;
        % T.accSubspaceM_prct= T.accSubspaceM_prct/76.5466;
        for r=reg
            %----first left than right hemisphere
            for c=1:2
                lH= getrow(T,T.hemisphere==1 & T.region==r & T.seqCondition==c);
                rH= getrow(T,T.hemisphere==2 & T.region==r & T.seqCondition==c);
                
                %----do table and report results
                fprintf('\n \n');
                fprintf('REGION: \t %s Condition: \t %s \n', regionName{r}, condition{c});
                fprintf('subj \t left \t r \t difference \t untr: #Voxel \t tr: #Voxel\n');
                for s=1:max(T.sn)
                    fprintf('%d \t %6.4f \t %6.4f \t %6.4f \t %6.4f \t %6.4f  \n',s, rH.(dep)(s,1), lH.(dep)(s,1), rH.(dep)(s,1)-lH.(dep)(s,1), 1, 1);%uC.numVoxel(s,1), tC.numVoxel(s,1))
                end
                fprintf('----------------------------------------------------------------\n')
                fprintf('mean \t %6.4f \t %6.4f \t %6.4f \n', mean(rH.(dep)), mean(lH.(dep)), mean(rH.(dep))-mean(lH.(dep)))
                fprintf('----------------------------------------------------------------\n')
                [t,p]=ttest(rH.(dep),lh.(dep),1,'paired');
                fprintf('================================================================\n')
            end
        end;
        
        %----can be deleted
    case 'ROI_plotAsymmetry'
        %df1_imana('ROI_plotAsymetrie',D, 1, 5)
        D= varargin{1};                     % ROI structure
        reg= varargin{2};                   % the region that you like to ananlyse
        cat_threshold= varargin{3};    % the quantiles that you like to see => there are 10 quantiles, input 9 would mean the 20% highest
        %reduce struct to the region that we like to analyse
        C= getrow(D,D.region==reg);
        
        % You can use Split data with SPLIT, so you do not have to loop
        % over subjects or anything
        % For BOLD use: mean(C.mean,2)
        C.accSubspace_cat = split_data(C.accSubspace,'numquant',10,'split',[C.sn C.seqCondition]);
        
        subplot(4,1,2);
        X= tapply(C,{'sn' 'hemisphere', 'accSubspace_cat' ,'seqCondition'},...
            {'accSubspace','length','name','accSubspace_num'},...
            {'accSubspace','mean','name','accSubspace_mean'});
        
        
        F= tapply(X,{'sn','seqCondition'},...
            {'accSubspace_num','sum','name','accSubspace_numL','subset',X.hemisphere==1 & X.accSubspace_cat>=cat_threshold},...
            {'accSubspace_num','sum','name','accSubspace_numR','subset',X.hemisphere==2 & X.accSubspace_cat>=cat_threshold},...
            {'accSubspace_num','sum','name','numL','subset',X.hemisphere==1},...
            {'accSubspace_num','sum','name','numR','subset',X.hemisphere==2});
        
        % Get asymmetry scores, normalized by the size of each region
        F.asymN=(F.accSubspace_numR./F.numR)./(F.accSubspace_numR./F.numR+F.accSubspace_numL./F.numL);
        
        % Plot asym scores for trained and untrained
        subplot(2,1,1);
        myboxplot(F.seqCondition,F.asymN);
        set(gca,'YLim',[0 1]);
        drawline(0.5,'dir','horz');
        
        % Now against each other
        subplot(2,1,2);
        Y=pivottable(F.sn,F.seqCondition,F.asymN,'mean');
        scatterplot(Y(:,1),Y(:,2));
        xlabel('Trained');ylabel('Untrained');
        
        % Test against symmetry
        fprintf('Trained\n');
        ttest(Y(:,1)-0.5,[],2,'onesample');
        fprintf('Untrained\n');
        ttest(Y(:,2)-0.5,[],2,'onesample');
        fprintf('Train_untrained\n');
        ttest(Y(:,1),Y(:,2),2,'paired');
    case 'mv_files'%___________________________________________
        % df1_imana('mv_files',9:11, 'untrained')
        sn=varargin{1}; seqCondition=varargin{2};
        for i=sn
            %----source
            %s1= fullfile(baseDir,'glm_firstlevel_1', subj_name{i}, seqCondition, 'volsearch.mat');
            %s2= fullfile(baseDir,'glm_firstlevel_1', subj_name{i}, seqCondition, 'mask_noskull.nii');
            s3= fullfile(baseDir,'glm_firstlevel_1', subj_name{i}, seqCondition, 'lda_radius12.nii');
            s4= fullfile(baseDir,'glm_firstlevel_4', subj_name{i}, seqCondition, 'lda_radius12.nii');
            %---target
            %t1= fullfile(baseDir,'glm_firstlevel_4', subj_name{i}, seqCondition, 'volsearch.mat');
            %t2= fullfile(baseDir,'glm_firstlevel_4', subj_name{i}, seqCondition, 'mask_noskull.nii');
            t3= fullfile(baseDir,'glm_firstlevel_1', subj_name{i}, seqCondition, 'withPW_lda_radius12.nii');
            t4= fullfile(baseDir,'glm_firstlevel_4', subj_name{i}, seqCondition, 'withPW_lda_radius12.nii');
            %---copy
            %copyfile(s1, t1);
            copyfile(s3, t3);
            copyfile(s4, t4);
        end;    case 'variance_components_MT'
        T=varargin{1};
        B=varargin{2};
        %----substract the mean response of mt
        mt= B.MT_median- mean(B.MT_median);
        %----make design matrix
        Zgood= [kron(ones(4*8,1),eye(4))];
        Zfact= mt;
        Ac{1}=blockdiag(eye(4),0);
        Ac{2}=blockdiag(zeros(4,4),1);
        %----Decompose the data into variance components
        Z=[Zgood Zfact];
        [G,h,u,l,n,jumpI,b]=mvpattern_covcomp(T.beta',Z,'Ac',Ac,'meanS',0,'num_iter',1000,'Ac',Ac,'TolL',1e-4,'accel_method','Aitken');
        D.decompCoeff_seq= h(1,end);
        D.decompCoeff_mt= h(2,end);
        D.decompCoeff_noise= h(3,end);
        varargout={D};
    case 'regress_behav'
        
        %         case 'ROI_corr'%
        %         %T= df1_imana('ROI_corr')
        %
        %         R= load(fullfile(regDir,sprintf('reg_data.mat')));
        %
        %         D=[];
        %         for h=1:2
        %             for r=1:numel(unique(R.region(R.hemisphere==h)))
        %                 for s=1:max(R.sn)
        %
        %                     C=[];
        %                     %----region subset
        %                     rindx=find(R.sn==s & R.region==r  & R.hemisphere==h & ~R.isNaN);
        %                     X=getrow(R,rindx);
        %                     %----normal corr
        %                     %C= df1_imana('ROI_calc_covVar', X.mean)
        %                     %----decomp
        %                     numblocks=8;
        %                     cc=repmat([1:4 1:4],1,numblocks);
        %                     hh=repmat([1 1 1 1  2 2 2 2 ],1,numblocks);
        %                     b=[];
        %                     for i=1:numblocks
        %                         b=[b ones(1,8)*i];
        %                     end;
        %                     C= df1_decomp(X.beta,cc,hh,b);
        %                     %----
        %                     C.sn= s;
        %                     %----calc correlations
        %                     C.corrO= C.covO/sqrt(C.varO);
        %                     C.corrE= C.covE/sqrt(C.varE+0.00000001);
        %                     C.corrI= C.covI/sqrt(C.varI+0.00000001);
        %
        %                     D= addstruct(D, C);
        %
        %                     fprintf('hem: %i \t subj: %i \t region: %s \n ',h, s,char(unique(X.regionName)));
        %                 end
        %                 figure
        %                 %myboxplot([ones(size(D.corrO)); ones(size(D.corrO))*2], [D.corrI; D.corrE]-[D.corrO; D.corrO])
        %                 myboxplot([ones(size(D.corrO)); ones(size(D.corrO))*2], [D.corrI; D.corrE])
        %                 title(char(unique(X.regionName)));
        %                 fprintf('intrinsic \n')
        %                 ttest(D.corrI, [], 1, 'onesample')
        %                 fprintf('extrinsic \n')
        %                 ttest(D.corrE, [], 1, 'onesample')
        %                 waitforbuttonpress
        %                %keyboard
        %                 D=[];
        %             end;
        %         end;
        %
    otherwise
        error('no such case!')
end


function COG=sl1_COG(data,coord)
i=find(~isnan(data));
data=data(i,:);coord=coord(i,:);
low=0; % prctile(data,10);
d=data-low;
d(d<0)=0;
COG(:,1)=sum(coord(:,1).*d)./sum(d);
COG(:,2)=sum(coord(:,2).*d)./sum(d);
COG(:,3)=sum(coord(:,3).*d)./sum(d);


function C= df1_decomp_struct(y,c,h,b) %var=

Kb=8;K=4;B=8;
a=0.95;

[P,N]=size(y);
tt=c+(h-1)*4;
y=y(~isnan(y(:,1)),:);
[y,Sw]=mva_prewhiten(y,tt,'diagonal',1);

Zc=[kron(eye(2),ones(4,1)), [eye(8)]];
Z=repmat(Zc,8,1);
Z=[Z kron(eye(8),ones(8,1))];
Ac={};
Ac{1}=blockdiag([1 0;0 0],zeros(Kb,Kb),zeros(B,B));                                         % Var_a_left
Ac{2}=blockdiag([0 0;0 1],zeros(Kb,Kb),zeros(B,B));                                         % Var_a_right
Ac{3}=blockdiag(zeros(2,2),eye(K),zeros(K,K),zeros(B,B));                                   % var_b_left
Ac{4}=blockdiag(zeros(2,2),zeros(K,K),eye(K),zeros(B,B));                                   % var_b_right
Ac{5}=Ac{1}+Ac{2}+a*blockdiag((1-eye(2)),zeros(Kb,Kb),zeros(B,B));                          % Cov_a
Ac{6}=blockdiag(zeros(2,2),eye(Kb)+a*(diag(ones(K,1),K)+diag(ones(K,1),-K)),zeros(B,B));    % Cov_b_intrinsic
Ac{7}=blockdiag(zeros(2,2),[[eye(K);a*[0 1 0 0 ; 1 0 0 0; 0 0 0 1; 0 0 1 0]],[a*[0 1 0 0 ; 1 0 0 0; 0 0 0 1; 0 0 1 0]; eye(K)]],zeros(B,B));  % Cov_b_extrinsic
Ac{8}=blockdiag(zeros(2,2),zeros(Kb,Kb),eye(B));


[G,h,u,l,n,jumpI]=mvpattern_covcomp(y',Z,'Ac',Ac,'meanS',1,'num_iter',1000,'Ac',Ac,'TolL',1e-4,'accel_method','Aitken');
C.d_var_left= G(1,1);
C.d_var_right= G(2,2);
C.d_cov_condition= G(1,2);
C.d_var_sequences_left= G(3,3);
C.d_var_sequences_right= G(7,7);
C.d_covI= G(7,3);
C.d_covE= G(8,3);
C.d_noise= h(end, end);
% var(1,8)=h(end,end);
% var(1,1)=G(1,1)./var(1,8)*100;
% var(1,2)=G(2,2)./var(1,8)*100;
% var(1,3)=G(1,2)./var(1,8)*100;
% var(1,4)=G(3,3)./var(1,8)*100;
% var(1,5)=G(8,8)./var(1,8)*100;
% var(1,6)=G(3,8)./var(1,8)*100;
% var(1,7)=G(13,13)./var(1,8)*100;

function  out = df1_calcAcc_pilot(LDAinput,r,h,s,d)
% function  out = sl1_calAcc(LDAinput, varargin)
% LDAinput is a P x N input matrix (P voxels, N trials)
LDAinput= LDAinput(~isnan(sum(LDAinput, 2)), :);

if isempty(LDAinput)
    accuracy= NaN;
    numVox= NaN;   
else
    [U,S,V]=svd(LDAinput,0);
    for i=1:max(r)
        testI{i}=find(r==i);%find(b(indx)==i);
        trainI{i}=find(r~=i);%find(b(indx)~=i);
    end;
    accuracy =crossval_takemultipleout(@classify_lda_KclassesQuicker,S*V',d',trainI,testI);
end
out= [accuracy*100; size(LDAinput, 1)];

