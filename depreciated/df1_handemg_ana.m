function varargout=df1_handemg_ana(what,varargin)
% 28/Jan/2013 - Created by Naveed Ejaz

% Directories
baseDir = '/Users/naveed/Documents/data/FingerPattern_dystonia/Individuation_EMG/data';
% baseDir = '/Volumes/MotorControl/data/FingerPattern_dystonia/Individuation_EMG/data';
analysisDir = '/Users/naveed/Documents/data/FingerPattern_dystonia/Individuation_EMG/analysis';

% Names
% s04 emg readings are corrupted
subj_name={'d01','d02','d03','d04','d05','d06','d07','d08','d09','d10','d11','s01','s02','s03','s05'};      
id = {'021112', '031212', '021112', '021112', '071112', '071112', '151112', '161112', '031212', '191112', '191112', '141212', '101212', '011112', '151112'};
numBlocks=[10 10 10 10 11 10 10 10 10 10 10 10 10 10 10 10];

numTrials = 30;        % number of runs per trial
mvcTrials = 20;
numDigits = 5;
nChan = 8;            % number of emg channels recorded from one hand (x2)
nForce = [0.25, 0.50, 0.75];
emgChan = {'R APB','R 2Lum','R 3Lum','R 4Lum','R ADM','R FDI','R Ext','R Flex', ...
           'L APB','L 2Lum','L 3Lum','L 4Lum','L ADM','L FDI','L Ext','L Flex'};

% analysis procedures
switch(what)
    % ----------------------------------------------------------------------------------------------
    % Pre-processing steps
    case 'process'
        sn=varargin{1};
        df1_handemg_ana('smr2mat',sn);
        df1_handemg_ana('process_emg_data',sn);
        df1_handemg_ana('behavioural_data',sn);
        
    case 'smr2mat' % rename and convert the CED smr files (EMG data) into matlab files
        % this procedure requires the sigTOOL toolbox
        s=varargin{1};
        
        for sn = s
            cd(fullfile(baseDir,subj_name{sn}));    % subject data directory
            
            % finding the latest mvc smr file
            mvc = dir('*MVC*.smr');
            mvc = sort({mvc(:).name});
            mvcFile = mvc{end};
            disp(mvcFile);

            % convert and rename the mvc file
            mvcFile = fullfile(mvcFile);    
            fid = fopen(mvcFile);
            SONImport(fid);
            fclose(fid);
            disp('Mvc file converted');
            try
                [~,mvcFile,~] = fileparts(mvcFile);
                source = fullfile(baseDir,subj_name{sn},[mvcFile '.mat']);    
                dest = fullfile(baseDir,subj_name{sn},[subj_name{sn} '_emg_mvc.mat']);  
                movefile(source,dest);
            catch
                disp(['MVC file could not be moved']);
            end;

            % move and rename the individual runs
            for i = 1:numBlocks(sn)
                dataFile = fullfile(baseDir,subj_name{sn}, [subj_name{sn} '_' id{sn} '_run' num2str(i) '.smr']);
                fid = fopen(dataFile);
                SONImport(fid);
                fclose(fid);
                disp(['Block ' num2str(i) ' converted']);

                try
                    source = fullfile(baseDir,subj_name{sn}, [subj_name{sn} '_' id{sn} '_run' num2str(i) '.mat']);
                    dest = fullfile(baseDir,subj_name{sn}, [subj_name{sn} '_r' num2str(i) '.mat']);
                    movefile(source,dest);
                catch exception
                    disp(['Block ' num2str(i) ' could not be moved']);
                end;

            end;
        end;
    
    case 'process_emg_data' % make structures for the signals captured form the EMG data
        s=varargin{1};

        nChans = 2*nChan+1;            % one additional channel for the TTL pulse
        Fs = 1000;
        threshold = 2.5e4;
        T = 5.5;
        plotFlag = 0;

        for sn=s 

            cd([baseDir '/' subj_name{sn}]);

            % process the mvc file
            D = df1_GetEMGChannelData([subj_name{sn} '_emg_mvc.mat'],plotFlag,'nChans',nChans,'Fs',Fs,'threshold',threshold,'T',T,'numTrials',mvcTrials);
            outFile = ['emg_' subj_name{sn} '_mvc.mat'];
            save(outFile,'-struct','D'); 
            
            D = struct;
            disp('');
            for i = 1:numBlocks(sn)
                block = df1_GetEMGChannelData([subj_name{sn} '_r' num2str(i) '.mat'],plotFlag,'nChans',nChans,'Fs',Fs,'threshold',threshold,'T',T,'numTrials',numTrials);
                disp([subj_name{sn} '_r' num2str(i) ' B: ' num2str(size(block.chan01,1))]); 
                D = addstruct(D,block);
            end;        
            outFile = ['emg_' subj_name{sn} '_ana.mat'];
            save(outFile,'-struct','D'); 
        end; 

    case 'behavioural_data' % make structures for the force levels calculated from the behavioural data
        s=varargin{1};

        trial = 1;
        block = 1;
        fig = 0;
        threshold = 0.0015;      % the single finger (non pilot) requires a slightly larger threshold for detecting finger movements
        emgChans = nChan*2;

        for sn = s

            % make structures for the behavioural data
            cd([baseDir '/' subj_name{sn}]);
            df1_handemg_subj(subj_name{sn},fig,block,trial,'threshold',threshold,'nChan',emgChans);

            % moving the analysis file
            src = ['IN2b_' subj_name{sn} '_ana.mat'];
            dest = [subj_name{sn} '_ana.mat'];
            movefile(src,dest);
            copyfile(dest,fullfile(analysisDir,dest));

            % moving the mvc file
            src = ['IN2b_' subj_name{sn} '_mvc.mat'];
            dest = [subj_name{sn} '_mvc.mat'];
            movefile(src,dest);
            copyfile(dest,fullfile(analysisDir,dest));
        end;
        
    case 'make_alldat' % make an alldat structure for the combined information across all subjects
%         s = 1:length(subj_name);
        s = varargin{1};
        
        R = struct;
        for sn = s
            src = fullfile(analysisDir,[subj_name{sn} '_ana.mat']);
            D = load(src);
            
            % adding subject number
            D.SN = sn + zeros(length(D.BN),1);
            
            R = addstruct(R,D,'row','force'); 
        end;
        outFile = fullfile(analysisDir,'alldat.mat');
        save(outFile,'-struct','R'); 

    % ----------------------------------------------------------------------------------------------
    % Analysis section    
    case 'emg_force_scaling_avg'  % for each muscle, look at the scaling for each of the three force levels
        % dystonics:    df1_handemg_ana('emg_force_scaling_avg',[1:4 7:9 11])
        % musicians:    df1_handemg_ana('emg_force_scaling_avg',[13 15])
        
        sn=varargin{1};
        scaling = [-5 -0.5];
        if size(varargin) < 2 
            fig = 1;
        else
            fig = varargin{2};
        end;
        
        MS.f1 = zeros(nChan*2,numDigits);
        MS.f2 = zeros(nChan*2,numDigits);
        MS.f3 = zeros(nChan*2,numDigits);
        
        for s = sn
            
            cd(analysisDir);
        
            src = [subj_name{s} '_ana.mat'];
            D = load(src);
            
            for f = 1:length(nForce)

                muscleScaling = [];
                for d = 1:5
                    i = find((D.digit==d) & (D.targetForce==nForce(f)));
                    R = getrow(D,i);
                    muscleScaling = [muscleScaling, mean(R.emgMean)'];
                end;
                
                ms = getfield(MS,['f' num2str(f)]);
                MS = setfield(MS,['f' num2str(f)],ms+muscleScaling);
            end;
        end;
        
        % getting heatmaps for muscle-force mapping
        f1 = getfield(MS,'f1');
        f2 = getfield(MS,'f2');
        f3 = getfield(MS,'f3');
        f1 = f1/length(sn);
        f2 = f2/length(sn);
        f3 = f3/length(sn);
        MS = setfield(MS,'f1',f1);
        MS = setfield(MS,'f2',f2);
        MS = setfield(MS,'f3',f3);
        
        if fig == 1
%             scrsz = get(0,'ScreenSize');
%             figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
            subplot(1,3,1);
            imagesc(log(f1),scaling);
            title(['Force = ' num2str(nForce(1))]);
            colorbar;
            set(gca,'YTick',1:nChan*2);
            set(gca,'YTickLabel',emgChan);
            set(gca,'XTick',1:numDigits);
            subplot(1,3,2);
            imagesc(log(f2),scaling);
            title(['Force = ' num2str(nForce(2))]);
            colorbar;
            set(gca,'YTick',1:nChan*2);
            set(gca,'YTickLabel',emgChan);
            set(gca,'XTick',1:numDigits);
            subplot(1,3,3);
            imagesc(log(f3),scaling);
            title(['Force = ' num2str(nForce(3))]);
            colorbar;
            set(gca,'YTick',1:nChan*2);
            set(gca,'YTickLabel',emgChan);
            set(gca,'XTick',1:numDigits);
        end;
        
        varargout = {MS};
    
    case 'mahalanobis_distance_comp'  % calculates and plots the mahalanobis distance
        sn=varargin{1};
        scaling = [0 15];
        
        Cavg = zeros(10,1);

        for s = sn
            
            src = fullfile(analysisDir,[subj_name{s} '_ana.mat']);
            D = load(src);
            idx = find(D.ttlFlag == 1);
            D = getrow(D,idx);
                
                
            % joerns functions
%             C = distance_mahalanobis(D.emgMean',[(D.hand-1)*5 + D.digit]');
            C = distance_mahalanobis(D.emgMean',D.digit');
%             C = squareform(C);

            Cavg = Cavg + C;
        end;

        % averaging over subjects
        Cavg = Cavg/length(sn);
        
        % for plotting barplot, putting nan along the main diagonal of the
        % covariance matrix
        Cavgx = SetMatrixValues(Cavg,'trace');

%         subplot(2,4,s+1)
        barplot([],Cavgx,'facecolor',[0.3 0.3 0.3]);
        set(gca,'XTickLabel',num2str([1:2*numDigits]'));
        title('Avg MD');
        
        varargout = {Cavgx};
        

end;

