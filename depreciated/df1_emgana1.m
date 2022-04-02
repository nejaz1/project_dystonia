function varargout = df1_emgana1 ( what, varargin )
% Suit of functions to extract information from the EMG experiment

nChans      =   16;                           % total number of channels (1 for trigger)
trigChan    =   17;                           % index for trigger channel
hand        =   [ones(1,8) ones(1,8)+1];      % channel location 
Fs          =   1000;                         % sampling rate
threshold   =   3e4;                          % threshold to detect trigger

baseDir=        fullfile('/Volumes/ANNA/data/FingerPattern_dystonia');
emgDir=         fullfile(baseDir, 'EMG/data');

emg_subj={      'd01','d02','d03','d04','d05',...
                'd06','d07','d08','d09','d10',...
                'd11',...
                's01','s02','s03','s04','s05',...
                's06','s08','s09'};

emg_NumRuns=[   10 10 10 10 11 ...
                10 10 10 10 10 ...
                10 ...
                10 10 10 10 10 ...
                10 11 10];
switch(what)
    
    case 'channelLabel'         % list out name of emg channels
        
        for s = 1:length(emg_subj);
        for sn = s
        
            % (1) check mvc file
            mvcFile   = fullfile(emgDir,emg_subj{sn},[emg_subj{sn} '_mvc.mat']);
            D=load(mvcFile);
                % channels per hand
                for i=1:nChans
                    T=getfield(D,['head' num2str(i)]);
                    fprintf('Hand%d %s\n',hand(i),T.title);
                end;
                % trigger channel
                if  isfield(D,['head' num2str(trigChan)])
                    fprintf('Trigger - yes\n');
                else 
                    fprintf('Trigger - no!\n');
                end;
        
            % (2) check individual runs 
            for i = 1:emg_NumRuns(sn)
            runFile = fullfile(emgDir,emg_subj{sn},[emg_subj{sn} '_run' num2str(i) '.mat']);
            D = load(runFile);
                % channels per hand
                for i=1:nChans
                    T=getfield(D,['head' num2str(i)]);
                    fprintf('Hand%d %s\n',hand(i),T.title);
                end;
                % trigger channel
                if  isfield(D,['head' num2str(trigChan)])
                    fprintf('Trigger - yes\n');
                else 
                    fprintf('Trigger - no!\n');
                end;
            end
        end
        end
        
   
        
    case 'countTrigger'         % count a list of triggers
        % df1_emgana('countTrigger','d01');
        subj = varargin{1};
        runFiles=dir([subj '_run*.mat']);   
        mvcFile =dir([subj '_mvc.mat']);   

        % (1) Triggers for mvc
        D=load(mvcFile(1).name);
        T=getfield(D,['chan' num2str(trigChan)]);       % get trigger channel
        ttl = RemoveTails(T >= threshold);              % detect trigger
        fprintf('MVC %d\n',sum(ttl));
        
        % (2) Triggers for runs
        for i=1:size(runFiles,1)
            D=load(runFiles(i).name);
            T=getfield(D,['chan' num2str(trigChan)]);       % get trigger channel
            ttl = RemoveTails(T >= threshold);              % detect trigger
            fprintf('Run%d %d\n',i,sum(ttl));
        end;        
    
    
    case 'plotTriggerSync'                                      % visualize emg data synchronization with trigger pulse
        % df1_emgana('plotTriggerSync','d01');
        subj = varargin{1};
        runFiles=dir([subj '_run*.mat']);   

        scrsz = get(0,'ScreenSize');
        figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/4])
        
        % (1) Triggers for runs
        for i=1% :size(runFiles,1)
            D=load(runFiles(i).name);
            ttl = getfield(D,['chan' num2str(trigChan)])';          % get trigger channel
            ttl = RemoveTails(ttl >= threshold);                      % detect trigger
            for j=1:nChans
                c         = getfield(D,['chan' num2str(j)])';
                chan(j,:) = double(c(1:length(ttl)));                  % get channel activity
            end;
            chan    = df1_emgana('preprocess',double(chan));           % emg preprocessing
            chan    = bsxfun(@rdivide,chan,std(chan,[],2));            % standardizing to unit variance
            rmsEmg  = sqrt(sum(chan.^2,1));                            % root-mean squared emg activity
            idx     = rmsEmg>4*std(rmsEmg);                            % outlier estimation
            rmsEmg(idx) = nan;                                         % removing outliers
            rmsEmg  = rmsEmg/max(rmsEmg);                              % scaling
            t       =(1:1:length(ttl))/Fs;                             % time
            ttlTime = t(ttl);                                          % time TTL pulse occurs
            
            % plotting ttl vector against rms EMG activity
            subplot(1,10,1:9);
            plot(t,ttl,t,rmsEmg);
            title(sprintf('Run %d',i));
            axis([t(1) t(end) -0.5 1.5]);
            subplot(1,10,10);
            plot(diff(ttlTime));
            drawline(2*median(diff(ttlTime)),'dir','horz','color','r')
            ylim([3 15]);
            title('Time (s) b/w trigger');
            clear chan
            keyboard;
        end;                
        
    case 'preprocess'           % emg preprocessing, mean subtraction, lowpass filter
        chans = varargin{1};        % emg channels N chans x T time points
        f_c   = 10;                 % cutoff frequency

        for i=1:nChans
            % lowpass filtering the emg signals at 40 Hz with a 4th
            % order butterworth filter
            c = chans(i,:);
            c = abs(c - mean(c));

            FsHalf = Fs/2;
            [b,a] = butter(4,10/FsHalf,'low');
            filt_c(i,:) = filtfilt(b,a,c);
        end;
        varargout={filt_c};
end;

