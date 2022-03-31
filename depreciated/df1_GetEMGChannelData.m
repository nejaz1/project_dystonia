function D = df1_GetEMGChannelData ( fileName, fig, varargin )
% Fucntion uses the subj_name to return the EMG traces for given experiment
% The channels returned are the ones specified in the emgChan string vector
% nChans are returned and depending on the fig flag, the results are
% plotted
% indexes of channels to be used is [1-nChans-1], where the channel index
% nChan is always the TTL pulse to be synchronized to

nChans = 15;
Fs = 1000;
threshold = 3e4;
T = 3.5;
numTrials = 15;
vararginoptions(varargin,{'nChans', 'Fs', 'threshold','T','numTrials'}); 


M = load(fileName);

fields = fieldnames(M);
for i = 1:nChans-1
    D.chanName{i} = eval(['M.head' num2str(i) '.title']);
end;
D.chanName = D.chanName';

L = length(eval(['M.chan' num2str(nChans)]));       % length of the TTL vector
t = (1:1:L)/Fs;                                     % calculating time vector


ttl = RemoveTails(eval(['M.chan' num2str(nChans)]) >= threshold);   % times of the TTL pulses sent to 
                                                                    % start emg acquisition
                                                                    
ttlIndex = find(ttl==1);                            % find indexes where there is an onset of a TTL signal
nTTL = length(ttlIndex);

% creating fields in the structure object
for i = 1:nChans-1
    D = setfield(D,['chan' sprintf('%2.2d',i)],[]);
end;
D = setfield(D,'chan_t',[]);
D.ttlFlag = [];               % specifies whether all TTL's were obtained in this block

% for each ttl pulse, pick up 0s preceeding and window s post emg signal in each
% channel
%T = 5;
window = T*Fs;
for it = 1:nTTL
    
    % looping over channels
    for nc = 1:nChans-1
        c = eval(['M.chan' num2str(nc)]);
        
        try
            t = 1:1:length(c);
            c = double(c);
            c = c(ttlIndex(it):ttlIndex(it)+window-1);
%             disp([num2str(it) '-' num2str(nc) ': ' num2str(length(c))]);
        catch
%             disp([num2str(it) '-' num2str(nc) ': Too short, padding with zeros']);
            c = c(ttlIndex(it):end);
            short = abs(length(c) - window);
            c = [c; zeros(short,1)];
        end;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EMG signal filtering
        
        % bandpass filtering the emg signals between 20-250 Hz with a 4th
        % order butterworth filter
        c = c - mean(c);                % dc correction
        c = abs(c);
        
        FsHalf = Fs/2;
        [b,a] = butter(4,40/FsHalf,'low');
        filt_c = filtfilt(b,a,c);
        
        % downsampling the emg signal into 5 ms resolution
        NewFs = 200;
        filt_c = downsample(filt_c,Fs/NewFs);
        t = 0:1/NewFs:T-1/NewFs;
 
%         figure;
%         plot(c,'b');
%         hold on;
%         plot(filt_c,'r');
%         pause;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        cD = getfield(D,['chan' sprintf('%2.2d',nc)]);
        D = setfield(D,['chan' sprintf('%2.2d',nc)],[cD; filt_c']);
    end;

    cT = getfield(D,['chan_t']);
    D = setfield(D,['chan_t'],[cT; t]);
end;

% -- HACK - ASK TOBI
% to fix strange error in the dystonia project data
% if there are fewer ttl pulses in the channels, then assume that the
% problem is as a result of the last few trials and set them to zero
diff = numTrials-nTTL;
if diff > 0
    D.ttlFlag = zeros(numTrials,1);
    for iD = 1:diff
        for nc = 1:nChans-1
            cD = getfield(D,['chan' sprintf('%2.2d',nc)]);
            D = setfield(D,['chan' sprintf('%2.2d',nc)],[cD; t*0]);
        end;

        cT = getfield(D,['chan_t']);
        D = setfield(D,['chan_t'],[cT; t]);
    end;
else
    D.ttlFlag = ones(numTrials,1);
end;



if fig > 0
    figure;
    for i = 1:nChans-1
        subplot(3,5,i);
        c = eval(['M.chan' num2str(i)]);
        hold on;
        plot(t,ttl*threshold,'r');
        plot(t,c(1:L));
        title(D.chanName{i});
        xlabel('time (s)');
        ylabel('Arbitrary units');
    end;
end;


function D = FixEMGChannelNames ( M, nChans )
% Function fixes the names of the channels when they are not in enumerate
% order
f = fieldnames(M);

chanNum = [];
for i = 1:length(f)
    if strcmp(f{i}(1:4),'head')     % header channel
        chanNum = [chanNum str2double(strtok(f{i},'head'))];
    end;
end;

% fixing the fields names with numbering in the correct enumerate order
D = struct;
for i = 1:nChans
    tempChan = getfield(M,['chan' sprintf('%0.0d',chanNum(i))]);
    tempHead = getfield(M,['head' sprintf('%0.0d',chanNum(i))]);
    
    D = setfield(D,['chan' sprintf('%0.0d',i)],tempChan);
    D = setfield(D,['head' sprintf('%0.0d',i)],tempHead);
end;





