function R=df1_handemg_subj (subjIDName, fig, block, trial, varargin)
% Function analyzes the subject/session from the single finger EMG finger
% individuation task
% USAGE:
%       handemg1_subj('th270612',1,3,1);        For the subject, shows trial 1 of block 3
%       handemg1_subj('th270612',1);            For the subject, shows trial 1

vararginoptions(varargin,{'threshold','nChan'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing for input arguments
if nargin < 2                                       % no plotting
    fig=0;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R=[];           % output object
dataFile = ['IN2b_' subjIDName '.mat'];               % contains the block runs information
mvcFile = ['IN2mvc_' subjIDName '.mat'];              % the mvc file for the subject

% for the pilot data
% mvcDataFile = ['IN2mvc_' subjIDName '_01.mat'];              % the mvc file for the subject

% for the single finger data
mvcDataFile = ['IN2mvc_' subjIDName '_00.mat'];              % the mvc file for the subject

outFile = ['IN2b_' subjIDName '_ana.mat'];            % file to write out
mvcOutFile = ['IN2b_' subjIDName '_mvc.mat'];            % file to write out

D = load(dataFile);                                 % load the block information
D = test_block_integrity(D);                        % testing out the integrity of the 
                                                    % experiment and removing the bad trials

% Add the mvc into the data structure 
M = load(mvcFile);                                  % load the MVC file data
MD = load(mvcDataFile);                             % load the MVC file with force measurements
D.mvc = repmat(M.mvcf',length(D.BN),1);             % repeat the mvc calculated over the 
                                                    % two runs for each trial

                                                    
% loading the emg data
E = load(['emg_' subjIDName '_ana.mat']);
EMVC = load(['emg_' subjIDName '_mvc.mat']);
                                                    
% loading block data
if nargin < 3                                       % pick all blocks, all trials
    s = 1;
else                                                
    if nargin < 4
        s = find(D.BN==block & D.TN==1);            % if [Trial] is not provided,
                                                    % start at the [Block],[Trial=1] level
    else
        s = find((D.BN==block) & (D.TN==trial));    % if both [Block] and [Trial] are specified
                                                    % start at the [Block],[Trial] level
    end;
end;

% -------------------------------------------------------------------------
% mvc calculations for emg activations
M.mvc = repmat(M.mvcf',length(M.TN),1);
T = struct;

M = rmfield(M,'mvcf');
EMVC = rmfield(EMVC,'chanName');
E = rmfield(E,'chanName');

for m = 1:length(M.TN)
    mvc_C = df1_handemg_trial(MD.MOV{M.TN(m)},getrow(M,m),fig,'E',getrow(EMVC,m),'threshold',threshold,'nChan',nChan);
    T = addstruct(T,mvc_C,'row','force'); 
end;

mvc_E = struct;
mvc_E.emgMean = [];
mvc_E.emgMax = [];

for d = 1:length(unique(M.digit))
    di = find(M.digit == d);
    mvc_E.emgMean = [mvc_E.emgMean; mean(T.emgMean(di,:))];
    mvc_E.emgMax = [mvc_E.emgMax; max(T.emgMax(di,:))];
end;
D.emgMvc = repmat(max(mvc_E.emgMax),length(D.BN),1);
save(mvcOutFile,'-struct','mvc_E'); 


% -------------------------------------------------------------------------
% data calculations
oldblock=-1;            
disp(sprintf('Reading\n-------'));
for i = s:length(D.BN)                              % picking up the trials starting from 
                                                    % trial s to the last trial    
   if oldblock ~= D.BN(i)
        oldblock = D.BN(i);
        load(['IN2b_' subjIDName '_' num2str(D.BN(i),'%0.2d') '.mat']); % pick up data from fileName
        disp(['IN2b_' subjIDName '_' num2str(D.BN(i),'%0.2d') '.mat']);
   end;
%    keyboard;
   C = df1_handemg_trial(MOV{D.TN(i)},getrow(D,i),fig,'E',getrow(E,i),'threshold',threshold,'nChan',nChan);
   R = addstruct(R,C,'row','force'); 
end;

% emgMvc = zeros(1,nChan);
% for nc = 1:nChan
%     chan = getfield(R,['chan' sprintf('%2.2d',nc)]);
%     emgMvc(nc) = max(max(chan));
%     chan = chan/emgMvc(nc);
%     R = setfield(R,['chan' sprintf('%2.2d',nc)],chan);
% end;
% R.emgMvc = repmat(emgMvc,length(D.BN),1);

D = R;
save(outFile,'-struct','D'); 
