function D=df1_handemg_trial(MOV,D,fig,varargin)
% Function analyzes the subject/block from the single finger EMG finger
% individuation task

vararginoptions(varargin,{'startthres','endthres','minstart','maxstop','E','threshold','nChan'});

% Extract data:
if isempty(MOV)
    return;
end;


% defining the chords that are used for the experiment
chords=[eye(5); 1,1,0,0,0; 1,0,1,0,0; 1,0,0,1,0; 1,0,0,0,1; 0,1,1,0,0;...
    0,1,0,1,0; 0,1,0,0,1; 0,0,1,1,0; 0,0,1,0,1; 0,0,0,1,1;...
    0,0,1,1,1; 0,1,0,1,1; 0,1,1,0,1; 0,1,1,1,0; 1,0,0,1,1;...
    1,0,1,0,1; 1,0,1,1,0; 1,1,0,0,1; 1,1,0,1,0; 1,1,1,0,0;...
    0,1,1,1,1; 1,0,1,1,1; 1,1,0,1,1; 1,1,1,0,1; 1,1,1,1,0;...
    1,1,1,1,1;];

% (1) Saving forces and state of the experiment
state = MOV(:,1);                       % state of the stimulus presentation
t = MOV(:,2);                           % time for the measurement of force
screen_t = MOV(:,3);                    % time for which screen flips are obtained
Force = MOV(:,4:13);                    % force data corresponding to L1-5 and R1-5 fingers
D.ttlFlag = E.ttlFlag;

Force = smooth_kernel(Force,5);         % Smoothing with Gaussian kernel 
D.maxforce = max(Force);                  % calculating the maximum force over the trial for all fingers
D.relmaxforce = D.maxforce./D.mvc;        % scaling the maximum force over the trial with the mvc
D.goalForce = chords(D.digit,:);



% (2) Calculating the active and the passive finger indexes
ActFin = getchordindex(D.hand,D.digit);
if D.hand == 1
    ActHand = 1:5;
else
    ActHand = 6:10;
end;
% PassFin = setdiff(1:10,ActFin);
% PassHand = 1:5;


% (3) Calculating the required forces in newtons
D.lowForceN = D.lowForce * D.mvc;
D.highForceN = D.highForce * D.mvc;
D.targetForceN = D.targetForce * D.mvc;


% (4) calculating the instance when the finger starts to move
mvc = repmat(D.mvc,length(Force),1);
dF = velocity_discr(Force./mvc,3);           % Derivative change in force, smoothed with a 3 point gaussian
dF_digit = dF(:,getfingerindex(D.hand,D.digit));     % derivative for the hand/digit combination being tested for in the trial

stim_i = find(state>1,1,'first');               % calculate when the stimulus was presented
press_i = find(dF_digit>=threshold,1,'first');  % calculate when the finger is first pressed
release_i = size(Force,1);


% (5) updating the times in the D object to be returned
D.pressTime = t(press_i);
D.signalTime = t(stim_i);
D.releaseTime = t(release_i);
D.baselineForce = Force(stim_i,:); 

%(6) Define how many fingers were pressed in the instructed chord
if (D.digit>0 && D.digit<6)
    D.chordType=1; %1-finger press
elseif (D.digit>5 && D.digit<16)
    D.chordType=2; %2-finger chords
elseif (D.digit>15 && D.digit<26)
    D.chordType=3; %3-finger chords
elseif (D.digit >25 && D.digit<31)
    D.chordType=4; %4-finger chords
elseif (D.digit==31)
    D.chordType=5; %5-finger chord
end;

% (7) Compute different individuation measures
baselineF = Force(stim_i,ActHand);
stdF = bsxfun(@minus,Force(:,ActHand),baselineF);
thresholdFPress = min(D.lowForceN(ActHand));  % in the task, we use %mvc, so using the minimum of the absolute forces for the instructed fingers
stdF = bsxfun(@rdivide,stdF,thresholdFPress-baselineF);
stdF(stdF>1) = 1;

projL = stdF * D.goalForce'; % Scalar projection of the instructed fingers 5 dim space of finger activations
projL = projL./D.chordType;  % divide by sum of total active fingers used
res = stdF - projL*D.goalForce; % total residuals 

resL=sqrt(sum(res(stim_i:release_i,:).^2,2));   % residuals calculated over the press and release events
D.meanDev=nanmean(resL);
D.meanDevS=nanmean(resL)./sqrt(D.goalForce*D.goalForce');

resL=sqrt(sum(res(stim_i:release_i,find(D.goalForce==1)).^2,2));    % residuals over active fingers
D.meanDevA=nanmean(resL);

resL=sqrt(sum(res(stim_i:release_i,find(D.goalForce==0)).^2,2));    % residuals over passive fingers
D.meanDevP=nanmean(resL);

% ----- Marking those meanDev calculations above 3 as bad
if D.meanDev > 3
    D.meanDevFlag = 0;
else
    D.meanDevFlag = 1;
end;


% (8) calculating the eigenvalues/eigenvectors of the square normalized force matrix
% saving the highest eigenvalue (as a percentage of all other eigenvalues)
% and the corresponding eigenvector
S = Force(stim_i:release_i,:)'*Force(stim_i:release_i,:);
[~,EVal] = eig(S);
EVal = diag(EVal);
D.EVal = max(EVal)/sum(EVal);

% (9) cutting the emg and the force vectors and saving them
Fs = 1/mean(diff(t));
if isempty(press_i)
    D.pressFlag = 0;
    press_i = stim_i;
else
    D.pressFlag = 1;
end;
at = press_i;
pre = 0;                % align to when the finger is pressed
post = round(3.5 * Fs)-1; % 2  s post data

% post_step = round((D.releaseTime-D.pressTime-0.5)* Fs)-1; % 2  s post data
% Force_step = Force(:,getfingerindex(D.hand,D.digit));     % active hand
% Force_step = cut(Force_step,pre,at,post_step,'padding','last');
% t_step = cut(t,pre,at,post_step,'padding','last');
% s = stepinfo(Force_step',t_step,'SettlingTimeThreshold',0.2);
% D.RiseTime = s.RiseTime;
% D.SettlingTime = s.SettlingTime-t_step(1);
% D.Overshoot = s.Overshoot;

t_c = cut(t,pre,at,post,'padding','nan')';

% UNCOMMENT TO SAVE FORCE TRACES
% Force_c = cut(Force,pre,at,post,'padding','nan');
% D.r1 = Force_c(:,getfingerindex(2,1))';
% D.r2 = Force_c(:,getfingerindex(2,2))';
% D.r3 = Force_c(:,getfingerindex(2,3))';
% D.r4 = Force_c(:,getfingerindex(2,4))';
% D.r5 = Force_c(:,getfingerindex(2,5))';
% D.t = t_c;

D.emgMean = zeros(1,nChan);
D.emgMax = zeros(1,nChan);

chanT = cut(getfield(E,['chan_t'])',pre,press_i-stim_i,post);       % getting channel time so that I can trim emg signals                                                                    
for nc = 1:nChan
    chan = cut(getfield(E,['chan' sprintf('%2.2d',nc)])',pre,press_i-stim_i,post);
    
    % estimating the indexes that need to be trimmed
    trimIndex = isnan(t_c);
    lastVal = find(trimIndex==0,1,'last'); % last value for the emg
    chan(trimIndex) = chan(lastVal);
    
    % scaling the channel by the maximum force observed in that channel
    % during mvc
    if isfield(D,'emgMvc')
        chan = chan/D.emgMvc(nc);
    end;

    % UNCOMMENT TO SAVE EMG CHANNELS    
%     D = setfield(D,['chan' sprintf('%2.2d',nc)],chan');

    D.emgMean(nc) = mean(chan(~trimIndex));
    D.emgMax(nc) = max(chan(~trimIndex));
end;

    
                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display plot for trial (if required)
if (fig > 0)
    subplot(3,1,[1 2]); 
    plot(t,Force(:,1),'k-',...
        t,Force(:,2),'b-',...
        t,Force(:,3),'g-',...
        t,Force(:,4),'y-',...
        t,Force(:,5),'r-',...
        t,Force(:,6),'k--',...
        t,Force(:,7),'b--',...
        t,Force(:,8),'g--',...
        t,Force(:,9),'y--',...
        t,Force(:,10),'r--','LineWidth',3);    
    legend('l1','l2','l3','l4','l5','r1','r2','r3','r4','r5','Location','West');
    drawline(D.pressTime,'color',[0 0 0]); 
    drawline(D.signalTime,'color',[1 0 0]);
    drawline(D.releaseTime,'color',[0 0 1]);
    title(['hand ' num2str(D.hand) '; finger ' num2str(D.digit) '; VAF ' num2str(D.EVal)]);
    v = axis;
    subplot(3,1,3); 
    plot(t,state); 
    
    keyboard; 
end;
