function D=df1_trial(MOV,D,fig,varargin)
% Analysis of a subject from Individuation2b - bedside version task 
% 
forceThreshold=2.5; % This is release and start threshold; 

% -------------------------------------------------------------------------
% OPTIONS:
vararginoptions(varargin,{'startthres','endthres','minstart','maxstop'});

% -------------------------------------------------------------------------
% Extract data:
if (isempty(MOV))
    return;
end;
state=MOV(:,1);
t=MOV(:,2);
screen_t=MOV(:,3);      % Screentime 

Force=MOV(:,4:13);              %cols 4-13 are the Force data
Force=smooth_kernel(Force,5);   % Smoothing with Gaussian kernel 

% correcting for older experiments which have force instead of target force as a field
if isfield(D,'force')
    tForce = getfield(D,'force');
    D = rmfield(D,'force');
    D = setfield(D,'targetForce',tForce);
    clear tForce;
end;


% -------------------------------------------------------------------------
% cleaning the data
% 1) calculating the start of the trial
% 2) only using Force measurements from the start for analysis
% 3) subtracting the baseline from the Force measurements
%       - baseline is defined as the force at the start of the trial
start=find(state==2,1,'first'); 
D.baselineForce=Force(start,:); 
normF=bsxfun(@minus,Force,D.baselineForce);  % Get the force from the active phase - baseline 
mvc = max(D.mvc,1);                          % if mvc for any finger is less than one, it is set to one
% normF=bsxfun(@rdivide,normF,mvc);          % Do NOT standardize to MVC 
normF = normF(start:end,:);                  % picking out a subset of the dataset

% calculating the eigenvalues/eigenvectors of the square normalized force matrix
% saving the highest eigenvalue (as a percentage of all other eigenvalues)
% and the corresponding eigenvector
S = Force(start:end,:)'*Force(start:end,:);
[EVec,EVal] = eig(S);
EVal = diag(EVal);
[EVal,sIndx] = sort(EVal,'descend');
EVec = EVec(sIndx);
D.EVal = EVal(1)/sum(EVal);
D.EVec = EVec(:,1)';

D.goodTrial = (D.EVal > 0.95);     % currently a hard threshold is used to distinguish
                                   % between good and bad trials

% Define indices into the force matrix 
IndActHand=(D.hand-1)*5+[1:5]; 
IndPassHand=(1-D.hand)*5+[1:5]; 
IndActDigit=(D.hand-1)*5+D.digit; 
   
% -------------------------------------------------------------------------
% Do mixing matrix calculation 
dF = velocity_discr(normF,3);           % Derivative change in force, smoothed with a 3 point gaussian
CovF=cov(dF);       % Calculate variance-covariance matrix 
CorrF=corrcoef(dF);
D.CovF=CovF(:)';    % Spread out those 25 numbers into 1 row vector we can store in the data file
% [~,~,pcomp] = princomp(normF);      % calculating the principal components
% D.princomp = (pcomp/sum(pcomp))';

D.maxforce=max(Force);                  % saving peak force
D.actForce=D.maxforce(IndActDigit);     % Max force on active finger 
D.relmaxforce=D.maxforce./mvc;          % force scaled by the mvc for every finger

% -------------------------------------------------------------------------
% Calculate measures of digit independence: deviation from straight path 
normFAct=normF(:,IndActHand);   % Active hand only 
forceA=normFAct(:,D.digit); 
passIndx=[1:5];
passIndx(D.digit)=[]; 
forceP=normFAct(:,passIndx);
D.meanDevF=sqrt(mean(normFAct.^2)); 
D.meanDevP=mean(sqrt(sum(forceP.^2,2))); % Mean Eucledian deviation from straight 
D.meanDevA=mean(sqrt(sum(forceA.^2,2))); % Mean movement along straight line 
D.meanDevDiff=D.meanDevA-D.meanDevP;

for i=1:size(mvc,2)
    relnormF(:,i)=normF(:,i)./mvc(i);
end
relnormFAct=relnormF(:,IndActHand);   % Active hand only 
relforceA=relnormFAct(:,D.digit); 
relforceP=relnormFAct(:,passIndx);
D.relmeanDevP=mean(sqrt(sum(relforceP.^2,2))); % Mean Eucledian deviation from straight 
D.relmeanDevA=mean(sqrt(sum(relforceA.^2,2))); % Mean movement along straight line 
D.relmeanDevDiff=D.relmeanDevA-D.relmeanDevP;

% Recalculate the timing of the presses 
c=D.mvc((D.hand-1)*5+D.digit);
D.targetForceAbs=0.8*c*D.targetForce; 
D.lowForceAbs=0.8*c*D.lowForce; 
D.highForceAbs=0.8*c*D.highForce; 


pressIndx=find(Force(:,IndActDigit)>D.lowForceAbs,1,'first');
startIndx=find(Force(:,IndActDigit)>forceThreshold,1,'first');
endholdIndx=find(state==4,1,'first');
endIndx=findend(Force(:,IndActDigit)>D.lowForceAbs,pressIndx);
releaseIndx=findend(Force(:,IndActDigit)>forceThreshold,pressIndx);

if (isempty(pressIndx))
    D.pressTime=NaN;
else
    D.pressTime=t(pressIndx); 
end
if (isempty(startIndx))
    D.startTime=NaN;
else
    D.startTime=t(startIndx); 
end
if (isnan(endIndx))
    D.endTime=NaN;
else
    D.endTime=t(endIndx); 
end
if (isnan(releaseIndx)) 
    D.releaseTime=NaN; 
else
    D.releaseTime=t(releaseIndx);
end; 

% Accuracy and speed measures 
D.riseTime=D.pressTime-D.startTime; 
D.dropTime=D.releaseTime-D.endTime; 
D.rmsAct=mean((Force(pressIndx:endholdIndx,IndActDigit)-D.targetForceAbs).^2);

% -------------------------------------------------------------------------
% Display trial
if (fig>0)
    subplot(4,1,[1 2]); 
    plot(t,Force(:,1),'k-',...
        t,Force(:,2),'b-',...
        t,Force(:,3),'g-',...
        t,Force(:,4),'r-',...
        t,Force(:,5),'y-',...
        t,Force(:,6),'k--',...
        t,Force(:,7),'b--',...
        t,Force(:,8),'g--',...
        t,Force(:,9),'r--',...
        t,Force(:,10),'y--','LineWidth',3);    
    leg = {'l1','l2','l3','l4','l5','r1','r2','r3','r4','r5'};
    legend(leg,'Location','West');
    drawline(D.signalTime,'color',[1 0 0]); 
    drawline(D.startTime,'color',[0 0 0]); 
    drawline(D.pressTime,'color',[0 0 0]); 
    drawline(D.releaseTime,'color',[0 0 1]);
    drawline(D.endTime,'color',[0 0 1]);
    
    drawline([2.5 D.targetForceAbs D.lowForceAbs D.highForceAbs],'dir','horz'); 
    
    title(['hand ' num2str(D.hand) '; finger ' num2str(D.digit) '; force ' num2str(D.targetForce) '; VAF ' num2str(D.EVal) '; GT ' num2str(D.goodTrial)]);
    subplot(4,1,3); 
    plot(t,state); 
    subplot(4,1,4);
	plot(t(start:end),abs(forceA),'r',t(start:end),sqrt(sum(forceP.^2,2)),'b');
    legend({'Active digit','Passive digits'},'Location','West');
    pause;
%    keyboard;
end;