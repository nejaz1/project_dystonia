function D=Individjhu2b_trial_dystonia(MOV,D,fig,varargin)
% Analysis of a subject from Individuation2b - bedside version task 

vararginoptions(varargin,{'startthres','endthres','minstart','maxstop','forceThreshold'});


% 0. make variables consistent with older experiments which have force instead of target force as a field
if isfield(D,'force')
    tForce = getfield(D,'force');
    D = rmfield(D,'force');
    D = setfield(D,'targetForce',tForce);
    D = setfield(D,'lowForce',tForce.*0.75);
    D = setfield(D,'highForce',tForce.*1.25);
    clear tForce;
end;


% 1. Extract data:
if (isempty(MOV))
    return;
end;
state   = MOV(:,1);
t       = MOV(:,2);
Fs      = mean(diff(t));
screen_t= MOV(:,3);                 % Screentime 
Force   = MOV(:,4:13);              % cols 4-13 are the Force data
Force   = smooth_kernel(Force,5);   % Smoothing with Gaussian kernel 
mvc     = max(D.mvc,1);             % if mvc for any finger is less than one, it is set to one


% 2. Indices for important experimental variables
% - state timing events
% - active/passive hands
I2          = find(state==2,1,'first');
I5          = find(state==5,1,'first'); 
if isempty(I5) I5 = length(state); end;

ts          = t(I2:I5);
D.dt        = mean(diff(ts));
IndActHand  = (D.hand-1)*5+[1:5]; 
IndPassHand = (2-D.hand)*5+[1:5]; 
IndActDigit = (D.hand-1)*5+D.digit; % instructed finger
IndHomDigit = (2-D.hand)*5+D.digit;
PassFin     = setdiff(IndActHand,IndActDigit);  % indices of four passive fingers in the instructed hand

% 3. Metrics used for determine trial exclusion criteria
% - total length of the trial
% - ration between rms active vs rms passive finger movement
% - variance account for by first principle component of movement

D.baselineForce = Force(I2,:); 
normF           = bsxfun(@minus,Force(I2:I5,:),D.baselineForce);
relnormF        = bsxfun(@rdivide,normF,mvc);

S           = normF'*normF;
[EVec,EVal] = eig(S);
EVal        = diag(EVal);
[EVal,sIndx]= sort(EVal,'descend');
EVec        = EVec(sIndx);
D.EVal      = EVal(1)/sum(EVal);
D.EVec      = EVec(:,1)';


% 4. Storing forces and covariance amongst digits
dF      = velocity_discr(normF,3);    
CovF    = cov(dF);       
D.CovF  = CovF(:)';    

D.maxforce    = max(normF);              % saving peak forces on all 10 digits
D.actHForce   = D.maxforce(:,IndActHand);  % max force on fingers in the instructed hand
D.relmaxforce = D.maxforce./mvc;         % force scaled by the mvc for every finger
D.passHForce  = D.maxforce(:,IndPassHand); % max force on fingers in the passive hand
D.actForce    = max(normF(:,IndActDigit));

% 5. Deviation from optimal trajectory
% - unnormalized force
% - force normalized by mvc
forceA          = normF(:,IndActDigit);     
forceP          = normF(:,PassFin);
D.actForce      = max(forceA);
D.APDiff        = max(forceA)-max(max(forceP),[],2);
D.meanDevP      = mean(sqrt(sum(forceP.^2,2)));     % Mean Eucledian deviation from straight 
D.meanDevA      = mean(sqrt(sum(forceA.^2,2)));     % Mean movement along straight line 
D.meanDevDiff   = D.meanDevA-D.meanDevP;

relforceA       = relnormF(:,IndActDigit); 
relforceP       = relnormF(:,PassFin);
D.relmeanDevP   = mean(sqrt(sum(relforceP.^2,2)));  % Mean Eucledian deviation from straight 
D.relmeanDevA   = mean(sqrt(sum(relforceA.^2,2)));  % Mean movement along straight line 
D.relmeanDevDiff= D.relmeanDevA-D.relmeanDevP;


% 7. Calculating degree of mirrored movements
% - mean deviation over all the passive fingers of the other hand
% - RMS force on all fingers of either hand
D.meanDevPMM    = mean(sqrt(sum(normF(:,IndPassHand).^2,2)));
D.MeanDevPDigit = sqrt(nanmean(normF.^2));

% 8. lag/corr between moving and non-moving hand
% fPassHand       = sqrt(mean(normF(:,IndPassHand).^2,2));
% fActHand        = sqrt(mean(normF(:,IndActHand).^2,2));
% [a,l]           = xcov(fPassHand,fActHand);
% [~,i]           = max(a);
% D.lag           = l(i) * Fs;
% D.rlag          = a(i)/(length(fPassHand)*std(fPassHand)*std(fActHand));
fPassHand       = normF(:,IndActDigit);
fActHand        = sqrt(mean(normF(:,PassFin).^2,2));
D.rlag          = corr(fPassHand,fActHand);
[a,l]           = xcov(fPassHand,fActHand,[],'unbiased');
[~,i]           = max(a);
D.lag           = l(i) * D.dt;


% 9. Peak forces on moving and non-moving hands
D.actPForce     = max(abs(normF(:,IndActDigit)));
D.ensPForce     = max(max(abs(normF(:,PassFin))));
D.mmPForce      = max(max(abs(normF(:,IndPassHand))));

D.peakA         = max(abs(normF(:,IndActDigit)));
D.peakENS       = max(mean(abs(normF(:,PassFin)),2));
D.peakMM        = max(mean(abs(normF(:,IndPassHand)),2));
D.peakF         = max(abs(normF),[],1);


% -------------------------------------------------------------------------
% Display trial

if (fig>0)
    h=figure(1);
    subplot(3,2,[1 3]); 
    plot(t,Force(:,1),'k-',...
        t,Force(:,2),'b-',...
        t,Force(:,3),'g-',...
        t,Force(:,4),'r-',...
        t,Force(:,5),'y-',...
        t,Force(:,6),'k--',...
        t,Force(:,7),'b--',...
        t,Force(:,8),'g--',...
        t,Force(:,9),'r--',...
        t,Force(:,10),'y--','LineWidth',2);
    leg = {'l1','l2','l3','l4','l5','r1','r2','r3','r4','r5'};
    legend(leg,'Location','NorthWest');
    drawline(t(I2),'dir','vert','color','r');
    drawline(t(I5),'dir','vert','color','r');

    title(['hand ' num2str(D.hand) '; finger ' num2str(D.digit) '; force '...
        num2str(D.targetForce) '; VAF ' num2str(D.EVal)]);

    subplot(3,2,5); 
    plot(t,state,'k');
    legend('States','Location','NorthWest');
    

    subplot(3,2,[2 4]);
    plot(t,Force(:,IndPassHand(1)),'k--',...
        t,Force(:,IndPassHand(2)),'b--',...
        t,Force(:,IndPassHand(3)),'g--',...
        t,Force(:,IndPassHand(4)),'r--',...
        t,Force(:,IndPassHand(5)),'y--','LineWidth',2);
    title('Passive hand');
    drawline(t(I2),'dir','vert','color','r');
    drawline(t(I5),'dir','vert','color','r');
    
    keyboard;
end;