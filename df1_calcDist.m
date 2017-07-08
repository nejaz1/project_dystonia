function  [dist] = df1_calcDist(Y,SPM,D)% data,condition,sequence,block
% Caluclate Accuracies (Nearest Neighbour) and Squared Eucledian distamces
% based on on the prewhitened Beta weigths from the raw data
beta=mva_prewhiten_beta(Y',SPM);
C=indicatorMatrix('allpairs',[1:5]);

HAND=[0 1 1]; 
ST=[0 0 1]; 
dist=nan(3,1); 
 
for i=1:length(HAND)
    indx    = (D.stimType==ST(i) &  D.hand==HAND(i)); 
    X       = indicatorMatrix('identity_p',D.digit.*indx);
    c       = D.digit(indx); 
    d       = distance_ldc(beta,X,C,D.run);
    dist(i) = mean(d); 
end;
