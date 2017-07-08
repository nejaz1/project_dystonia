function varargout=df1_tgt(what,varargin) 

switch (what)
        
    case 'make_run'
        vararginoptions(varargin,{'numMiniBlock','pushTime', 'blockNum', 'numberTrials'})
        index= randperm(5*3*numMiniBlock);   
        digit= repmat(1:5,1,3*numMiniBlock);
        hand= [zeros(1,5*numMiniBlock) ones(1,5*2*numMiniBlock)];
        stimType= [zeros(1,5*numMiniBlock) ones(1,5*numMiniBlock) zeros(1,5*numMiniBlock)];
        % (23+ 18*numberTrials)*5*3*numMiniBlock; 
        % 6688 % gives us 209TR
        sliceNumber= 1:131:6688;
        ind= randperm(5*3*numMiniBlock+6); % plus six restphases
        sliceNumber(ind(1:6))= [];
        D=[];
        for i=1:length(digit)
            D= make_Trials(D, sliceNumber(i), digit(index(i)), hand(index(i)), stimType(index(i)), pushTime, numberTrials);
        end
        varargout={D};
        
    case 'make_tgt'
       vararginoptions(varargin,{'subNum'})
       numMiniBlock= 3; 
       numberTrials= 6;
       blockNum= 1:10;
       pushTime= 1000; 
       for i=blockNum
            D= df1_tgt('make_run','numMiniBlock', numMiniBlock,'pushTime', pushTime, 'blockNum', blockNum, 'numberTrials', numberTrials); 
            dsave (['s', num2str(subNum),'_b',num2str(blockNum(i)), '.tgt'], D);
       end
end; 



        
function C=make_Trials(C, sliceNumber, digit, hand,stimType, pushTime, numberTrials)
    startSlices= [1,[1:18:18*numberTrials]+23]-1;
    endSlices= [23,[1:18:18*numberTrials]+23+18]-1;
	D.startSlice = startSlices'+ sliceNumber; 
    D.startTime=zeros(numberTrials+1,1);
	D.endSlice = endSlices'+ sliceNumber; 
	D.endTime=[1940; ones(numberTrials,1)*1500];
	D.digit=repmat(digit,numberTrials+1,1); 
    D.hand=repmat(hand,numberTrials+1,1); 
    D.stimType=repmat(stimType,numberTrials+1,1); 
    D.pushTime= repmat(pushTime,numberTrials+1,1); 
	D.announce=[1; zeros(numberTrials,1)];
    D.feedback= [0; ones(numberTrials,1)];
    
 	C= addstruct(C, D);
       
    
    



