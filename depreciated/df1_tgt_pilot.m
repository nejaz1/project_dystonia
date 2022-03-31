function varargout=df1_tgt_pilot(what,varargin) 

switch (what)
        
    case 'make_run_pilot'
        vararginoptions(varargin,{'numMiniBlock','pushTime', 'blockNum', 'numTrials', 'sliceNumAnnounce', 'sliceNumTrial', 'numSlicesTR', 'numRestphases', 'dummySlices'})
        numFingers= 5;
        %----define the run conditions
        index= randperm(numFingers*numMiniBlock);   
        digit= repmat(1:numFingers,1,numMiniBlock);
        hand= [ones(1,numFingers*numMiniBlock)];
        stimType= [ones(1,numFingers*numMiniBlock)];
         
        %generate miniblock onset vector
        %----number of slices in a miniBlock
        sliceNumMiniblock= sliceNumAnnounce+ sliceNumTrial*numTrials;
        %----number of slices that are missing to end at a full TR
        restSlices= numSlicesTR- mod(sliceNumMiniblock*(length(digit)+numRestphases-1), numSlicesTR); 
        endSlice= sliceNumMiniblock*(length(digit)+numRestphases-1) + restSlices;      
        %----dummy slices
        %----add them only for the normal seq
        if dummySlices== 1
            dummySlices=4*numSlicesTR;
        end
        %----print number of TR
        fprintf('number of TR: %4.1f \n', (endSlice+dummySlices)/numSlicesTR);
        %----vector of mini blocks onset and add slices to complete the TR
        sliceNumber= [1:sliceNumMiniblock:endSlice, endSlice]+dummySlices;
        %----define restphases
        ind= randperm(numFingers*numMiniBlock+numRestphases); % plus six restphases
        %----delete restphases in the onset vector
        sliceNumber(ind(1:numRestphases))= [];
        %generate runs
        D=[];
        for i=1:length(digit)
            
%             if i==length(digit)
%                 D= make_Trials_end(D, sliceNumber(end));
%             else
                D= make_Trials(D, sliceNumber(i), digit(index(i)), hand(index(i)), stimType(index(i)), pushTime, numTrials, sliceNumAnnounce, sliceNumTrial);
%             end
        end
        %%%%
        D= make_Trials_end(D, sliceNumber(end));
        %%%%
        pivottable(D.digit, [], D.digit, 'size');
        varargout={D};

    case 'make_run'
        vararginoptions(varargin,{'numMiniBlock','pushTime', 'blockNum', 'numTrials', 'sliceNumAnnounce', 'sliceNumTrial', 'numSlicesTR', 'numRestphases', 'dummySlices'})
        numFingers= 5;
        %----define the run conditions
        index= randperm(5*3*numMiniBlock);   
        digit= repmat(1:5,1,3*numMiniBlock);
        hand= [zeros(1,5*numMiniBlock) ones(1,5*2*numMiniBlock)];
        stimType= [zeros(1,5*numMiniBlock) ones(1,5*numMiniBlock) zeros(1,5*numMiniBlock)];
        
        %generate miniblock onset vector
        %----number of slices in a miniBlock
        sliceNumMiniblock= sliceNumAnnounce+ sliceNumTrial*numTrials;
        %----number of slices that are missing to end at a full TR
        restSlices= numSlicesTR- mod(sliceNumMiniblock*(length(digit)+numRestphases), numSlicesTR);
        endSlice= sliceNumMiniblock*(length(digit)+numRestphases) + restSlices;
        %----dummy slices
        %----add them only for the normal seq
        if dummySlices== 1
            dummySlices=4*numSlicesTR;
        end
        %----print number of TR
        fprintf('number of TR: %4.1f \n', (endSlice+dummySlices)/numSlicesTR);
        %----vector of mini blocks onset and add slices to complete the TR
        sliceNumber= [1:sliceNumMiniblock:endSlice, endSlice]+dummySlices;
        %----define restphases
        ind= randperm(numFingers*numMiniBlock+numRestphases); % plus six restphases
        %----delete restphases in the onset vector
        sliceNumber(ind(1:numRestphases))= [];
        %generate runs
        D=[];
        for i=1:length(digit)
            
            %             if i==length(digit)
            %                 D= make_Trials_end(D, sliceNumber(end));
            %             else
            D= make_Trials(D, sliceNumber(i), digit(index(i)), hand(index(i)), stimType(index(i)), pushTime, numTrials, sliceNumAnnounce, sliceNumTrial);
            %             end
        end
        %%%%
        D= make_Trials_end(D, sliceNumber(end));
        %%%%
        pivottable(D.digit, [D.hand D.stimType], D.digit, 'size', 'subset', D.announce==1);
        varargout={D};
        
    case 'make_tgt'
       vararginoptions(varargin,{'subNam'})
       numSlicesTR=32;
       numMiniBlock= 2; 
       numTrials= 6;
       blockNum= 1:10;
       pushTime= 1000; 
       announceTime= 1300; 
       trialTime= 1500; 
       numRestphases= 6;
       %----
       %normal scanner sequence
       sliceNumAnnounce= ceil(announceTime/85);
       sliceNumTrial= ceil(trialTime/85);
       for i=blockNum
            D= df1_tgt_pilot('make_run','numMiniBlock', numMiniBlock,'pushTime', pushTime, 'blockNum', blockNum, 'numTrials', numTrials, ...
                'sliceNumAnnounce', sliceNumAnnounce, 'sliceNumTrial', sliceNumTrial, 'numSlicesTR', numSlicesTR, 'numRestphases', numRestphases, ...
                'dummySlices', 1);  
            dsave ([subNam,'_b',num2str(blockNum(i)), '.tgt'], D);
       end
%        %----
%        %3D scanner sequence
%        numSlicesTR=40;
%        %----
%        sliceNumAnnounce= ceil(announceTime/70);
%        sliceNumTrial= ceil(trialTime/70);
%        for i=blockNum
%             D= df1_tgt_pilot('make_run','numMiniBlock', numMiniBlock,'pushTime', pushTime, 'blockNum', blockNum, 'numTrials', numTrials, ...
%                 'sliceNumAnnounce', sliceNumAnnounce, 'sliceNumTrial', sliceNumTrial, 'numSlicesTR', numSlicesTR, 'numRestphases', numRestphases, ...
%                 'dummySlices', 0); 
%             dsave (['s', num2str(subNum),'_3D_b',num2str(blockNum(i)), '.tgt'], D);
%        end
end; 



        
function C=make_Trials(C, sliceNumber, digit, hand,stimType, pushTime, numTrials, sliceNumAnnounce, sliceNumTrial)
    startSlices= [1,[1:sliceNumTrial:sliceNumTrial*numTrials]+sliceNumAnnounce]-1;
    endSlices= [sliceNumAnnounce,[1:sliceNumTrial:sliceNumTrial*numTrials]+sliceNumAnnounce+sliceNumTrial]-1;
	D.startSlice = startSlices'+ sliceNumber; 
    D.startTime=zeros(numTrials+1,1);
	D.endSlice = endSlices'+ sliceNumber; 
	D.endTime=[1940; ones(numTrials,1)*1500];
	D.digit=repmat(digit,numTrials+1,1); 
    D.hand=repmat(hand,numTrials+1,1); 
    D.stimType=repmat(stimType,numTrials+1,1); 
    D.pushTime= repmat(pushTime,numTrials+1,1); 
	D.announce=[1; zeros(numTrials,1)];
    D.feedback= [0; ones(numTrials,1)];
    
 	C= addstruct(C, D);
    
    function C=make_Trials_end(C, sliceNumber)
	D.startSlice = sliceNumber; 
    D.startTime=0;
	D.endSlice = sliceNumber; 
	D.endTime=10000000;
	D.digit=0; 
    D.hand=0; 
    D.stimType=1; 
    D.pushTime= 0; 
	D.announce=0;
    D.feedback= 0;
    
 	C= addstruct(C, D);   
    
    



