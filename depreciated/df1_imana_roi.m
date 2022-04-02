function varargout = df1_imana_roi(what,varargin)

% Sadnicka to be added to larger df1_imana once finalised

rootDir     = '/Volumes/ANNA/data/FingerPattern_dystonia/';
roiDir      = [rootDir 'RegionOfInterest'];

hemName         = {'LeftHem','RightHem'};
hemNum          = {1,1,1,1,1,2,2,2,2,2};
regionName      = {'S1','M1','PMd','Cant','Cinf','S1','M1','PMd','Cant','Cinf'};
regionNum       = {1,2,3,9,10,11,12,13,19,20};

cd(roiDir)

% Analysis Functions

switch (what)
     
 %%% CONTRALATERAL DATA
    case 'LeftHem_RightMotor' % ttest(psc), ttest(mean distance), rmANOVA(ind distance)
            D = load('reg_distance_raw.mat');
            for i = 1:5; % each ROI (left hemisphere)
                r = getrow(D,ismember(D.region,[regionNum{i}])); % select data from ROI
                r = getrow(r, r.hand~=r.regSide & r.stimtype==0); % hand/regSide not equal.  stimtype motor.         

                fprintf(hemName{hemNum{i}}), fprintf(regionName{i}) % label command window 

                fprintf('\nPSC\n')
                [t,p] = ttest(r.psc(r.group==1), r.psc(r.group==2),2,'independent') % t-test psc
                subplot(3,5,i), barplot(r.group, r.psc); title(sprintf(regionName{i}))  %optional plot

                fprintf('\nmean distance\n')
                d = mean(r.dist')'; % calculate mean distance
                [t,p]= ttest(d(r.group==1), d(r.group==2), 2,'independent')
                subplot(3,5,i+5), barplot(r.group, d); title(sprintf(regionName{i}))   %optional plot

                fprintf('\ndistance rmANOVA\n')
                subplot(3,5,i+10), traceplot(1:10,r.dist,'split',r.group,'errorfcn','stderr','leg','auto'), title(sprintf(regionName{i}))  
                MANOVA1(r.group,r.dist)
            end

    case 'RightHem_LeftMotor' % ttest(psc), ttest(mean distance), rmANOVA(ind distance)
            D = load('reg_distance_raw.mat');
            for i = 6:10; % each ROI (right hemisphere)  
                r = getrow(D,ismember(D.region,[regionNum{i}])); 
                r = getrow(r, r.hand~=r.regSide & r.stimtype==0); % hand/regSide not equal.  stimtype motor.           

                fprintf(hemName{hemNum{i}}), fprintf(regionName{i}) % label command window 

                fprintf('\nPSC\n')
                [t,p] = ttest(r.psc(r.group==1), r.psc(r.group==2),2,'independent') % t-test psc
                subplot(3,5,i-5), barplot(r.group, r.psc); title(sprintf(regionName{i}))  %optional plot

                fprintf('\nmean distance\n')
                d = mean(r.dist')';
                [t,p]= ttest(d(r.group==1), d(r.group==2), 2,'independent')
                subplot(3,5,i), barplot(r.group, d); title(sprintf(regionName{i}))   %optional plot

                fprintf('\ndistance rmANOVA\n')
                subplot(3,5,i+5), traceplot(1:10,r.dist,'split',r.group,'errorfcn','stderr','leg','auto'), title(sprintf(regionName{i}))  
                MANOVA1(r.group,r.dist)
            end

    case 'LeftHem_RightSensory'% ttest(psc), ttest(mean distance), rmANOVA(ind distance)
            D = load('reg_distance_raw.mat');
            for i = 1:5; % each ROI (left hemisphere),
                r = getrow(D,ismember(D.region,[regionNum{i}])); 
                r = getrow(r, r.hand~=r.regSide & r.stimtype==1); % hand/regSide not equal.  stimtype sensory.           

                fprintf(hemName{hemNum{i}}), fprintf(regionName{i}) % label command window 

                fprintf('\nPSC\n')
                [t,p] = ttest(r.psc(r.group==1), r.psc(r.group==2),2,'independent') % t-test psc
                subplot(3,5,i), myboxplot(r.group, r.psc); title(sprintf(regionName{i}))  %optional plot

                fprintf('\nmean distance\n')
                d = mean(r.dist')';
                [t,p]= ttest(d(r.group==1), d(r.group==2), 2,'independent')
                subplot(3,5,i+5), myboxplot(r.group, d); title(sprintf(regionName{i}))   %optional plot

                fprintf('\ndistance rmANOVA\n')
                subplot(3,5,i+10), traceplot(1:10,r.dist,'split',r.group,'errorfcn','stderr','leg','auto'), title(sprintf(regionName{i}))  
                MANOVA1(r.group,r.dist)
            end

        %%% IPSILATERAL DATA
    case 'LeftHem_LeftMotor' % ttest(psc), ttest(mean distance), rmANOVA(ind distance)
            D = load('reg_distance_raw.mat');
            for i = 1:5; % each ROI, 
                r = getrow(D,ismember(D.region,[regionNum{i}])); 
                r = getrow(r, r.hand==r.regSide & r.stimtype==0); 

                fprintf(hemName{hemNum{i}}), fprintf(regionName{i}) 

                fprintf('\nPSC\n')
                [t,p] = ttest(r.psc(r.group==1), r.psc(r.group==2),2,'independent') % t-test psc
                subplot(3,5,i), myboxplot(r.group, r.psc); title(sprintf(regionName{i}))  %optional plot

                fprintf('\nmean distance\n')
                d = mean(r.dist')';
                [t,p]= ttest(d(r.group==1), d(r.group==2), 2,'independent')

                fprintf('\ndistance rmANOVA\n')
                MANOVA1(r.group,r.dist)

                subplot(3,5,i),     myboxplot(r.group, r.psc); title(sprintf(regionName{i}))  %optional plot
                subplot(3,5,i+5),   myboxplot(r.group, d); title(sprintf(regionName{i}))   %optional plot
                subplot(3,5,i+10),  traceplot(1:10,r.dist,'split',r.group,'errorfcn','stderr','leg','auto'); title(sprintf(regionName{i}))  
            end

    case 'RightHem_RightMotor'% ttest(psc), ttest(mean distance), rmANOVA(ind distance)
            D = load('reg_distance_raw.mat');
            for i = 6:10; % each ROI,  
                r = getrow(D,ismember(D.region,[regionNum{i}])); 
                r = getrow(r, r.hand==r.regSide & r.stimtype==0); 
                fprintf(hemName{hemNum{i}}), fprintf(regionName{i}) 

                fprintf('\nPSC\n')
                [t,p] = ttest(r.psc(r.group==1), r.psc(r.group==2),2,'independent') % t-test psc

                fprintf('\nmean distance\n')
                d = mean(r.dist')';
                [t,p] = ttest(d(r.group==1), d(r.group==2), 2,'independent')

                fprintf('\ndistance rmANOVA\n')
                MANOVA1(r.group,r.dist) 

                subplot(3,5,i-5),   barplot(r.group, r.psc); title(sprintf(regionName{i}))  %optional plot
                subplot(3,5,i),     barplot(r.group, d); title(sprintf(regionName{i}))   
                subplot(3,5,i+5),   traceplot(1:10,r.dist,'split',r.group,'errorfcn','stderr','leg','auto'); title(sprintf(regionName{i}))  
         end

    case 'RightHem_RightSensory' % ttest(psc), ttest(mean distance), rmANOVA(ind distance)
            D = load('reg_distance_raw.mat');
            for i = 6:10; % each ROI,  
                r = getrow(D,ismember(D.region,[regionNum{i}])); 
                r = getrow(r, r.hand==r.regSide & r.stimtype==1); % sensory

                fprintf(hemName{hemNum{i}}), fprintf(regionName{i}) 

                fprintf('\nPSC\n')
                [t,p] = ttest(r.psc(r.group==1), r.psc(r.group==2),2,'independent') % t-test psc

                fprintf('\nmean distance\n')
                d = mean(r.dist')';
                [t,p] = ttest(d(r.group==1), d(r.group==2), 2,'independent')

                fprintf('\ndistance rmANOVA\n')
                MANOVA1(r.group,r.dist)

                subplot(3,5,i-5),   barplot(r.group, r.psc); title(sprintf(regionName{i}))  %optional plot
                subplot(3,5,i),     barplot(r.group, d); title(sprintf(regionName{i}))   
                subplot(3,5,i+5),   traceplot(1:10,r.dist,'split',r.group,'errorfcn','stderr','leg','auto'); title(sprintf(regionName{i}));  
            end

    case 'REL_betweenSubj'  % corr distance measures between subjects
            region      = [2 12];
            stimType    = 0;
            group       = 1;
            vararginoptions(varargin,{'region','stimType','group'})


            % 1. Get subset of data 
            D = load(fullfile(roiDir,'reg_distance_raw.mat'));
            D = getrow(D,ismember(D.region,region) & D.stimtype==stimType & ...
                         D.group==group & D.hand~=D.regSide);

            % 2. Get one distance structure for each subject
            T = tapply(D,{'SN'},{'dist','mean(x,1)'});

            % 3. Get correlation between subjects
            i = 1;
            r = zeros(length(T.SN),1);
            for s=unique(T.SN)'
                d1      = T.dist(T.SN==s,:);
                d2      = T.dist(T.SN~=s,:);
                r(i)    = mean(corr(d1',d2'));
                i       = i+1;
            end;

            % 4. Display intersubject reliabilities
            fprintf('Mean: %2.3f, SE: %2.3f\n',mean(r),stderr(r))
            varargout = {r};

    case 'REL_betweenSubj_ttest'

        
            for i = 1:5 % motor condition (combining both left and right hemispheres)

                r=unique(regionNum);
                region = [r(i), r(i+5)];
                g1 = roi_dystonia('REL_betweenSubj', 'region', region, 'group', 1);
                g2 = roi_dystonia('REL_betweenSubj', 'region', region, 'group', 2);
                fprintf(regionName{i})
                ttest(g1,g2,2,'independent')
            end

            for i = 1:5 % repeat for sensory condition (left hemisphere only)
                r=unique(regionNum);
                region = [r(i)];
                g1 = roi_dystonia('REL_betweenSubj', 'region', region, 'group', 1, 'stimType', 1);
                g2 = roi_dystonia('REL_betweenSubj', 'region', region, 'group', 2, 'stimType', 1);
                fprintf(regionName{i})
                ttest(g1,g2,2,'independent')
            end

    case 'REL_intraHem' % 
            dim = 5;
            region      = [2 12];
            stimType    = 0;
            group       = 1;
            vararginoptions(varargin,{'region','stimType','group'})

            % 1. Get subset of data 
            D = load(fullfile(roiDir,'reg_distance_raw.mat'));
            D = getrow(D,ismember(D.region,region) & D.stimtype==stimType & ...
                         D.group==group & D.hand~=D.regSide);

            % 2. MANOVA fixed factor hemisphere (rmANOVA not enough data)
            [u,s,v] = svd(D.dist);
            D.distRed = u(:,1:dim)*s(1:dim,1:dim);
    %         MANOVArp(D.regSide, D.SN, D.distRed);
            MANOVA1(D.regSide, D.distRed);
    
    case 'ROI_psc_meanDistance_plot'                  % Bar plots of psc & mean distance values
        % df1_imana_roi('ROI_psc_meanDistance_plot',2); select ROI
        
        D=load(fullfile(roiDir,'reg_distance_raw.mat'));
        regType=varargin{1};
        HAND=[0 1 1];
        STIMTYPE=[0 0 1];
        plotname={'Left Motor','Right Motor','Right Sensory'};
        facecolor={[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1],[1 1 1]};
        edgecolor={[0 0 1],[1 0 0],[0 0 1],[1 0 0],[0 0 1],[1 0 0]};
        for i = 1:2;
            subplot(4,1,i*2-1)
            d=getrow(D,(D.regType==regType & D.regSide==(i-1)));
            barplot([d.hand d.stimtype],d.psc,'split', [d.group], 'facecolor', facecolor, 'edgecolor', edgecolor, 'gapwidth', [0.1 0.1 0.1]);
            ylabel({(regionName{regType}); (hemName{i})})
            set(gca, 'XTickLabel', {'LM' , 'LM', 'RM', 'RM', 'RS', 'RS'}, 'fontSize', 8)
            title('PSC')
        end
        
        for i = 1:2;
            subplot(4,1,i*2) 
            d=getrow(D,(D.regType==regType & D.regSide==(i-1)));
            barplot([d.hand d.stimtype], mean(d.dist,2),'split', [d.group], 'facecolor', facecolor, 'edgecolor', edgecolor, 'gapwidth', [0.1 0.1 0.1]);
            ylabel({(regionName{regType}); (hemName{i})})
            set(gca, 'XTickLabel', {'LM' , 'LM', 'RM', 'RM', 'RS', 'RS'}, 'fontSize', 8)
            title('Mean Distance')
        end

        
    
    case 'ROI_distance_plot'                  % Plot distance values
        % df1_imana_roi('ROI_distance_plot',1,'dist');
        
        color={[0 0 1],[1 0 0],[0 0 1],[1 0 0]};
        style={'-','-',':',':'};
        D=load(fullfile(roiDir,'reg_distance_raw.mat'));
        regType=varargin{1};
        type = varargin{2}; 
        HAND=[0 1 1];
        STIMTYPE=[0 0 1];
        plotname={'Left Motor','Right Motor','Right Sensory'};
        xticklabels={'1/2', '1/3', '1/4', '1/5', '2/3', '2/4', '2/5', '3/4', '3/5', '4/5'};
        % D.ndist=bsxfun(@rdivide,D.dist,sqrt(sum(D.dist.^2,2)));
        
        set(gcf,'Name',regionName{regType})
        set(gca, 'fontsize', 16);

        D=getrow(D,D.regType==regType);
        
        for h=1:3
            for s=0:1
                
                subplot(2,3,h+s*3);
                traceplot(1:10,D.(type),'split',[D.group], ...
                    'subset',D.regType==regType & D.hand==HAND(h) & D.stimtype==STIMTYPE(h) & D.regSide==s,...
                    'linecolor',color,'patchcolor',color,'linestyle',style,'errorfcn','stderr');
               
                d = getrow(D, D.regType==regType & D.hand==HAND(h) & D.stimtype==STIMTYPE(h) & D.regSide==s);
                maxim=max(mean(d.(type))+8*stderr(d.(type)));
                set(gca,'YLim',[-maxim/10 maxim*0.8], 'fontsize', 14, ...
                        'XTick', [1:10], 'XTickLabel', {'1/2', '1/3', '1/4', '1/5', '2/3', '2/4', '2/5', '3/4', '3/5', '4/5'});
               
                if(s==0)
                    title(plotname{h});
                end;
                if(h==1 & s==0)
                    ylabel({(regionName{regType}); 'Left Hemisphere'}, 'FontWeight', 'bold');
                end;
                if(h==1 & s==1)
                    ylabel({(regionName{regType}); 'Right Hemisphere'}, 'FontWeight', 'bold');
                end;
               
            end;
        end;  
end

