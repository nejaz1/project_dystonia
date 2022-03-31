function varargout=df1_IndividAnalyze_anita(what,varargin)

% Analysis of individuation data. 
% Anna Sadnicka May 2014

% baseDir=            '/Users/asadnicka/Documents/data/FingerPattern_dystonia'; % Sadnicks laptop
% baseDir=          'C:/Projects/FingerPattern_dystonia'; % Sadnicks Sobell desktop
baseDir=          '/Users/naveed/Documents/data/FingerPattern_dystonia';
% baseDir =         '/Users/joern/Projects/fingerPattern_dystonia';

behaviourEMGDir =   fullfile(baseDir, 'Individuation_EMG');

% Load text file to generate labels
sub = dload(fullfile(behaviourEMGDir,'subject_list.txt'));

switch (what)
    
    case 'subj' %loops through analysis for available subjects defined in subject_list.txt
        cd (fullfile(behaviourEMGDir,'analysis')); 
        for i=1:length(sub.name); 
            df1_subj(sub.name{i},0);
        end; 
             
    case 'make_alldat'     
        S=[];
        cd(fullfile(behaviourEMGDir,'analysis')); 
        for i=1:length(sub.name)
            Si = load(fullfile(sprintf('IN2b_%s_ana.mat',sub.name{i})));
            Si.SN =         zeros(length(Si.BN),1)+i;
            Si.group =      zeros(length(Si.BN),1)+sub.group(i);
            Si.instrum =    zeros(length(Si.BN),1)+sub.instrum(i);
            Si.scan =       zeros(length(Si.BN),1)+sub.scan(i);
            Si.L1 =         zeros(length(Si.BN),1)+sub.L1(i);
            Si.L2 =         zeros(length(Si.BN),1)+sub.L2(i);
            Si.L3 =         zeros(length(Si.BN),1)+sub.L3(i);
            Si.L4 =         zeros(length(Si.BN),1)+sub.L4(i);
            Si.L5 =         zeros(length(Si.BN),1)+sub.L5(i);
            Si.R1 =         zeros(length(Si.BN),1)+sub.R1(i);
            Si.R2 =         zeros(length(Si.BN),1)+sub.R2(i);
            Si.R3 =         zeros(length(Si.BN),1)+sub.R3(i);
            Si.R4 =         zeros(length(Si.BN),1)+sub.R4(i);
            Si.R5 =         zeros(length(Si.BN),1)+sub.R5(i);
            S=addstruct(S,Si);
        end;
        save(fullfile(behaviourEMGDir,'analysis','IN2b_alldat.mat'),'-struct','S');
    
    case 'devF' % Generates matrix & image of devF variable for each subject.  
        cd(fullfile(baseDir,'Individuation_EMG','analysis'));
        D = load('IN2b_alldat.mat');
        
        dev_L = zeros(5,5);
        dev_R = zeros(5,5);
        for i=unique(D.SN)' % need help writing loop with this ... do I need to use eval so that I can change the ...
            % variable name on each loop?
            dL = getrow(D,D.SN==i & D.hand==1); % data for subj i and hand 1
            dR = getrow(D,D.SN==i & D.hand==2); % data for subj i and hand 2
            
            dev_L = pivottablerow(dL.digit,dL.meanDevF,'mean'); % get meanDevF for each finger
            dev_R = pivottablerow(dR.digit,dR.meanDevF,'mean');
            
            dev_L = eye(5) + triu(dev_L,1) + tril(dev_L,-1);    % removing main diagonal
            dev_R = eye(5) + triu(dev_R,1) + tril(dev_R,-1);    

%             %% Saves each matrix. Subject number given at end.
%             filename = sprintf('devL_%d.mat', i);
%             save(filename, 'dev_L')
%             
%             filename = sprintf('devR_%d.mat', i);
%             save(filename, 'dev_R')
%             
%             %% Graphically demonstrates devF for each individual
%             %% subject and saves a copy
%             
%             devF_fig = figure;
%             hold on
%             subplot(1,2,1)
%             imagesc_rectangle(dev_L, 'YDir', 'reverse', 'scale', [0 1]);
%             colorbar
%             title('LH')
%             subplot(1,2,2)
%             imagesc_rectangle(dev_R, 'YDir', 'reverse', 'scale', [0 1]);
%             colorbar
%             title('RH')
%             saveas (devF_fig, sprintf('sn%d.bmp',i))
            
            %% Calculate Individuation measures
 
            
            % replace diagonal of '1' with NaN
            dev_L = diag(nan(1,5)) + dev_L;
            dev_R = diag(nan(1,5)) + dev_R;

            % indA = Individuation A = row mean (i.e. when finger x moves how much did other fingers move
            indA_left = nanmean(dev_L');
            indA_right = nanmean(dev_R');
            
            % IndB = Individuation B = column mean (i.e. when other four fingers move what was the mean amount of activation for finger x)
            indB_left= nanmean(dev_L);
            indB_right= nanmean(dev_R);
            
            
            %% Add individuation measures to subject text file
            sub.indA_L1(i,1) = indA_left(1,1);
            sub.indA_L2(i,1) = indA_left(1,2);
            sub.indA_L3(i,1) = indA_left(1,3);
            sub.indA_L4(i,1) = indA_left(1,4);
            sub.indA_L5(i,1) = indA_left(1,5);
            sub.indA_R1(i,1) = indA_right(1,1);
            sub.indA_R2(i,1) = indA_right(1,2);
            sub.indA_R3(i,1) = indA_right(1,3);
            sub.indA_R4(i,1) = indA_right(1,4);
            sub.indA_R5(i,1) = indA_right(1,5);
          
            sub.indB_L1(i,1) = indB_left(1,1);
            sub.indB_L2(i,1) = indB_left(1,2);
            sub.indB_L3(i,1) = indB_left(1,3);
            sub.indB_L4(i,1) = indB_left(1,4);
            sub.indB_L5(i,1) = indB_left(1,5);
            sub.indB_R1(i,1) = indB_right(1,1);
            sub.indB_R2(i,1) = indB_right(1,2);
            sub.indB_R3(i,1) = indB_right(1,3);
            sub.indB_R4(i,1) = indB_right(1,4);
            sub.indB_R5(i,1) = indB_right(1,5);
        end
        dsave('subject_devA&B.txt', sub)

    case 'ind_summ' % Compares individuation measures 
            
            %% (A) dystonia vs control
            D = importdata('subject_devA&B.txt');
            
            % dystonia vs control index 
            dys_i = find(D.data(:,1) ~=1);
            cont_i = find(D.data(:,1) ==1);
            
            %%  Matrix with clinically affected fingers and graph with hot shading
            clinAffect = D.data(dys_i, 4:13);
            figure
            imagesc(clinAffect)
            axis([0.5 10.5 0.5 11.5])
            title ('Dystonia patients - symptomatic digits red')
            set(gca, 'XTickLabel', str2mat('1_left','2_left','3_left','4_left','5_left','1_right','2_right','3_right','4_right','5_right'))
            xlabel('digit (1=thumb:5=little finger)')
            set(gca, 'YTickLabel', str2mat('piano','piano','piano','guitar','guitar','guitar','piano','guitar','guitar','guitar','piano'))
            ylabel('dystonic subject')
            keyboard;
            
            %% LEFT HAND - pull data, perform t-test, plot bar
           
            % all digits
            control (1,1) = mean2(D.data(cont_i, 14:23));
            dystonia (1,1) = mean2(D.data(dys_i, 14:23));
            [t_allA, p_allA] = ttest((D.data(cont_i, 14:23)), (D.data(dys_i, 14:23)), 2, 'independent');
            [t_allB, p_allB] = ttest((D.data(cont_i, 24:33)), (D.data(dys_i, 24:33)), 2, 'independent');
            % left hand
            control (1,2) = mean2(D.data(cont_i, 14:18));
            dystonia (1,2) = mean2(D.data(dys_i, 14:18));
            [t_leftA, p_leftA] = ttest((D.data(cont_i, 14:18)), (D.data(dys_i, 14:18)), 2, 'independent');
            [t_leftB, p_leftB] = ttest((D.data(cont_i, 24:28)), (D.data(dys_i, 24:28)), 2, 'independent');
            % right hand
            control (1,3) = mean2(D.data(cont_i, 19:23));
            dystonia (1,3) = mean2(D.data(dys_i, 19:23));
            [t_rightA, p_rightA] = ttest((D.data(cont_i, 19:23)), (D.data(dys_i, 19:23)), 2, 'independent');
            [t_rightB, p_rightB] = ttest((D.data(cont_i, 29:33)), (D.data(dys_i, 29:33)), 2, 'independent');
            
            figure(2)% difference of mean individuation measures by hand and all digits
            bar(1:3, [control' dystonia'], 1)
            title({'Comparison of Individuation Measures'; 'All digits: IndA p = 0.0048 IndB p = 0.012'; 'Left hand: IndA p = 0.14 IndB p = 0.17'; 'Right Hand: IndA p = 0.015 IndB p = 0.033'})
            set(gca, 'XTickLabel', str2mat('all digits','left hand','right hand'));
            ylabel('Mean Individuation Measure (IndA==IndB)')
            legend('control', 'dystonia')
            
            figure(3) %IndA by digit
            subplot(2,1,1)
            control_digit = mean(D.data(cont_i, 14:23));
            dystonia_digit = mean(D.data(dys_i, 14:23));    
            bar(1:10, [control_digit' dystonia_digit'], 1)
            title('Individuation A')
            legend('control', 'dystonia')
            set(gca, 'XTickLabel', str2mat('1_left','2_left','3_left','4_left','5_left','1_right','2_right','3_right','4_right','5_right'))
            
            subplot(2,1,2)
            control_digit = mean(D.data(cont_i, 24:33));
            dystonia_digit = mean(D.data(dys_i, 24:33));    
            bar(1:10, [control_digit' dystonia_digit'], 1)
            title('Individuation B')
            legend('control', 'dystonia')
            set(gca, 'XTickLabel', str2mat('1_left','2_left','3_left','4_left','5_left','1_right','2_right','3_right','4_right','5_right'))
            
            %% (B) within dystonia group symptomatic fingers vs non symptomatic fingers
            
            % Symptomatic fingers
            %  x y co-ordinates for symptomatic fingers
            [x y] = find(clinAffect ==1);
            
            % pulls indA (s_A) and indB (s_B) values for symptomatic fingers
            symp_indA = [];
            for i = 1:length(x)
                ai = dys_indA(x(i,1), y(i,1));
                symp_indA = [symp_indA,ai];
            end  
            
            symp_indB = [];
            for i = 1:length(x)
                bi = dys_indB(x(i,1), y(i,1));
                symp_indB = [symp_indB,bi];
            end  
 
            %% Assymptomatic fingers
            % gets x y co-ordinates for unaffected fingers
            [x y] = find(clinAffect ==0);  
            
            % pulls indA (a_A)and indB (a_B) values for asymptomatic fingers
            asymp_indA = [];
            for i = 1:length(x)
                ai = dys_indA(x(i,1), y(i,1));
                asymp_indA = [asymp_indA,ai];
            end  
            
            asymp_indB = [];
            for i = 1:length(x)
                bi = dys_indB(x(i,1), y(i,1));
                asymp_indB = [asymp_indB,bi];
            end  
           
            
            %% T-TEST
            % Explore difference between symptomatic and
            % asymptomatic fingers for both IndA and IndB
            
           ttest(symp_indA, asymp_indA, 2, 'independent');
           ttest(symp_indB, asymp_indB, 2, 'independent');
        
    case 'plot_meanDevF' % For each subject this plots actForce vs meanDevP for individual fingers...
%        % showing LH and RH as separate plots'
%        
        D = [];
        cd(fullfile(baseDir,'Individuation_EMG','analysis'));
        D = load('IN2b_alldat.mat');
        
        for i=1:unique(D.SN);
           
        data = pivottable([D.hand D.digit], D.targetForce, D.meanDevF, 'mean', 'subset', D.SN == i);
        result = [];
        result(((i-1))*10+i:i*10, 1) = i; % labels subject number
        
            subplot(2,1,1)

            for j=1:5;

                %CAT.linecolor = {'r', 'b', 'g', 'k', 'm'};
                %CAT.markercolor = {'r', 'b', 'g', 'k', 'm'};
                %CAT.markerfill = {'r', 'b', 'g', 'k', 'm'};

                result((i-1)*10+1:(i-1)*10+5, 2) = 1; % labels hand

                hold on
                x = 25:25:75;
                y = data(j,1:3);
                [r2, b, t, p] = scatterplot(x',y', 'regression', 'linear', 'intercept', 0);  
                result((i-1)*10+j, 3) = r2; 
                result((i-1)*10+j, 4) = b;
                result((i-1)*10+j, 5) = t;
                result((i-1)*10+j, 6) = p;

             end

            % legend('1', '2', '3', '4', '5')
            title('Left Hand')
            axis([0 75 -5 30])
            
            
            subplot(2,1,2)

            for j=1:5;

                %CAT.linecolor = {'r', 'b', 'g', 'k', 'm'};
                %CAT.markercolor = {'r', 'b', 'g', 'k', 'm'};
                %CAT.markerfill = {'r', 'b', 'g', 'k', 'm'};

                result((i-1)*10+6:(i-1)*10+10, 2) = 2;  % labels hand

                hold on
                x = 25:25:75;
                y = data(j+5,1:3);
                [r2, b, t, p] = scatterplot(x',y', 'regression', 'linear', 'intercept', 0);
                result((i-1)*10+j+5,3) = r2; 
                result((i-1)*10+j+5,4) = b;
                result((i-1)*10+j+5,5) = t;
                result((i-1)*10+j+5,6) = p;

            end
            title('Right Hand')
            axis([0 75 -5 30])
            % legend('1', '2', '3', '4', '5')

        end
    save('reg_DevF', result);
    

 %       CAT.linecolor = {'r', 'b', 'g', 'k', 'm'};
%         CAT.markercolor = {'r', 'b', 'g', 'k', 'm'};
%         CAT.markerfill = {'r', 'b', 'g', 'k', 'm'};for i=unique(D.SN)'
%             d = getrow(D, D.SN==i);
%             h = figure;
%             subplot(1,2,1),
%             xyplot(d.actForce, d.meanDevP, d.targetForce, 'split', d.digit, 'subset', d.hand==1, 'CAT', CAT, 'leg', 'auto');
%             title('LH'), axis([0 40 0 2.5]), xlabel('actForce'), ylabel('meanDevP')
%             subplot(1,2,2)
%             xyplot(d.actForce, d.meanDevP, d.targetForce, 'split', d.digit, 'subset', d.hand==2, 'CAT', CAT, 'leg', 'auto');
%             title('RH'), axis([0 40 0 2.5]), xlabel('actForce'), ylabel('meanDevP')
%         end
        
    
        
    otherwise
        error('no such case!')
end

end


