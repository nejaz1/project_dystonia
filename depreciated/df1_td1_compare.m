function varargout = df1_td1_compare(what,varargin)
% Direct comparision between fingerPattern_dystonia and tendigit1


df1Dir= fullfile('/Volumes/MacintoshHD2/fingerPattern_dystonia');
td1Dir= fullfile('/Users/jdiedrichsen/Projects/FingerPattern/tendigit1');
df1Dir= fullfile('/Users/joern/Projects/fingerPattern_dystonia');
td1Dir= fullfile('/Users/joern/Projects/FingerPattern/tendigit1');
baseDir={td1Dir,df1Dir};
caretDir={fullfile(td1Dir,'surfaceCaret'),fullfile(df1Dir,'surfaceCaret')};
regDir={fullfile(td1Dir,'RegionOfInterest'),fullfile(df1Dir,'RegionOfInterest')};
subj_name={{'s01','s02','s03','s04','s05','s06'},...
    { 'd01', 's01', 'd02', 's02', 'd04', ...
            's04', 'd06', 's03', 'd07', 'd08', ...
            'd09', 'd10', 'd11'}};
hemName={'LeftHem','RightHem'};
hem={'lh','rh'};
subj_group=[2 1 2 1 2 ...
            1 2 1 2 2 ...
            2 2 2];

switch (what)
    case 'distance'         % Compare mean finger distances across groups
        
        field_dist={'dist_mahal_all','dist_all','dist_sub'};
        field_oth={'SN','hand','region','regSide','regType'};
        pair_name={'12','13','14','15','23','24','25','34','35','45'}'; 
        
        TD=load(fullfile(df1Dir,'RegionOfInterest','reg_distance_comp.mat'));
        T1=load(fullfile(td1Dir,'RegionOfInterest','reg_distance_comp.mat'));
        % Select contralateral Motor activity
        TD.hand=TD.hand+1; % Make hands 1,2: DO THIS IN GENERAL??
        TD=getrow(TD,TD.hand~=TD.regSide & TD.stimtype==0); % RegType seems to be not correctly assigned 
        T1=getrow(T1,T1.hand~=T1.regSide & T1.regType==6);
        % wrangle the data into a new format: In general it may be useful
        % to have this type of structure
        for i=1:length(field_dist)
            x=TD.(field_dist{i}); % Unnormalized
            xs=bsxfun(@rdivide,TD.(field_dist{i}),mean(TD.(field_dist{i}),2));
            TDa.(field_dist{i})=x (:);
            TDa.([field_dist{i} 's'])=xs (:);
            
            x=T1.(field_dist{i}); % Unnormalized
            xs=bsxfun(@rdivide,T1.(field_dist{i}),mean(T1.(field_dist{i}),2));
            T1a.(field_dist{i})=x(:);
            T1a.([field_dist{i} 's'])=xs(:);
        end;
        for i=1:length(field_oth)
            TDa.(field_oth{i})=repmat(TD.(field_oth{i}),10,1);
            T1a.(field_oth{i})=repmat(T1.(field_oth{i}),10,1);
        end;
        TDa.pair=kron([1:10]',ones(length(TD.SN),1));
        T1a.pair=kron([1:10]',ones(length(T1.SN),1));
        
        TDa.group=subj_group(TDa.SN)'+1; % Pianists and dysntonics 
        
        T1a.group=ones(length(T1a.SN),1)*1; % Control
        T=addstruct(TDa,T1a);
        T.pair_name={pair_name{T.pair}}';
        lineplot(T.pair,T.dist_sub,'split',[T.group],'leg','auto','style_thickline');
        set(gca,'XTickLabel',pair_name); 
        varargout={T};  
    case 'test_hand'
        T=varargin{1}; 
        g=varargin{2}; 
        T=getrow(T,isincluded(g,T.group)); 
        lineplot(T.pair,T.dist_mahal_all,'split',T.hand,'leg','auto','style_thickline');
        X=pivottable([T.group T.SN],T.hand,T.dist_mahal_all,'mean'); 
        ttest(X(:,1),X(:,2),2,'paired'); 
        varargout={X}; 
    case 'plot_fingerpatterns' % Plots finger patterns of the two experiments right in the same figure
        sn=varargin{1};     % String of subject numbers
        exp=varargin{2};    % String of experiment numbers
        h=varargin{3};      % Hemisphere
        N=length(sn);       % Number of subjects
        for i=1:N
            groupDir=fullfile(caretDir{exp(i)},'fsaverage_sym',hemName{h});
            cd(groupDir);
            border=fullfile(groupDir,['CS.border']);
            switch(h)
                case 1
                    coord='lh.FLAT.coord';
                    topo='lh.CUT.topo';
                    data='lh.surface_shape';
                    xlims=[-20 10];
                    ylims=[-15 30];
                case 2
                    coord='rh.FLAT.coord';
                    topo='rh.CUT.topo';
                    data='rh.surface_shape';
                    xlims=[-10 20];
                    ylims=[-15 30];
            end;
            
            B=caret_load(border);
            
            data=fullfile(caretDir{exp(i)},['x' subj_name{exp(i)}{sn(i)}],hemName{h},[subj_name{exp(i)}{sn(i)} '_finger.metric']);
            for j=1:5 %----left motor----right motor----right sensory
                subplot(N,5,j+(i-1)*5);
                [M,d]=caret_plotflatmap('col',(2-h)*5+j,'data',data,'cscale',[-6 12],...
                    'border',B.Border,'topo',topo,'coord',coord,'xlims',xlims,'ylims',ylims);
                maxT(j)=max(d(:));
            end;
            mm=max(maxT)
            for j=1:5
                subplot(N,5,j+(i-1)*5);
                caxis([-mm/2 mm]);
            end;
        end;
        set(gcf,'PaperPosition',[1 1 8 1.5*N+0.5]);
        wysiwyg;
    case 'commonROIstruct'         % Make common region of interest structure
        T1=load(fullfile(regDir{1},'reg_data_run.mat'));
        T1=getrow(T1,isincluded([6 12],T1.regNum));
        T2=load(fullfile(regDir{2},'reg_data_BA.mat'));
        T2=getrow(T2,isincluded([6 12],T2.regNum));
        % Add mean activation for the T2: Add this to ROI data at some point
        c=kron(ones(1,8),[1:15]);
        for i=1:15
            T2.meanAct(:,i)=mean(T2.beta(:,find(c==i)),2);
        end;   
        % Make a common structure 
        T.SN=[T1.SN;T2.SN]; 
        T.group=[ones(length(T1.SN),1);ones(length(T2.SN),1)*2]; 
        T.regNum=[T1.regNum;T2.regNum];
        T.xyz=[T1.xyz;T2.xyz];
        T.meanAct=[T1.meanAct(:,6:10);T2.meanAct(:,6:10)]; % Mean activation contra lateral hand 
        T.meanAct(T.regNum==12,:)=[T1.meanAct(T1.regNum==12,1:5);T2.meanAct(T2.regNum==12,1:5)];
        T.ResMs=[T1.ResMs;T2.ResMs]; 
        
        varargout={T};     
    case 'compare_acc'
        subplot(1,2,1); 
        T1=load(fullfile(regDir{1},'reg_acc_8.mat'));
        barplot(T1.regType,T1.accR,'split',T1.regSide); 
        set(gca,'YLim',[0 0.9]); 
        subplot(1,2,2); 
        T2=load(fullfile(regDir{2},'reg_acc_8.mat'));
        barplot(T2.regType,T2.accRm,'split',T2.regSide); 
        set(gca,'YLim',[0 0.9]); 
    case 'compare_acc_2'
        subplot(1,2,1); 
        T1=load(fullfile(regDir{1},'reg_acc_8.mat'));
        barplot(T1.SN,T1.accR,'subset',T1.regSide==1 & T1.regType==1); 
        set(gca,'YLim',[0.25 1]); 
        subplot(1,2,2); 
        T2=load(fullfile(regDir{2},'reg_acc_8.mat'));
        barplot(T2.SN,T2.accRm,'subset',T2.regSide==1 & T2.regType==1); 
        set(gca,'YLim',[0.25 1]); 
    case 'quantiles'         % makes a structure with quantiles for signal/noise values, makes qqplot that shows sparse representation  
        T=df1_td1_compare('commonROIstruct'); 
        T.act=mean(T.meanAct,2)./sqrt(T.ResMs);         % pseudo-tvalue (SNR)
        R=unique([T.group T.SN T.regNum],'rows'); 
        D=[]; 
        for i=1:size(R,1)
            indx=T.group==R(i,1) & T.SN==R(i,2) & T.regNum==R(i,3);
            E.group=R(i,1); 
            E.SN=R(i,2); 
            E.regNum=R(i,3); 
            E.quant=quantile(T.act(indx,:),[0.01:0.01:0.99]); 
            E.quant_std=E.quant./(E.quant(:,75) -E.quant(25));  % This is a relative arbitrary normalization 
                                                 % Although quite standard
                                                 % for qqplots 
            D=addstruct(D,E); 
        end; 
        varargout={D}; 
        
        %Make qqplot with error bars 
        m1=mean(D.quant_std(D.group==1,:)); 
        m2=mean(D.quant_std(D.group==2,:)); 
        traceplot(m1,D.quant,'split',D.group,'errorfcn','stderr','leg','auto'); 
    case 'meanActivation' 
        T=df1_td1_compare('commonROIstruct'); 
        T.act=bsxfun(@rdivide,T.meanAct,sqrt(T.ResMs));         % pseudo-tvalue (SNR)
        barplot(T.group,T.act);
    case 'spatial_kernel' 
        T=df1_td1_compare('commonROIstruct'); 
        T.act=mean(T.meanAct,2)./sqrt(T.ResMs);         % pseudo-tvalue (SNR)
        T.actF=bsxfun(@rdivide,T.meanAct,sqrt(T.ResMs)); % Pseudo t-value for each Finger
        T.actFd=bsxfun(@minus,T.actF,mean(T.actF,2));    % Mean subtracted for each voxel (
        
        R=unique([T.group T.SN T.regNum],'rows'); 
        D=[]; 
        border=[0 2 3 4 6 8 10 12 16 20 ]; 
        for i=1:size(R,1)
            indx=T.group==R(i,1) & T.SN==R(i,2) & T.regNum==R(i,3);
            E.group=R(i,1); 
            E.SN=R(i,2); 
            E.regNum=R(i,3); 
            actFdd=bsxfun(@minus,T.actFd(indx,:),mean(T.actFd(indx,:))); 
            E.spatial=mva_spatial_kernel(T.xyz(indx,:),T.actF(indx,:)',border); 
            [E.spatial_f,E.sig,E.dist]=mva_spatial_kernel(T.xyz(indx,:),T.actFd(indx,:)',border); 
            [E.spatial_m,E.sigm,E.dist]=mva_spatial_kernel(T.xyz(indx,:),actFdd',border); 
            D=addstruct(D,E); 
        end;         
        varargout={D}; 
        traceplot(mean(D.dist),D.spatial_f,'split',[D.group],'errorfcn','stderr','subset',~[D.SN==3 & D.group==2 & D.regNum==12]);
 
end;