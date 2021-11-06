function R1_Artur_Trajectory_glut
clc; close all;
% data related settings
match = {'AC'};
CCTh = 1;
cols = getColors;
selSCC = 3;

[~, name] = system('hostname');
if strcmp(strtrim(name),'mohajem3-psy13')
    pathPrefix = 'Y:\homes\navvab.afrashteh\CloudStation\Navvab-Sam\NeuroPhotonic Paper\Test\First Revision\';
end
pathPrefix = 'Y:\homes\navvab.afrashteh\CloudStation\Navvab-Sam\NeuroPhotonic Paper\Test\First Revision\';

TrajMethod = 'FE_StretchContour';

ppType = 'BPF';
% ppType = 'HPF';
ss = 1;
calcTrajFlag = 0;
if calcTrajFlag
    % calculate and save the trajectories for FL and HL
    FindSaveTraj(ss,ppType);
end

plotTrajFlag = 1;
if plotTrajFlag
    pdfFlag = 1;
    pdfOpen = 1;
    fineSp = 0.001;
    annotationFontSize = 8; axesFontSize = 6; legendFontSize = 8; xlabelFontSize = 9; ylabelFontSize = 9;
    titleFontSize = 8;
    plotTraj(ss,ppType)
end

return;

    function FindSaveTraj(ss,ppType)
        % load data from base work space where all the data has to be pre loaded
        HPF = 0;
        if strcmpi(ppType,'HPF')
            mdataE = evalin('base','mdataE_glut_HPF');
            mdataS = evalin('base','mdataS_glut_HPF');
            HPF = 1;
        elseif strcmpi(ppType,'BPF')
            mdataE = evalin('base','mdataE_glut_BPF_050Hz_6Hz');
            mdataS = evalin('base','mdataS_glut_BPF_050Hz_6Hz');
        end
        meta = evalin('base','meta');
        animalNumber = mdataE{ss}.animalNumber;
        selLMH_all = mdataE{ss}.selLMH.selLMH_all;
        animals = [1 2 3 4 5];
        amplitude_threshold = meta.amplitude_threshold;
        roi = 1;
        dfth = 2;
        
        for an = 1:length(animals)
            if HPF
                valsE{an} = getDistributions_simpleAmps_artur_2_glut_HPF(mdataE{ss},animals(an),roi(ss),ss,amplitude_threshold);
                valsS{an} = getDistributionsSpon_simple_artur_glut(mdataS{1},animals(an),roi(ss),ss,dfth,amplitude_threshold);
                for LMH = 3%length(ssLMH)
                    valsSCC{an,LMH} = getDistributionsSpon_motifs_simple_artur_glut_2_HPF (mdataE,mdataS{1},animals(an),roi(ss),'AC','Max',CCTh,LMH,amplitude_threshold);
                end
            else
                valsE{an} = getDistributions_simpleAmps_artur_2_glut(mdataE{ss},animals(an),roi(ss),ss,amplitude_threshold);
%                 valsS{an} = getDistributionsSpon_simple_artur_glut(mdataS{1},animals(an),roi(ss),ss,dfth,amplitude_threshold);
                for LMH = 3%length(ssLMH)
                    valsSCC{an,LMH} = getDistributionsSpon_motifs_simple_artur_glut_2 (mdataE,mdataS{1},animals(an),roi(ss),'AC','Max',CCTh,LMH,amplitude_threshold);
                end
            end
        end
        
        multiplier = 0.5;
        r = 20; % 15 pixels = 1mm. radius for roi around center of stim modality. used to evaluate if the start point of activity falls in an roi around the stim modality.
        rM = 15*3;
        
        for an = 1:length(animals)
            % Masks
            mask = mdataE{ss}.masks{an}{1};
            MaskStim = mdataE{ss}.masks{an}{ss+1};
            [xMask,yMask] = findMaskBorder(mask);
            center = regionprops(MaskStim, MaskStim, 'Centroid'); center = center.Centroid;
            xc = center(1); yc = center(2);
            
            anPathName = makeName(animalNumber(an),getMainDataFolderGlut);
            sPathName = makeName('pSpon',anPathName);
            ePathName = makeName('pEvoked',anPathName);
            
            % Traj for Evoked
            [lists, names] = getEvokedListsDiffStimAmps_glut(animalNumber(an),1);
            for level = 1:length(names{4})
                onsets = valsE{an}.onsets(level,:);
                pfs = valsE{an}.pfs(level,:);
                trials = getEvokedListsDiffStimAmpsTrials_glut(animalNumber(an),names{4}{level});
                fPathName = makeName(names{4}(level),ePathName);
                selTrials = 1:length(trials);
                TrialsTrajInfo = cell(1,length(selTrials));
                ne = 0;
                eX = {};
                eY = {};
                eXinterp = [];
                eYinterp = [];
                Trials = [];
                for jj = 1:length(selTrials)
                    try
                        f1 = onsets(jj)+109;
                        f2 = pfs(jj)+109;
                        if isnan(f1) || isnan(f2) || ~and(f1,f2)
                            continue;
                        end
                    catch
                        continue;
                    end
                    trialFolder = makeName(trials{selTrials(jj)},fPathName);
                    if HPF
                        fileName = makeName(sprintf('ImgSeq_HPF.mat'),trialFolder);
                    else
                        fileName = makeName(sprintf('ImgSeq_BPF_050Hz_6Hz.mat'),trialFolder);
                    end
                    
                    ImgSeq = load(fileName);
                    ImgSeq = ImgSeq.ImgSeq;
                    ImgSeq = ImgSeq(:,:,f1:f2);
                    
                    
%                     TrialsTrajInfo{jj} = calcAllCenterOfMassTrajectory6(ImgSeq,mask,multiplier,xc,yc,rM);
%                     TrialsTrajInfo{jj}.x = inpaint_nans(TrialsTrajInfo{jj}.x,3);
%                     TrialsTrajInfo{jj}.y = inpaint_nans(TrialsTrajInfo{jj}.y,3);
%                     figure(1);clf;imagesc(mask);axis image;hold on
%                     plot(TrialsTrajInfo{jj}.x,TrialsTrajInfo{jj}.y,'b','linewidth',2)
%                     TrialsTrajInfo{jj} = calcAllCenterOfMassTrajectory7(ImgSeq,mask,multiplier,xc,yc,rM);
%                     TrialsTrajInfo{jj}.x = inpaint_nans(TrialsTrajInfo{jj}.x,3);
%                     TrialsTrajInfo{jj}.y = inpaint_nans(TrialsTrajInfo{jj}.y,3);
%                     plot(TrialsTrajInfo{jj}.x,TrialsTrajInfo{jj}.y,'r','linewidth',2)
%                     TrialsTrajInfo{jj} = calcAllCenterOfMassTrajectory8(ImgSeq,mask,multiplier,xc,yc,rM);
%                     TrialsTrajInfo{jj}.x = inpaint_nans(TrialsTrajInfo{jj}.x,3);
%                     TrialsTrajInfo{jj}.y = inpaint_nans(TrialsTrajInfo{jj}.y,3);
%                     plot(TrialsTrajInfo{jj}.x,TrialsTrajInfo{jj}.y,'k','linewidth',2)
%                     TrialsTrajInfo{jj} = calcAllCenterOfMassTrajectory9(ImgSeq,mask,multiplier,xc,yc,rM);
%                     TrialsTrajInfo{jj}.x = inpaint_nans(TrialsTrajInfo{jj}.x,3);
%                     TrialsTrajInfo{jj}.y = inpaint_nans(TrialsTrajInfo{jj}.y,3);
%                     plot(TrialsTrajInfo{jj}.x,TrialsTrajInfo{jj}.y,'w','linewidth',2)
%                     
                    TrialsTrajInfo{jj} = calcAllCenterOfMassTrajectory7(ImgSeq,mask,multiplier,xc,yc,rM);
                    TrialsTrajInfo{jj}.x = inpaint_nans(TrialsTrajInfo{jj}.x,3);
                    TrialsTrajInfo{jj}.y = inpaint_nans(TrialsTrajInfo{jj}.y,3);
                    xS = TrialsTrajInfo{jj}.x(1);
                    yS = TrialsTrajInfo{jj}.y(1);
                    in = isInOnCircle(xS,yS,xc,yc,r);
                    in = 1;
                    if (sum(TrialsTrajInfo{jj}.x <max(xMask))/length(TrialsTrajInfo{jj}.x) < 1) || (sum(TrialsTrajInfo{jj}.y <max(yMask))/length(TrialsTrajInfo{jj}.y) < 1)
                        in = 0;
                    end
                    if (sum(TrialsTrajInfo{jj}.x >min(xMask))/length(TrialsTrajInfo{jj}.x) < 1) || (sum(TrialsTrajInfo{jj}.y >min(yMask))/length(TrialsTrajInfo{jj}.y) < 1)
                        in = 0;
                    end
                    if in
                        ne = ne+1;
                        eX{ne} = TrialsTrajInfo{jj}.x;
                        eY{ne} = TrialsTrajInfo{jj}.y;
                        
                        idx = 1:length(TrialsTrajInfo{jj}.x);
                        idxM = max(idx);
                        idxInterp = linspace(1,idxM,50);
                        eXinterp(ne,:) = interp1(idx,TrialsTrajInfo{jj}.x,idxInterp);
                        eYinterp(ne,:) = interp1(idx,TrialsTrajInfo{jj}.y,idxInterp);
                        Trials(ne) = jj;
                    end
                end
                Etraj{an}.X{level} = eX;
                Etraj{an}.Y{level} = eY;
                Etraj{an}.Xinterp{level} = eXinterp;
                Etraj{an}.Yinterp{level} = eYinterp;
                Etraj{an}.Trials{level} = Trials;
            end
            
            % Traj for Spontaneous
            LMH = 3;
            listR = listOfRecordingsSpon_glut(animalNumber(an));
            ns = 0;
            sX = [];
            sY = [];
            sXinterp = [];
            sYinterp = [];
            TrajFrames = [];
            TrialsTrajInfo = cell(1,size(valsSCC{an,LMH}.allFrames,1));
            jj = 1;
            
            preFN = 100;
            postFN = 100;
            ImgSeq_avg = zeros(128,128,preFN+postFN+1);
            mmTot = 0;
            for ii = 1%:length(listR)
                if HPF
                    DF_F0PathName = makeName(sprintf('%s.mat',listR{3}),anPathName);
                    DF_F0 = load(DF_F0PathName);
                    DF_F0 = DF_F0.DF_F0;
                else
                    DF_F0PathName = makeName(sprintf('%s.mat',listR{2}),anPathName);
                    DF_F0 = load(DF_F0PathName);
                    DF_F0 = DF_F0.DF_F0;
                end
                mm = 0;
                CorrMaxDF_frames_temp = [];
                for motif = 1:size(valsSCC{an,LMH}.allFrames,1)
%                     if valsSCC{an,LMH}.allFrames(motif,2)~= ii
%                         continue;
%                     end
                    f1 = valsSCC{an,LMH}.allFrames(motif,4);
                    f2 = valsSCC{an,LMH}.allFrames(motif,3);
                    ImgSeq = DF_F0(:,:,f1:f2);
                    
                    TrialsTrajInfo{jj} = calcAllCenterOfMassTrajectory7(ImgSeq,mask,multiplier,xc,yc,rM);
                    TrialsTrajInfo{jj}.x = inpaint_nans(TrialsTrajInfo{jj}.x,3);
                    TrialsTrajInfo{jj}.y = inpaint_nans(TrialsTrajInfo{jj}.y,3);
                    xS = TrialsTrajInfo{jj}.x(1);
                    yS = TrialsTrajInfo{jj}.y(1);
                    [in, xyCirc] = isInOnCircle(xS,yS,xc,yc,r);
                    in = 1;
                    if (sum(TrialsTrajInfo{jj}.x <max(xMask))/length(TrialsTrajInfo{jj}.x) < 1) || (sum(TrialsTrajInfo{jj}.y <max(yMask))/length(TrialsTrajInfo{jj}.y) < 1)
                        in = 0;
                    end
                    if (sum(TrialsTrajInfo{jj}.x >min(xMask))/length(TrialsTrajInfo{jj}.x) < 1) || (sum(TrialsTrajInfo{jj}.y >min(yMask))/length(TrialsTrajInfo{jj}.y) < 1)
                        in = 0;
                    end
                    if in
                        ns = ns+1;
                        sX{ns} = TrialsTrajInfo{jj}.x;
                        sY{ns} = TrialsTrajInfo{jj}.y;
                        idx = 1:length(TrialsTrajInfo{jj}.x);
                        idxM = max(idx);
                        idxInterp = linspace(1,idxM,50);
                        sXinterp(ns,:) = interp1(idx,TrialsTrajInfo{jj}.x,idxInterp);
                        sYinterp(ns,:) = interp1(idx,TrialsTrajInfo{jj}.y,idxInterp);
                        TrajFrames(ns,:) = valsSCC{an,LMH}.allFrames(motif,:);
                        try
                            ImgSeq_avg = ImgSeq_avg + DF_F0(:,:,f2-preFN:f2+postFN);
                            mm = mm+1;
                            CorrMaxDF_frames_temp(mm) = f2;
                        catch
                        end
                    end
                    jj = jj+1;
                end
                CorrMaxDF_frames{an,ii} = CorrMaxDF_frames_temp;
                mmTot = mmTot + mm;
            end
            Straj{an}.X = sX;
            Straj{an}.Y = sY;
            Straj{an}.Xinterp = sXinterp;
            Straj{an}.Yinterp = sYinterp;
            Straj{an}.TrajFrames = TrajFrames;
            
            %{
            ImgSeq = ImgSeq_avg/mmTot;
            if HPF
                fileName =  ['ImgSeq_HPF_avg_an', num2str(an),'_' ,match{ss}, '_glut_all_04corr.raw'];
            else
                fileName =  ['ImgSeq_BPF_05_6Hz_avg_an', num2str(an),'_' ,match{ss}, '_glut_all_04corr.raw'];
            end
            imwriteallraw([pathPrefix '\' fileName],ImgSeq,'*float32')
            
            if HPF
                fileName =  ['ImgSeq_HPF_avg_an', num2str(an),'_' ,match{ss}, '_glut_all_04corr.mat'];
            else
                fileName =  ['ImgSeq_BPF_05_6Hz_avg_an', num2str(an),'_' ,match{ss}, '_glut_all_04corr.mat'];
            end
            save([pathPrefix '\' fileName],'ImgSeq','-v7.3')
            
            mask1.bigMask = getMask(animalNumber(an),'ec',1.1,'glut',1,'ac',0);
            mask1.mask = getMask(animalNumber(an),'glut',1,'ac',0);
            CLG_on_file(pathPrefix,fileName, 'ImgSeq', mask1)
            
            if HPF
                CorrMaxDF_frames_fileName = [pathPrefix '\CorrMaxDF_HPF_frames','_', match{ss}, '_glut_all_04corr.mat'];
                fileName =  ['ImgSeq_HPF_avg_an', num2str(an),'_' ,match{ss}, '_glut_all_04corr.mat'];
            else
                CorrMaxDF_frames_fileName = [pathPrefix '\CorrMaxDF_BPF_05_6Hz_frames','_', match{ss}, '_glut_all_04corr.mat'];
                fileName =  ['ImgSeq_BPF_05_6Hz_avg_an', num2str(an),'_' ,match{ss}, '_glut_all_04corr.mat'];
            end
            save(CorrMaxDF_frames_fileName,'CorrMaxDF_frames')
            %}
        end
        % save results
        if HPF
            EtrajName = 'Etraj_HPF_';
        else
            EtrajName = 'Etraj_BPF_05_6Hz_';
        end
        saveEtraj_FileName = sprintf('%s\\%s%s_allAn_04corr.mat',pathPrefix ,EtrajName,match{ss});
        save(saveEtraj_FileName,'Etraj')
        
        if HPF
            StrajName = 'Straj_HPF_';
        else
            StrajName = 'Straj_BPF_05_6Hz_';
        end
        saveStraj_FileName = sprintf('%s\\%s%s_allAn_04corr.mat',pathPrefix ,StrajName,match{ss});
        save(saveStraj_FileName,'Straj')
    end
    
    function plotTraj(ss,ppType)
        nRows = 3; nCols = 5; sprc = [10 10]; whUnits = [-15 -15]; rtup = [30 0]; rbUnits = [0 0];
        hf = makeFigureWindow(1002,6.9,nRows*1.5,1);
%         ff = makeFigureRowsCols(111,[6 2 5.9 1.75],'RowsCols',[1 4],'spaceRowsCols',[0.1 0.03],'rightUpShifts',[0.06 0.01],'widthHeightAdjustment',[-40 -80]);
        ff = makeFigureRowsCols(1112,[3,3 4.5135 1.75],'RowsCols',[1 3],'spaceRowsCols',[0.1 0.0395],'rightUpShifts',[0.0782 0.01],'widthHeightAdjustment',[-59 -80]);
        rr = 1; 
        [pL pB pW pH] = getPanelProps(nRows,nCols,sprc,whUnits,rtup,rbUnits);
        % load data from base work space where all the data has to be pre loaded
        HPF = 0;
        if strcmpi(ppType,'HPF')
            mdataE = evalin('base','mdataE_glut_HPF');
            mdataS = evalin('base','mdataS_glut_HPF');
            HPF = 1;
        elseif strcmpi(ppType,'BPF')
            mdataE = evalin('base','mdataE_glut_BPF_050Hz_6Hz');
            mdataS = evalin('base','mdataS_glut_BPF_050Hz_6Hz');
        end
        meta = evalin('base','meta');
%         mdataS = evalin('base','mdataS_glut_HPF');
%         meta = evalin('base','meta');
        Cmap = evalin('base','Cmap');
        animalNumber = mdataE{ss}.animalNumber;
        animals = [1 2 3 4 5];
        sel_an = 4;
        selLMH_all = mdataE{ss}.selLMH.selLMH_all;
        % load Trajectory results
        if HPF
            EtrajName = 'Etraj_HPF_';
        else
            EtrajName = 'Etraj_BPF_05_6Hz_';
        end
        saveEtraj_FileName = sprintf('%s\\%s%s_allAn_04corr.mat',pathPrefix ,EtrajName,match{ss});
        load(saveEtraj_FileName,'Etraj')
        if HPF
            StrajName = 'Straj_HPF_';
        else
            StrajName = 'Straj_BPF_05_6Hz_';
        end
        saveStraj_FileName = sprintf('%s\\%s%s_allAn_04corr.mat',pathPrefix ,StrajName,match{ss});
        load(saveStraj_FileName,'Straj')
        % Create axes
        for ann = 1:length(animals)
            an = animals(ann);
            selLMH = selLMH_all{ss}(an,:);            
            paPositionBoth(1,:) = [pL(1,an) pB(1,an) pW pH];
            paPositionBoth(2,:) = [pL(2,an) pB(2,an) pW pH];
            paPositionBoth(3,:) = [pL(3,an) pB(3,an) pW pH];
            anPathName = makeName(animalNumber(an),getMainDataFolderGlut);
            BregmaFileName = makeName('bregmaXY.txt',anPathName);
            if ~exist(BregmaFileName,'file')
            else
                bregma = load(BregmaFileName);
            end
            mPathName = makeName('Mask.mat',anPathName);
            mask = load(mPathName); mask = mask.mask;
            [xMask, yMask] = findMaskBorder(mask);
            
            method = 'bicubic';
            param.Sigma = 2;
            param.GaussWindow = 4*param.Sigma+1;
            Center  = fix([param.GaussWindow/2,param.GaussWindow/2])+1;
            [R,C]   = ndgrid(1:param.GaussWindow, 1:param.GaussWindow);
            Gauss2d = gauss2dC(R,C,param.Sigma,Center);
            % Evoked
            LMHnames = {'Lo', 'Med', 'Hi'};
            LH = 0;
            for LMH = [1,3]
                LH = LH+1;
                ax1 = axes('parent',hf, 'Position',paPositionBoth(LH+1,:));hold on;
                [dim1, dim2] = size(mask);
                imE = zeros(dim1, dim2);
                eX = Etraj{an}.Xinterp{selLMH(LMH)};
                eY = Etraj{an}.Yinterp{selLMH(LMH)};
                if LH==1
                    eX = [eX;Etraj{an}.Xinterp{selLMH(LMH+1)}];
                    eY = [eY;Etraj{an}.Yinterp{selLMH(LMH+1)}];
                end
                Ic = zeros(dim1,dim2);
                for jj = 1:size(eX,1)
                    I = zeros(dim1, dim2);
                    xi = eX(jj,:);
                    yi = eY(jj,:);
                    xd = diff(xi);
                    yd = diff(yi);
                    TrajLen = sqrt(xd.^2 + yd.^2);
                    TrajLen = sum(TrajLen);
                    np = round(TrajLen*3);
                    try
                        [cx,cy,~] = improfile(I,xi,yi,np,method);
                    catch
                        cx = [];cy = [];
                    end
                    cx = round(cx); cy = round(cy);
                    I((cx-1)*dim1+cy) = 1;
                    Ic = Ic + I;
                    I = conv2(I,Gauss2d,'same');
                    imE = imE + I;
                end
                imE = imE.*mask;
                im = imE;
%                 im = im./size(eX,1); %Sam added to normalize the trajectory plot by the number of trials or motifs
                maxI = max(im(:));
                maxI = max(maxI,1);
                level = graythresh(im);
                imbw = im2bw(im,level);
                
                axes(ax1); imagesc(im,[0,maxI]); colormap(Cmap); hold on; axis off; axis equal; xlim([1 128]);
                plot(xMask, yMask, 'w'); %plot(bregma(1), bregma(2), 'mo','markersize',2,'markerfacecolor','m');
                set(ax1,'ydir','reverse')
                if sel_an == an
                    cc = 1+LH;axes(ff.h_axes(rr,cc));
                    imagesc(im,[0,maxI]); colormap(Cmap); hold on; axis off; axis equal; xlim([1 128]);
                    plot(xMask, yMask, 'w');% plot(bregma(1), bregma(2), 'mo','markersize',2,'markerfacecolor','m');
                    set(ff.h_axes(rr,cc),'ydir','reverse')
                    if cc == 3
                        cc = 4;
                    end
                    textString = meta.stimNames{cc};
                    text(60,-10,textString,'fontsize',meta.legendFontSize+3,'Rotation',0,'FontWeight','Bold',...
                        'HorizontalAlignment','center','VerticalAlignment','middle');
                end
                axes(ax1);
                if an == 1
                    axes(ax1);
                    xt = -5; yt = 64.5;
                    if ss == 3
                        txtStr = sprintf('%s Evoked',match{ss});
                    else
%                         txtStr = sprintf('%s %s Evoked',LMHnames{LMH},match{ss});
                        txtStr = meta.stimNames{LMH+1};
                        thisCols = {cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3},cols.stimColors{4}};
                    end
                    ht = text(xt,yt,txtStr,'HorizontalAlignment','center','VerticalAlignment','middle','rotation',90,'color',thisCols{LMH+1});
                    
                end
                FD_E(an,LH) = BoxCountfracDim(imbw);
            end
            FD.FD_E = FD_E;
            % Spon
            ax2 = axes('parent',hf, 'Position',paPositionBoth(1,:));hold on;
            imS = zeros(dim1, dim2);
            sX = Straj{an}.Xinterp;
            sY = Straj{an}.Yinterp;
            Ic = zeros(dim1,dim2);
            for jj = 1:size(sX,1)
                I = zeros(dim1, dim2);
                xi = sX(jj,:);
                yi = sY(jj,:);
                xd = diff(xi);
                yd = diff(yi);
                TrajLen = sqrt(xd.^2 + yd.^2);
                TrajLen = sum(TrajLen);
                np = round(TrajLen*3);
                try
                    [cx,cy,~] = improfile(I,xi,yi,np,method);
                catch
                    cx = [];cy = [];
                end
                cx = round(cx); cy = round(cy);
                I((cx-1)*dim1+cy) = 1;
                Ic = Ic + I;
                I = conv2(I,Gauss2d,'same');
                imS = imS + I;
            end
            imS = imS.*mask;
            im = imS;
%             im = im./size(sX,1); %Sam added to normalize the trajectory plot by the number of trials or motifs
            maxI = max(im(:));
            maxI = max(maxI,1);
            level = graythresh(im);
            imbw = im2bw(im,level);
            FD_S(an) = BoxCountfracDim(imbw);
            
            axes(ax2); imagesc(im,[0,maxI]); colormap(Cmap); hold on; axis off; axis equal; xlim([1 128]);
            plot(xMask, yMask, 'w'); %plot(bregma(1), bregma(2), 'mo','markersize',2,'markerfacecolor','m');
            set(ax2,'ydir','reverse')
            xt = 64.5; yt = -5;
            text(xt,yt,['Animal #',num2str(an)],'HorizontalAlignment','center','VerticalAlignment','middle','fontsize',8);
            if sel_an == an
                cc = 1;axes(ff.h_axes(rr,cc));
                imagesc(im,[0,maxI]); colormap(Cmap); hold on; axis off; axis equal; xlim([1 128]);
                plot(xMask, yMask, 'w'); %plot(bregma(1), bregma(2), 'mo','markersize',2,'markerfacecolor','m');
                set(ff.h_axes(rr,cc),'ydir','reverse')

                textString = meta.stimNames{cc};
                text(60,-10,textString,'fontsize',meta.legendFontSize+3,'Rotation',0,'FontWeight','Bold',...
                            'HorizontalAlignment','center','VerticalAlignment','middle');
            end
            axes(ax2);
            if an == 1
                xt = -5; yt = 64.5;
                txtStr = 'Spon';
                text(xt,yt,txtStr,'HorizontalAlignment','center','VerticalAlignment','middle','rotation',90,'color',cols.sponColorsN{1})
            end
            if an == 2
                xL = 5;
                yL = 120;
                plot(xL+[0 30], yL*[1 1], 'w','linewidth',1.5)
                text(xL,yL-10,'2 mm','color','w','FontSize',axesFontSize)

                midX = 112;
                midY = 21;
                hL = 7; % half length of the cross
                sp = 3; % space between letters and the ends of cross
                plot(midX+[-hL hL], midY*[1 1], 'w','linewidth',1)
                plot(midX*[1 1], midY+[-hL hL], 'w','linewidth',1)
                text(midX,midY-hL-sp,'A','color','w','FontSize',axesFontSize,'HorizontalAlignment','center','VerticalAlignment','baseline');
                text(midX,midY+hL+sp,'P','color','w','FontSize',axesFontSize,'HorizontalAlignment','center','VerticalAlignment','cap');
                text(midX+hL+sp,midY,'L','color','w','FontSize',axesFontSize,'HorizontalAlignment','left','VerticalAlignment','middle');
                text(midX-hL-sp,midY,'M','color','w','FontSize',axesFontSize,'HorizontalAlignment','right','VerticalAlignment','middle');
            end
        end
        FD.FD_S = FD_S;
%         load(sprintf('Fractional Analysis_%s.mat',match{ss}),'FD')        
        pdfn = 2*[101 101];
        pfileName = deletePDF(pdfn(ss),ss);
        pfileName = makeName(pfileName,getpdffolder);
        save2pdf(pfileName,ff.hf,600);
        % plot fractal dimension bar graph
        ff = makeFigureWindow__one_axes_only(2,[1 4 1 1.75],[0.38 0.2 0.6 0.69]);
        plotBarGraphFDs (ff,[1 2 3 4 5],FD,meta.colors,meta)
        append_pdfs(pfileName,makeName('temp.pdf',getpdffolder))
        save2pdf(makeName('temp.pdf',getpdffolder),hf,600);
        append_pdfs(pfileName,makeName('temp.pdf',getpdffolder));
        n= 0;
    end
    function plotBarGraphFDs (ff,animals,FD,cols,meta)
        % figure related settings
        annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; xlabelFontSize = 6; ylabelFontSize = 6;
        titleFontSize = 8;
        LH = [1,3];
        for an = 1:length(animals)
            for iii = 1:length(LH)
                avgValsE(an,iii) = FD.FD_E(an,iii);
            end
            avgValsSCC(an) = FD.FD_S(an);
        end
        figure(1000);
        post_hoc = 'hsd';
        annovaVar = [avgValsSCC' avgValsE(:,1) avgValsE(:,2)];
        disp('Fractal Dimension')
        prT = ranovaGlut(annovaVar);
        thisCols = {cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{3},cols.stimColors{4}};
        
%         [p,tbl,stats] = anova1(annovaVar,{'SpontCC','Evok_1','Evok_3'},'off');
%         [c,~,~,gnames] = multcompare(stats,'CType',post_hoc);
%         prT = c(:,6);
        
        hrT = prT<0.05;
        annovaVar = [avgValsSCC' avgValsE(:,1) avgValsE(:,2)];
        axes(ff.ha);
        [mVals semVals] = findMeanAndStandardError(annovaVar,1);
        combs = nchoosek(1:size(annovaVar,2),2);
        plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols,'ySpacingFactor',10);
        xlim([0.4 0.6+size(annovaVar,2)]);
        YLim = get(gca,'ylim'); YLim(1) = 1; ylim(YLim)
        hy = ylabel('Fractal Dimension (A.U.)');
        pos = get(hy,'Position');pos(1) = pos(1)+0.1;set(hy,'Position',pos);
        set(ff.ha,'TickDir','out','FontSize',axesFontSize);
        set(ff.ha,'XTick',1:size(annovaVar,2));
        set(ff.ha,'XTickLabel',meta.stimNames([1,2,4]));
        xtickangle(ff.ha,45);
        set(ff.ha,'Fontweight','bold');
        save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);
    end
% Tracking front-end of wave. Direction based on how contour stretch between two consecutive frames
    function TrajInfo = calcAllCenterOfMassTrajectory6(ImgSeq,mask,multiplier,x0,y0,r)
        [dim1,dim2,~] = size(ImgSeq);
        xc = x0; yc = y0;
        ImgSeq = applyMask(ImgSeq,mask,0);
        
        idxS = 1; idxE = size(ImgSeq,3);
        idxs = idxS:idxE;
        for ii = idxs
            ROI = CircROI(xc,yc,r,dim1,dim2);
            
            ii1 = max(ii-1,1);
            im = ImgSeq(:,:,ii1);
            maskVals = getMaskedValues(im,ROI);
            threshold = mean(maskVals) + multiplier*std(maskVals);
            if isnan(threshold); threshold = 0; end
            if isinf(threshold); threshold = 1; end
            if threshold > 1; threshold = 1; end
            threshold = max(0,threshold);
            [BlargestP, xc, yc] = getlargestROI(ImgSeq(:,:,ii1).*ROI,threshold);
            
            ii2 = ii;
            im = ImgSeq(:,:,ii2);
            maskVals = getMaskedValues(im,ROI);
            threshold_s = mean(maskVals) + 0.25*multiplier*std(maskVals);
            if isnan(threshold_s); threshold_s = 0; end
            if isinf(threshold_s); threshold_s = 1; end
            if threshold_s > 1; threshold_s = 1; end
            threshold_s = max(0,threshold_s);
            [Blargest_s, x1, y1] = getlargestROI(ImgSeq(:,:,ii).*ROI,threshold_s);
            
            nDir = 4*5;
            flowdir = getFlowDir(Blargest_s{1}, x1, y1,BlargestP{1}, xc, yc,nDir);
            
            ii3 = ii;
            im = ImgSeq(:,:,ii3);
            maskVals = getMaskedValues(im,ROI);
            threshold_s = mean(maskVals) + 0.6*multiplier*std(maskVals);
            if isnan(threshold_s); threshold_s = 0; end
            if isinf(threshold_s); threshold_s = 1; end
            if threshold_s > 1; threshold_s = 1; end
            threshold_s = max(0,threshold_s);
            [Blargest, xc, yc] = getlargestROI(ImgSeq(:,:,ii).*ROI,threshold_s);
            
            [xe,ye] = FlowPoint(Blargest{1},x0,y0,flowdir);
            
            x0 = xc; y0 = yc;
            
            if ii == 1;
                TrajInfo.x(ii) = xc;
                TrajInfo.y(ii) = yc;
            end
            if isempty(xe) || isempty(ye)
                xe = nan; ye = nan;
            end
            TrajInfo.x(ii+1) = xe;
            TrajInfo.y(ii+1) = ye;
        end
        
        idx = find(~isnan(TrajInfo.x));
        if length(idx)>1
            TrajInfo.x = interp1(idx,TrajInfo.x(idx),1:max(idxs)+1);
        end
        idx = find(~isnan(TrajInfo.y));
        if length(idx)>1
            TrajInfo.y = interp1(idx,TrajInfo.y(idx),1:max(idxs)+1);
        end
    end

    function TrajInfo = calcAllCenterOfMassTrajectory7(ImgSeq,mask,multiplier,x0,y0,r)
        [dim1,dim2,~] = size(ImgSeq);
        xc = x0; yc = y0;
        ImgSeq = applyMask(ImgSeq,mask,0);
        idxS = 1; idxE = size(ImgSeq,3);
        idxs = idxS:idxE;
        for ii = idxs            
            ROI = mask;%CircROI(xc,yc,r,dim1,dim2).*mask;
            ii1 = max(ii-1,1);
            im = ImgSeq(:,:,ii1);
            maskVals = getMaskedValues(im,ROI);
            threshold = prctile(maskVals,40);%40
            if isnan(threshold); threshold = 0; end
            if isinf(threshold); threshold = 1; end
            if threshold > 1; threshold = 1; end
            threshold = max(0,threshold);
            [BlargestP, xc, yc] = getlargestROI(ImgSeq(:,:,ii1).*ROI,threshold);
            ii2 = ii;
            im = ImgSeq(:,:,ii2);
            maskVals = getMaskedValues(im,ROI);
            threshold_s = prctile(maskVals,25);%25
            if isnan(threshold_s); threshold_s = 0; end
            if isinf(threshold_s); threshold_s = 1; end
            if threshold_s > 1; threshold_s = 1; end
            threshold_s = max(0,threshold_s);
            [Blargest_s, x1, y1] = getlargestROI(ImgSeq(:,:,ii).*ROI,threshold_s);
            nDir = 4*5;
            flowdir = getFlowDir(Blargest_s{1}, x1, y1,BlargestP{1}, xc, yc,nDir);
            ii3 = ii;
            im = ImgSeq(:,:,ii3);
            maskVals = getMaskedValues(im,ROI);
            threshold_s = prctile(maskVals,50);%50
            if isnan(threshold_s); threshold_s = 0; end
            if isinf(threshold_s); threshold_s = 1; end
            if threshold_s > 1; threshold_s = 1; end
            threshold_s = max(0,threshold_s);
            [Blargest, xc, yc] = getlargestROI(ImgSeq(:,:,ii).*ROI,threshold_s);
            [xe,ye] = FlowPoint(Blargest{1},x0,y0,flowdir);
            x0 = xc; y0 = yc;            
            if ii == 1
                TrajInfo.x(ii) = xc;
                TrajInfo.y(ii) = yc;
            end
            if isempty(xe) || isempty(ye)
                xe = nan; ye = nan;
            end
            TrajInfo.x(ii+1) = xe;
            TrajInfo.y(ii+1) = ye;
        end
        idx = find(~isnan(TrajInfo.x));
        if length(idx)>1
            TrajInfo.x = interp1(idx,TrajInfo.x(idx),1:max(idxs)+1);
        end
        idx = find(~isnan(TrajInfo.y));
        if length(idx)>1
            TrajInfo.y = interp1(idx,TrajInfo.y(idx),1:max(idxs)+1);
        end
    end

    function TrajInfo = calcAllCenterOfMassTrajectory8(ImgSeq,mask,multiplier,x0,y0,r)
        [dim1,dim2,~] = size(ImgSeq);
        xc = x0; yc = y0;
        ImgSeq = applyMask(ImgSeq,mask,0);
        idxS = 1; idxE = size(ImgSeq,3);
        idxs = idxS:idxE;
        for ii = idxs            
            ROI = CircROI(xc,yc,r,dim1,dim2).*mask;
            ii1 = max(ii-1,1);
            im = ImgSeq(:,:,ii1);
            maskVals = getMaskedValues(im,ROI);
            threshold = prctile(maskVals,50);
            if isnan(threshold); threshold = 0; end
            if isinf(threshold); threshold = 1; end
            if threshold > 1; threshold = 1; end
            threshold = max(0,threshold);
            [BlargestP, xc, yc] = getlargestROI(ImgSeq(:,:,ii1).*ROI,threshold);
            ii2 = ii;
            im = ImgSeq(:,:,ii2);
            maskVals = getMaskedValues(im,ROI);
            threshold_s = prctile(maskVals,25);
            if isnan(threshold_s); threshold_s = 0; end
            if isinf(threshold_s); threshold_s = 1; end
            if threshold_s > 1; threshold_s = 1; end
            threshold_s = max(0,threshold_s);
            [Blargest_s, x1, y1] = getlargestROI(ImgSeq(:,:,ii).*ROI,threshold_s);
            nDir = 4*5;
            flowdir = getFlowDir(Blargest_s{1}, x1, y1,BlargestP{1}, xc, yc,nDir);
            ii3 = ii;
            im = ImgSeq(:,:,ii3);
            maskVals = getMaskedValues(im,ROI);
            threshold_s = prctile(maskVals,50);%50
            if isnan(threshold_s); threshold_s = 0; end
            if isinf(threshold_s); threshold_s = 1; end
            if threshold_s > 1; threshold_s = 1; end
            threshold_s = max(0,threshold_s);
            [Blargest, xc, yc] = getlargestROI(ImgSeq(:,:,ii).*ROI,threshold_s);
            [xe,ye] = FlowPoint(Blargest{1},x0,y0,flowdir);
            x0 = xc; y0 = yc;            
            if ii == 1
                TrajInfo.x(ii) = xc;
                TrajInfo.y(ii) = yc;
            end
            if isempty(xe) || isempty(ye)
                xe = nan; ye = nan;
            end
            TrajInfo.x(ii+1) = xe;
            TrajInfo.y(ii+1) = ye;
        end
        idx = find(~isnan(TrajInfo.x));
        if length(idx)>1
            TrajInfo.x = interp1(idx,TrajInfo.x(idx),1:max(idxs)+1);
        end
        idx = find(~isnan(TrajInfo.y));
        if length(idx)>1
            TrajInfo.y = interp1(idx,TrajInfo.y(idx),1:max(idxs)+1);
        end
    end

function TrajInfo = calcAllCenterOfMassTrajectory9(ImgSeq,mask,multiplier,x0,y0,r)
        [dim1,dim2,~] = size(ImgSeq);
        xc = x0; yc = y0;
        ImgSeq = applyMask(ImgSeq,mask,0);
        idxS = 1; idxE = size(ImgSeq,3);
        idxs = idxS:idxE;
        for ii = idxs            
            ROI = CircROI(xc,yc,r,dim1,dim2);%.*mask;
            ii1 = max(ii-1,1);
            im = ImgSeq(:,:,ii1);
            maskVals = getMaskedValues(im,ROI);
            threshold = prctile(maskVals,40);
            if isnan(threshold); threshold = 0; end
            if isinf(threshold); threshold = 1; end
            if threshold > 1; threshold = 1; end
            threshold = max(0,threshold);
            [BlargestP, xc, yc] = getlargestROI(ImgSeq(:,:,ii1).*ROI,threshold);
            ii2 = ii;
            im = ImgSeq(:,:,ii2);
            maskVals = getMaskedValues(im,ROI);
            threshold_s = prctile(maskVals,25);
            if isnan(threshold_s); threshold_s = 0; end
            if isinf(threshold_s); threshold_s = 1; end
            if threshold_s > 1; threshold_s = 1; end
            threshold_s = max(0,threshold_s);
            [Blargest_s, x1, y1] = getlargestROI(ImgSeq(:,:,ii).*ROI,threshold_s);
            nDir = 4*5;
            flowdir = getFlowDir(Blargest_s{1}, x1, y1,BlargestP{1}, xc, yc,nDir);
            ii3 = ii;
            im = ImgSeq(:,:,ii3);
            maskVals = getMaskedValues(im,ROI);
            threshold_s = prctile(maskVals,50);%50
            if isnan(threshold_s); threshold_s = 0; end
            if isinf(threshold_s); threshold_s = 1; end
            if threshold_s > 1; threshold_s = 1; end
            threshold_s = max(0,threshold_s);
            [Blargest, xc, yc] = getlargestROI(ImgSeq(:,:,ii).*ROI,threshold_s);
            [xe,ye] = FlowPoint(Blargest{1},x0,y0,flowdir);
            x0 = xc; y0 = yc;            
            if ii == 1
                TrajInfo.x(ii) = xc;
                TrajInfo.y(ii) = yc;
            end
            if isempty(xe) || isempty(ye)
                xe = nan; ye = nan;
            end
            TrajInfo.x(ii+1) = xe;
            TrajInfo.y(ii+1) = ye;
        end
        idx = find(~isnan(TrajInfo.x));
        if length(idx)>1
            TrajInfo.x = interp1(idx,TrajInfo.x(idx),1:max(idxs)+1);
        end
        idx = find(~isnan(TrajInfo.y));
        if length(idx)>1
            TrajInfo.y = interp1(idx,TrajInfo.y(idx),1:max(idxs)+1);
        end
    end


    function flowdir = getFlowDir(B1,x1,y1,B2,x2,y2,nDir)
        alltheta = linspace(0,2*pi,nDir);
        alldir = cos(alltheta)+1i*sin(alltheta);
        xe1 = zeros(nDir-1,1); ye1 = zeros(nDir-1,1); idx1 = zeros(nDir-1,1);
        xe2 = zeros(nDir-1,1); ye2 = zeros(nDir-1,1); idx2 = zeros(nDir-1,1);
        for k = 1:nDir-1
            [xe1(k),ye1(k),idx1(k)] = FlowPoint(B1,x1,y1,alldir(k));
            [xe2(k),ye2(k),idx2(k)] = FlowPoint(B2,x2,y2,alldir(k));
        end
        [idx1, I1] = sort(idx1);
        [idx2, I2] = sort(idx2);
        
        x1dir = zeros(nDir-1,1); y1dir = zeros(nDir-1,1); x2dir = zeros(nDir-1,1); y2dir = zeros(nDir-1,1);
        for k = 1:nDir-2
            if ~isnan(idx1(k)) && ~isnan(idx1(k+1))
                x1dir(I1(k)) = mean(B1(idx1(k):idx1(k+1),2));
                y1dir(I1(k)) = mean(B1(idx1(k):idx1(k+1),1));
            end
            if ~isnan(idx2(k)) && ~isnan(idx2(k+1))
                x2dir(I2(k)) = mean(B2(idx2(k):idx2(k+1),2));
                y2dir(I2(k)) = mean(B2(idx2(k):idx2(k+1),1));
            end
        end
        for k = nDir-1
            if ~isnan(idx1(k)) && ~isnan(idx1(1))
                x1dir(I1(k)) = mean([B1(1:idx1(1),2); B1(idx1(k):end,2)]);
                y1dir(I1(k)) = mean([B1(1:idx1(1),1); B1(idx1(k):end,1)]);
            end
            if ~isnan(idx2(k)) && ~isnan(idx2(1))
                x2dir(I2(k)) = mean([B2(1:idx2(1),2); B2(idx2(k):end,2)]);
                y2dir(I2(k)) = mean([B2(1:idx2(1),1); B2(idx2(k):end,1)]);
            end
        end
        orig1 = x1 - 1i*y1;
        vec1 = x1dir + 1i*y1dir - orig1;
        orig2 = x2 - 1i*y2;
        vec2 = x2dir + 1i*y2dir - orig2;
        vec = vec1 - vec2;
        flowdir = sum(vec);
        flowdir = flowdir/abs(flowdir);
        
        % plot(B1(:,2),B1(:,1),'w','linewidth',2);
        % plot(x1,y1,'*w');
        % plot(B2(:,2),B2(:,1),'k','linewidth',2);
        % plot(x2,y2,'*k');
        % plot(x1dir,y1dir,'r.')
        % plot(x2dir,y2dir,'g.')
        % quiver(x2,y2,real(flowdir),imag(flowdir),50)
        % quiver(x2dir,y2dir,real(vec),imag(vec),5)
        % n=0;
    end
    function [Blargest, xc, yc] = getlargestROI(im,threshold)
        BW = im2bw(im,threshold);
        CC = bwconncomp(BW,4);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [~,idxM] = max(numPixels);
        if isempty(idxM)
            centerOfMass = [];
            Blargest{1} = [];
            Blargest{2} = [];
        else
            largestROI = zeros(size(BW));
            largestROI(CC.PixelIdxList{idxM}) = 1;
            largestROI = imfill(largestROI,'holes');
            Blargest = bwboundaries(largestROI,4);
            grayImage = im.*largestROI;
            centerOfMass = regionprops(largestROI, grayImage, 'WeightedCentroid');
        end
        if isempty(centerOfMass)
            centerOfMass = [nan,nan];
        else
            centerOfMass = centerOfMass.WeightedCentroid;
        end
        xc = centerOfMass(1); yc = centerOfMass(2);
    end
    function [xe,ye,idx] = FlowPoint(B,xe,ye,avguv)
        if isempty(B)
            idx = nan;
            xe = nan;
            ye = nan;
            return;
        end
        m = imag(avguv)/real(avguv);
        y = m*(B(:,2)-xe)+ye;
        b = B(:,1)-y;
        b = b(1:end-1).*b(2:end);
        idx = find(b<=0);
        b = [B(idx,2), B(idx,1)];
        b = b - repmat([xe,ye],size(b,1),1);
        b = b * [real(avguv); imag(avguv)];
        b = b > 0;
        idx = idx(b);
        B = B(idx,:);
        d = B - repmat([ye,xe],size(B,1),1);
        d = sqrt(sum(d.^2,2));
        [~,idx1] = max(d);
        ye = B(idx1,1);
        xe = B(idx1,2);
        idx = idx(idx1);
        if isempty(idx)
            idx = nan;
            xe = nan;
            ye = nan;
        end
    end

    function pfileName = deletePDF (pdfn,pfn)
        if pdfFlag
            pfileName = sprintf('R1_Artur_Traj_%s_%d.pdf',match{pfn},pdfn);
            if nameExists(pfileName)
                delete(pfileName);
            end
            if nameExists('temp.pdf')
                delete('temp.pdf');
            end
        else
            pfileName = 'temp.pdf';
        end
    end
    function hf = makeFigureWindow (figNum,width,height,magFac)
        hf = figure(figNum);clf; columnWidth = magFac*width; columnHeight = magFac*height;
        set(hf,'Units','inches','Position',[1 1 columnWidth columnHeight],'MenuBar','none','ToolBar','none',...
            'NumberTitle','on','Color','w','Resize','off',...
            'NextPlot','add');
    end
    function [panelLeft, panelBottom, panelWidth, panelHeight] = getPanelProps(nRows,nCols,sprc,whUnits,rtup,rbUnits)
        spaceBetweenRows = sprc(1)*fineSp;
        spaceBetweenColumns = sprc(2) *fineSp;
        panelWidth = 1/nCols + whUnits(1) * fineSp;
        panelHeight = 1/nRows + whUnits(2) * fineSp;
        panelLefts(1) = rtup(1)*fineSp; % shift right for all columns
        for ii = 2:nCols
            panelLefts(ii) = panelLefts(ii-1) + panelWidth + spaceBetweenColumns;
        end
        panelBottoms(1) = rtup(2)*fineSp; % shift up for all rows
        for ii = 2:nRows
            panelBottoms(ii) = panelBottoms(ii-1) + panelHeight + spaceBetweenRows;
        end
        panelBottoms = fliplr(panelBottoms);
        for ii = 1:nCols
            for jj = 1:nRows
                panelLeft(jj,ii) = panelLefts(ii) + rbUnits(1) * fineSp;
                panelBottom(jj,ii) = panelBottoms(jj) + rbUnits(2) * fineSp;
            end
        end
    end
    function paPosition = adjustPosition(paPosition,adjustment)
        paPosition(1) = paPosition(1) + adjustment(1) * fineSp;
        paPosition(2) = paPosition(2) + adjustment(2) * fineSp;
        paPosition(3) = paPosition(3) + adjustment(3) * fineSp;
        paPosition(4) = paPosition(4) + adjustment(4) * fineSp;
    end
    function panelCheck(pos,paPosition)
        ha = axes('Position',paPosition);plot(0,0);xlabel('XLABEL','FontSize',xlabelFontSize);
        ylabel('YLABEL','FontSize',ylabelFontSize);
        titleText = sprintf('TITLE-%d-%d',pos(1),pos(2));
        title(titleText,'FontSize',titleFontSize);
        set(ha,'FontSize',axesFontSize);
    end
    function finalizepdf(pfileName)
        if nameExists('temp.pdf')
            delete('temp.pdf');
        end
        if pdfOpen & pdfFlag
            close(gcf);
            winopen(pfileName);
        end
    end
    function makeTempPDF
        fileName = 'temp.pdf';
        if nameExists(fileName)
            delete(fileName);
        end
        set(gcf, 'PaperPositionMode', 'Auto'); eval(sprintf('print -painters -dpdf -r600 ''%s''',fileName));
        %         set(gcf, 'PaperPositionMode', 'Auto'); eval(sprintf('print -opengl -dpdf -r1200 ''%s''',fileName));
    end
end