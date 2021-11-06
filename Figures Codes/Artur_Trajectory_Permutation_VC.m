function Artur_Trajectory_Permutation_VC
% data related settings
match = {'FL' 'HL' 'VC'};
CCTh = 1;
cols = getColors;
selSCC = 3;

[~, name] = system('hostname');
if strcmp(strtrim(name),'imaging1-psy16')  || strcmp(strtrim(name),'imaging4-psy16')
    pathPrefix = 'Y:\homes\navvab.afrashteh\CloudStation\Navvab-Sam\NeuroPhotonic Paper\Test\';
end
if strcmp(strtrim(name),'imaging3-psy16')
    pathPrefix = 'E:\Users\samsoon.inayat\CloudStation\Navvab-Sam\NeuroPhotonic Paper\Test\';
    if ~exist(pathPrefix,'dir')
        pathPrefix = 'Y:\homes\navvab.afrashteh\CloudStation\Navvab-Sam\NeuroPhotonic Paper\Test\';
    end
end
if strcmp(strtrim(name),'Sam-WS')
    pathPrefix = 'T:\CloudStation\Navvab-Sam\NeuroPhotonic Paper\Test\';
end
if strcmp(strtrim(name),'mohajem3-psy13')
    pathPrefix = 'Y:\homes\navvab.afrashteh\CloudStation\Navvab-Sam\NeuroPhotonic Paper\Test\';
end
TrajMethod = 'FE_StretchContour';

ss = 3;
calcPermFlag = 0;
if calcPermFlag
    runPerm(ss)
end
plotPermFlag = 1;
if plotPermFlag
    pdfFlag = 1;
    pdfOpen = 1;
    fineSp = 0.001;
    annotationFontSize = 8; axesFontSize = 6; legendFontSize = 8;
    xlabelFontSize = 9; ylabelFontSize = 9; titleFontSize = 8;
    plotPerm(ss)
end
return;

    function runPerm(ss)
        % load data from base work space where all the data has to be pre loaded
        mdataE = evalin('base','mdataE_VC');
        mdataS = evalin('base','mdataS_VC');
        meta = evalin('base','meta');
        animalNumber = mdataE{ss}.animalNumber;
        selLMH_all = mdataE{ss}.selLMH.selLMH_all;
        animals = [1 2 3 4];
        Cmap = evalin('base','Cmap');
        
        
        % load Trajectory results
        EtrajName = 'Etraj_';
        saveEtraj_FileName = sprintf('%sFigures_new_code\\%s%s_allAn_04corr.mat',pathPrefix ,EtrajName,match{ss});
        load(saveEtraj_FileName,'Etraj')
        
        StrajName = 'Straj_';
        saveStraj_FileName = sprintf('%sFigures_new_code\\%s%s_allAn_04corr.mat',pathPrefix ,StrajName,match{ss});
        load(saveStraj_FileName,'Straj')
        
        % Create axes
        for ann = 1:length(animals)
            an = animals(ann);
            selLMH = 1;           
            anPathName = makeName(animalNumber(an),getMainDataFolder);
            BregmaFileName = makeName('bregmaXY.txt',anPathName);
            bregma = load(BregmaFileName);
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
            LMH = 1;
            [dim1, dim2] = size(mask);
            ne = length(Etraj{an}.Trials{1});
            allXe = Etraj{an}.Xinterp{selLMH(LMH)}(1:ne,:);
            allYe = Etraj{an}.Yinterp{selLMH(LMH)}(1:ne,:);
            p = randperm(ne);
            allXe = allXe(p,:);
            allYe = allYe(p,:);
            % create permutation matrix
            n = 10; k = 500;
            Ke = zeros(k,n);
            i = 0;
            while i<k
                Ke(i+1,:) = sort(randperm(ne,n));
                i = i + 1;
                if i==k
                    Ke = unique(Ke,'rows','stable');
                    i = size(Ke,1);
                    Ke = [Ke;zeros(k-i,n)];
                end
            end
            
            for p = 1:k
                Ic = zeros(dim1,dim2);
                eX = allXe(Ke(p,:),:);
                eY = allYe(Ke(p,:),:);
                imE = zeros(dim1, dim2);
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
                level = graythresh(im);
                imbw = im2bw(im,level);
                
                FD_E(an,p) = BoxCountfracDim(imbw);
            end
            % Spon
            ns = size(Straj{an}.TrajFrames,1);
            allXs = Straj{an}.Xinterp(1:ns,:);
            allYs = Straj{an}.Yinterp(1:ns,:);
            % create permutation matrix
            Ks = zeros(k,n);
            i = 0;
            while i<k
                Ks(i+1,:) = sort(randperm(ns,n));
                i = i + 1;
                if i==k
                    Ks = unique(Ks,'rows','stable');
                    i = size(Ks,1);
                    Ks = [Ks;zeros(k-i,n)];
                end
            end
            % spon perm fractal number
            for p = 1:k
                sX = allXs(Ks(p,:),:);
                sY = allYs(Ks(p,:),:);
                Ic = zeros(dim1,dim2);
                imS = zeros(dim1, dim2);
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
                level = graythresh(im);
                imbw = im2bw(im,level);
                FD_S(an,p) = BoxCountfracDim(imbw);
            end
        end
        FDPermFileName = sprintf('%sFigures_new_code_AC\\Fractional Analysis_Permutation_%s.mat',pathPrefix ,match{ss});
        save(FDPermFileName,'FD_E','FD_S')
        
        FD.FD_E = mean(FD_E,2);
        FD.FD_S = mean(FD_S,2);
        save(sprintf('Fractional Analysis_%s.mat',match{ss}),'FD')
    end

    function plotPerm(ss)
        % suppress warnings
        [msg, id] = lastwarn;
        warning('off', id)
        close all
        nRows = 1; nCols = 5+1; sprc = [10 30]; whUnits = [-34 -300]; 
        rtup = [45 220]; rbUnits = [0 0];
        hf = makeFigureWindow(101,6.9,nRows*1.5,1.5);
        [pL, pB, pW, pH] = getPanelProps(nRows,nCols,sprc,whUnits,rtup,rbUnits);
        % load Trajectory Permutation results
        FDPermFileName = sprintf('%sFigures_new_code_AC\\Fractional Analysis_Permutation_%s.mat',pathPrefix ,match{ss});
        load(FDPermFileName,'FD_E','FD_S')
        animals = 1:size(FD_E,1);
        % calculate distribution and plot for evok and spon
        nbins = 20;
        ninterp = 60;
        minX = 1.30;
        maxX = 1.7;
        Lx = maxX - minX;
        maxY = 15;
        xbins = linspace(minX,maxX,nbins);
        cInterp = linspace(minX,maxX,ninterp);
        linewidth = 2;
        colS = 'b';
        colE = 'r';
        faceAlpha = 0.75;
        for an = 1:length(animals)
            paPositionBoth = [pL(1,an) pB(1,an) pW pH];
            ax1 = axes('parent',hf, 'Position',paPositionBoth);hold on;
            axes(ax1);
            % calc dist
            [cEi(an,:), cSCCi(an,:)] = calcCounts(FD_E(an,:), FD_S(an,:), xbins, cInterp);
            hold on;
            plot(cInterp,cEi(an,:),'col',colE,'linewidth',linewidth);
            plot(cInterp,cSCCi(an,:),'col',colS,'linewidth',linewidth)
            xlim([minX,maxX])
            set(gca,'xtick',[minX, (minX+maxX)/2, maxX])
            ylim([0 maxY])
            box off;
            set(ax1,'tickdir','out')
            tickLen = get(ax1,'ticklength');
            set(ax1,'ticklength',2*tickLen)
            if an == 1
                textString = 'Probability Density Function';
                xt = minX - 0.25*Lx;
                yt = maxY/2;
                text(xt,yt,textString,'fontsize',ylabelFontSize,...
                    'Rotation',90,'FontWeight','bold',...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','middle');
                textString = 'Fractal Dimension';
                xt = minX + 3*Lx;
                yt = -3.5;
                text(xt,yt,textString,'fontsize',xlabelFontSize,...
                    'Rotation',0,'FontWeight','bold',...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','middle');
            end            
            xt = (minX+maxX)/2;
            yt = maxY+1.0;
            text(xt,yt,['Animal #',num2str(an)],...
                'HorizontalAlignment','center','FontWeight','bold',...
                'VerticalAlignment','middle','fontsize',titleFontSize);
            % statistical test
%             [pSR(an),hSR(an)] = kstest2(FD_S(an,:),FD_E(an,:),'Tail','larger');
%             disp(['animal #',num2str(an),': K-S test: h = ', num2str(hSR(an)),', p = ', num2str(pSR(an))])
            [pRS(an),hRS(an),~] = ranksum(FD_S(an,:),FD_E(an,:),'tail','right');
            disp(['animal #',num2str(an),': ranksum test: h = ', num2str(hRS(an)),', p = ', num2str(pRS(an))])
        end
        
        figure(hf);
        paPositionBoth = [pL(1,an+1) pB(1,an+1) pW pH];
        ax1 = axes('parent',hf, 'Position',paPositionBoth);hold on;
        axes(ax1);
        % get avg dist for evok and spon
        meanE = mean(cEi);
        semE = std(cEi)/sqrt(length(animals));
        meanS = mean(cSCCi);
        semS = std(cSCCi)/sqrt(length(animals));
        shadedErrorBar(cInterp,meanE,semE,{'color',colE,'linewidth',linewidth},faceAlpha);
        shadedErrorBar(cInterp,meanS,semS,{'color',colS,'linewidth',linewidth},faceAlpha);
        xlim([minX,maxX])
        set(gca,'xtick',[minX, (minX+maxX)/2, maxX])
        ylim([0 maxY])
        box off;
        set(ax1,'tickdir','out')
        tickLen = get(ax1,'ticklength');
        set(ax1,'ticklength',2*tickLen)
        xt = (minX+maxX)/2;
        yt = maxY+1.0;
        text(xt,yt,'All Animals','HorizontalAlignment','center','FontWeight','bold',...
            'VerticalAlignment','middle','fontsize',titleFontSize);
        
        axes(ax1);
        xt = minX+0.22*Lx;
        yt = 0.9*maxY;
        plot([xt-0.12*Lx,xt-0.02*Lx],yt*[1,1],...
            'col',colS,'linewidth',linewidth)
        text(xt,yt,'Spon','fontsize',legendFontSize)
        yt = 0.8*maxY;
        plot([xt-0.12*Lx,xt-0.02*Lx],yt*[1,1],...
            'col',colE,'linewidth',linewidth)
        text(xt,yt,'Evok','fontsize',legendFontSize)
        % save as pdf
        pdfn = 2*[101 101 101];
        pfileName = deletePDF(pdfn(ss),ss);
        pfileName = makeName(pfileName,getpdffolder);
        save2pdf(pfileName,hf,600);
        % test on means
        mE = mean(FD_E,2);
        mS = mean(FD_S,2);
        [h,p] = ttest2(mS,mE,'tail','right');
        disp(['all animals: ttest: h = ', num2str(h),', p = ', num2str(p)])
        [p,h] = ranksum(mS,mE,'tail','right');
        fprintf('all animals: Wilcoxon rank sum test: h = %d; p = %0.10f\n', h,p)
        fprintf('\n\n');
    end

    function pfileName = deletePDF (pdfn,pfn)
        if pdfFlag
            pfileName = sprintf('Artur_Fractal Dimension PDF_%s_%d.pdf',match{pfn},pdfn);
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
    function [cEi, cSCCi] = calcCounts(bVarE, bVarSCC, xbins, cInterp)
        [cE,centers] = hist(bVarE,xbins);
        [cSCC,~] = hist(bVarSCC,xbins);
        
        cEi = interp1(centers, cE, cInterp);
        cSCCi = interp1(centers, cSCC, cInterp);
        dx = cInterp(2)-cInterp(1);
        cEi = cEi/sum(cEi)/dx;
        cSCCi = cSCCi/sum(cSCCi)/dx;
    end
end

