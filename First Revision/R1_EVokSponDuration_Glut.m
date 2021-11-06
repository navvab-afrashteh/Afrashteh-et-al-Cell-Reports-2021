function R1_EVokSponDuration_Glut(ss)
clc; close all;
if ~exist('ss','var')
    ss = 1;
end
stim = {'AC'};
pdfFileName = sprintf('R1_Evok_High_Spon_Duration_Glut_%s.pdf',stim{ss});
pdfFileName = makeName(pdfFileName,getpdffolder);
if nameExists(pdfFileName)
    delete(pdfFileName);
end
if nameExists(makeName('temp.pdf',getpdffolder))
    delete(makeName('temp.pdf',getpdffolder));
end
    
[m_evok,m_spon,cE_all,cSCC_all,cInterp] = calcDur(ss);
plotMedianDur(m_evok,m_spon,ss);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))
plotDistDur(cE_all,cSCC_all,cInterp,ss);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))


function [m_evok,m_spon,cE_all,cSCC_all,cInterp] = calcDur(ss)
% load data from base work space where all the data has to be pre loaded
mdataE = evalin('base','mdataE_glut_BPF_050Hz_6Hz');
mdataS = evalin('base','mdataS_glut_BPF_050Hz_6Hz');
            
mainDataFolder = getMainDataFolderGlut;
multLevels = 1;

d = mdataS{1};
thn = 1;
stim = {'AC'}; match = stim{ss};
roi = 1;
CrossOrMax = 'Max';
sf1 = 25; sf2 = 25;
xM = sf1+sf2+1;
x = 1:xM; m_interp = 10; % interpolation
xq = linspace(1,xM,m_interp*xM);

nbins = 10;
ninterp = nbins*2;
minX = 0;
maxX = (xM+1); maxX = maxX*1000/100;% frames to msec
xbins = linspace(minX,maxX,nbins); xbins = [0:1:xbins(2)/2, xbins(2:end-1),round(xbins(end-1)+xbins(2)/10):5:round(xbins(end)+xbins(2)/2)];
cInterp = linspace(minX,maxX,ninterp); cInterp = [cInterp(1:end-1),round(cInterp(end-1)+cInterp(2)/8):5:xbins(end)];

animalNumber = d.animalNumber; animalNumber = animalNumber(1:5);
for an = 1:length(animalNumber)
    anPathName = makeName(animalNumber(an),mainDataFolder);
    ePathName = makeName('pEvoked',anPathName);
    [~, names]= getEvokedListsDiffStimAmps_glut(animalNumber(an),multLevels);
    listR = listOfRecordingsSpon_glut(animalNumber(an));
    selLMH = 3;
    LMHname = names{4}{mdataE{ss}.anGroup{an}{selLMH}};
    fPathName = makeName(LMHname,ePathName);
    eTrajInfoPathName = makeName('TrajInfo1.mat',fPathName);
    load(eTrajInfoPathName,'TrajInfo')
    % calc sDT1 and sDT2 for evoked
    Edt = [];
    jjj=1;
    for jj = 1:length(TrajInfo.TrialsTrajInfo)
        try
            Edt(jjj) =  TrajInfo.TrialsTrajInfo{jj}(2,3) - TrajInfo.TrialsTrajInfo{jj}(1,3)+1;
            jjj = jjj+1;
        catch
        end
    end
    idxNaN = isnan(Edt);
    idxNaN = idxNaN > 0;
    Edt(idxNaN) = [];
    % calc DT1 and DT2 for spon
    P1 = TrajInfo.TrialsTrajInfo_Stat(1)/100;
    P2 = TrajInfo.TrialsTrajInfo_Stat(2)/100;
    P1 = 0.25; P2 = 0.25;
    ns = 0;
    dt1 = []; dt2 = [];
    durAllSpon(an) = 0;
    for ii = 1%:length(listR)
        motifs = getMotifsFromCCResultsSam_simple_glut(animalNumber(an),listR,ii,thn,match,CrossOrMax,LMHname,mdataE{1},multLevels);
        if ~isfield(motifs, 'frames')
            continue;
        end
        dfMatch = d.d{an,ii}{roi};
        uvMatch = abs(d.d1{an,ii}{roi});
        motifs = refineMotifs(motifs,dfMatch,uvMatch);
        frames = motifs.framesdf;
        durAllSpon(an) = durAllSpon(an) + length(dfMatch);
        for jj = 1:length(frames)
            f = frames(jj);
            try
                mImgSeq = dfMatch(f-sf1:f+sf2);
                mImgSeq = interp1(x,mImgSeq,xq,'spline');
                [maxVal,fm] = max(mImgSeq);
                fm = round((sf1+0.5)*m_interp);
                maxVal = mImgSeq(fm);
                [baseline,idxB] = min(mImgSeq(1:fm));
                valS = baseline + P1*(maxVal-baseline);
                temp = mImgSeq(idxB:fm)- valS;
                idxS = find(temp <= 0, 1, 'last');
                idxS = idxS + idxB; % start index
                valE = baseline + P2*(maxVal-baseline);
                temp = mImgSeq(fm:end) - valE;
                idxE = find(temp <= 0, 1 );
                idxE = idxE + fm - 2; % end index
                if isempty(idxE)
                    idxE = length(mImgSeq);
                end
                if fm>idxS && idxE>fm
                    ns = ns+1;
                    dt1(ns) = (fm - idxS)/m_interp;
                    dt2(ns) = (idxE - fm)/m_interp;
                end
            end
        end
    end
    Sdt = dt1+dt2;
    sponDurPercentage(an) = sum(Sdt)/durAllSpon(an);
    % calculate the histogram
    bVarE = Edt; bVarE = bVarE*1000/150; % sampling rate for evoked is 150
    bVarSCC = Sdt; bVarSCC = bVarSCC*1000/100; % sampling rate for evoked is 100
    [cEi, cSCCi] = calcCounts(bVarE, bVarSCC, xbins, cInterp);
    cE_all(an,:) = cEi;
    cSCC_all(an,:) = cSCCi;
    % stats
    [pan(an),han(an),stat] = ranksum(bVarE,bVarSCC,'tail','left');
    m_evok(an) = median(bVarE);
    m_spon(an) = median(bVarSCC);
end
sponDurPercentage = sponDurPercentage*100; % percentage
[mean(m_evok) mean(m_spon) mean(sponDurPercentage)]
folName = 'Y:\homes\navvab.afrashteh\CloudStation\Navvab-Sam\NeuroPhotonic Paper\Test\First Revision';
fileName = sprintf('%s\\sponDurPercentage_%s',folName,match);
save(fileName,'sponDurPercentage')


function plotDistDur(cE_all,cSCC_all,cInterp,ss)
hf = makeFigureWindow(101,2.2,2,1);
nRows = 1; nCols = 1; sprc = [100 100]; whUnits = [-300 -300]; rtup = [220 220]; rbUnits = [0 0];
[pL, pB, pW, pH] = getPanelProps(nRows,nCols,sprc,whUnits,rtup,rbUnits);

colorE = 'r';
colorSCC = 'b';
N = size(cE_all,1);
cEmean = mean(cE_all,1);
cSCCmean = mean(cSCC_all,1);
cEsem = std(cE_all,[],1)/sqrt(N);
cSCCsem = std(cSCC_all,[],1)/sqrt(N);

% Plot for duration of evoked and spon motifs
paPosition = [pL(1) pB(1) pW pH];
axhist = axes('Position',paPosition);hold on;
axes(axhist); shadedErrorBar(cInterp,cSCCmean,cSCCsem,{'color',colorSCC},0);
% axes(axhist); shadedErrorBar(cInterp,cEmean,cEsem,{'color',colorE},0);
plot(cInterp,cSCCmean,'color',colorSCC,'linewidth',1.);
% plot(cInterp,cEmean,'color',colorE,'linewidth',1.);

xlim([0 max(cInterp)]);
ylim([0 8*10^(-3)]);
axesFontSize = 8; legendFontSize = 8; xlabelFontSize = 10; ylabelFontSize = 10;
hy = ylabel('Probability Density Function','FontSize',ylabelFontSize);
pos = get(hy,'Position');pos(1) = pos(1)-10;set(hy,'Position',pos);
set(axhist,'TickDir','out','FontSize',axesFontSize);
xlabel('Duration (msec)','FontSize',xlabelFontSize); 
Ylim = get(gca,'ylim'); dy = diff(Ylim);
Xlim = get(gca,'xlim'); dx = diff(Xlim);
xticks = 0:max(cInterp)/5:max(cInterp); xticks = 10*round(0.1*xticks);
xticks = linspace(0,510,4); xticks = 10*round(0.1*xticks);
set(gca,'xtick',xticks)
yticks = linspace(Ylim(1),Ylim(2),3);
set(gca,'ytick',yticks)
thisCols = {colorSCC,colorE};
x1 = Xlim(1)+dx*0.8; x2 = x1+dx*0.1; y1 = Ylim(1):dy/7:Ylim(2); y1 = y1(end-2:end); y2 = y1;
legs = {'Spon','Evok'};
if ss==10
    for ii = 1%:length(legs)
        plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols{ii},'linewidth',2);
        text(x2+0.005*dx,y1(ii),legs{ii},'Color',thisCols{ii},'FontSize',legendFontSize);
    end
end
x = Xlim(1)+dx*0.5; y = y1(end);
txtStr = 'Auditory';
text(x,y,txtStr,'Color','k','FontSize',legendFontSize+2,'HorizontalAlignment','center');

save2pdf(makeName('temp.pdf',getpdffolder),hf,600);


function plotMedianDur(m_evok,m_spon,ss)
hf = makeFigureWindow(101,1,1.5,1);
nRows = 1; nCols = 1; sprc = [100 100]; whUnits = [-420 -300]; rtup = [360 220]; rbUnits = [0 0];
[pL, pB, pW, pH] = getPanelProps(nRows,nCols,sprc,whUnits,rtup,rbUnits);
paPosition = [pL pB pW pH];
axes('Position',paPosition);hold on;
hold on;
axesFontSize = 10; legendFontSize = 10; 
cols = getColors();

all_m_evok = mean(m_evok);
all_m_spon = mean(m_spon);
all_sem_evok = std(m_evok)/sqrt(length(m_evok));
all_sem_spon = std(m_spon)/sqrt(length(m_spon));
xse = ones(size(m_evok));
xss = ones(size(m_spon))*2;
for ii = 1:length(m_evok)
    plot(xse(ii),m_evok(ii),'.','color',cols.avgStimColors{1});
    plot(xss(ii),m_spon(ii),'.','color',cols.sponColors{1});
    plot([xse(ii) xss(ii)],[m_evok(ii) m_spon(ii)],'color','k');
end
thiscol = {cols.avgStimColors{1} cols.sponColors{1}};
xsem = xse(1)-0.25;
xssm = xss(1)+0.25;
plot(xsem,all_m_evok,'.','color',thiscol{1});
errorbar(xsem,all_m_evok(1),all_sem_evok(1),all_sem_evok(1),'color',thiscol{1});
plot(xssm,all_m_spon,'.','color',thiscol{2});
errorbar(xssm,all_m_spon(1),all_sem_spon(1),all_sem_spon(1),'color',thiscol{2});
xlim([0.5 2.5]);
[hrT,prT] = ttest(m_evok,m_spon,'Tail','left')
ylabel('Median Duration (ms)','FontSize',axesFontSize);
if ss == 1
   xlabel('Forelimb','FontSize',axesFontSize);
end
if ss == 2
   xlabel('Hindlimb','FontSize',axesFontSize);
end
if ss == 3
   xlabel('Visual','FontSize',axesFontSize);
end 
set(gca,'TickDir','out','FontSize',axesFontSize-2);
xlab = {'E','S'};
set(gca,'XTick',[1 2],'XTickLabel',xlab);

xl = xlim; yl = ylim; dxl = xl(2) - xl(1); dyl = yl(2) - yl(1);
if prT<0.001
    text(xl(1)+0.5*dxl,yl(2)-0*dyl/8,'***','fontsize',legendFontSize+1,'FontWeight','normal')
elseif prT<0.01
    text(xl(1)+0.5*dxl,yl(2)-0*dyl/8,'**','fontsize',legendFontSize+1,'FontWeight','normal')
elseif prT<0.05
    text(xl(1)+0.5*dxl,yl(2)-0*dyl/8,'*','fontsize',legendFontSize+1,'FontWeight','normal')
end
save2pdf(makeName('temp.pdf',getpdffolder),hf,600);


function [cEi, cSCCi] = calcCounts(bVarE, bVarSCC, xbins, cInterp)
[cE,centers] = hist(bVarE,xbins);
[cSCC,~] = hist(bVarSCC,xbins);
cEi = interp1(centers, cE, cInterp);
cSCCi = interp1(centers, cSCC, cInterp);
dx = cInterp(2)-cInterp(1);
cEi = cEi/sum(cEi)/dx; cEi(1) = 0;
cSCCi = cSCCi/sum(cSCCi)/dx; cSCCi(1) = 0;


function hf = makeFigureWindow (figNum,width,height,magFac)
hf = figure(figNum);clf; columnWidth = magFac*width; columnHeight = magFac*height;
set(hf,'Units','inches','Position',[1 1 columnWidth columnHeight],'MenuBar','none','ToolBar','none',...
    'NumberTitle','on','Color','w','Resize','off',...
    'NextPlot','add');


function [panelLeft, panelBottom, panelWidth, panelHeight] = getPanelProps(nRows,nCols,sprc,whUnits,rtup,rbUnits)
fineSp = 0.001;
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



