function R1_SponPSD_VSD(ss)
clc;close all;
if ~exist('ss','var')
    ss = 1; 
end
stim = {'FL' 'HL' 'VC'};
pdfFileName = sprintf('R1_Spon_PSD_%s.pdf',stim{ss});
pdfFileName = makeName(pdfFileName,getpdffolder);
if nameExists(pdfFileName)
    delete(pdfFileName);
end
if nameExists(makeName('temp.pdf',getpdffolder))
    delete(makeName('temp.pdf',getpdffolder));
end
    
[PSD,f] = calcPSD(ss);
plotMeanPSD(PSD,f,ss);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))


function [PSD,f] = calcPSD(ss)
% load data from base work space where all the data has to be pre loaded
if ss==3
    mdataS = evalin('base','mdataS_VC');
    d = mdataS{1};
else
    mdataS = evalin('base','mdataSDFth');
    d = mdataS{1};
end
stim = {'FL' 'HL' 'VC'}; match = stim{ss};
roi = 1; % the entire imaging window
% PSD estimator parameters
order = 200;
nfft = 2^10;
FsCam = 150;
animalNumber = d.animalNumber;
for an = length(animalNumber):-1:1
    % calc PSD for spon
    Nspon = length(d.d(an,:));
    ns = 0;
    Pxx = [];
    for ii = Nspon:-1:1
        try
            dfMatch = d.d{an,ii}{roi};
            ns = ns+1;
            [Pxx(ns,:),f] = pyulear(dfMatch,order,nfft,FsCam);
        end
    end
    PSD(an,:) = mean(10*log10(Pxx)); % mean power in dB
end
folName = 'Y:\homes\navvab.afrashteh\CloudStation\Navvab-Sam\NeuroPhotonic Paper\Test\First Revision';
fileName = sprintf('%s\\sponPSD_%s.mat',folName,match);
save(fileName,'PSD','f')


function plotMeanPSD(PSD,f,ss)
hf = makeFigureWindow(101,2.2,2,1);
nRows = 1; nCols = 1; sprc = [100 100]; whUnits = [-300 -300]; rtup = [220 220]; rbUnits = [0 0];
[pL, pB, pW, pH] = getPanelProps(nRows,nCols,sprc,whUnits,rtup,rbUnits);
colorSCC = 'k';
N = size(PSD,1);
cSCCmean = mean(PSD,1);
cSCCsem = std(PSD,[],1)/sqrt(N);

% Plot for duration of evoked and spon motifs
paPosition = [pL(1) pB(1) pW pH];
axhist = axes('Position',paPosition);hold on;
axes(axhist); shadedErrorBar(f,cSCCmean,cSCCsem,{'color',colorSCC},0);
plot(f,cSCCmean,'color',colorSCC,'linewidth',1.);

axesFontSize = 8; legendFontSize = 8; xlabelFontSize = 10; ylabelFontSize = 10;
ylabel('Power Spectral Density (dB)','FontSize',ylabelFontSize);
set(axhist,'TickDir','out','FontSize',axesFontSize,'xscale','log');
xlabel('Frequency (Hz)','FontSize',xlabelFontSize); 
Ylim = get(gca,'ylim'); dy = diff(Ylim);
xlim([0,75])
Xlim = get(gca,'xlim'); dx = diff(Xlim);
xticks = [0.5,2,5,10,30,75];
set(gca,'xtick',xticks)
yticks = fliplr(Ylim(2):-20:Ylim(1));
set(gca,'ytick',yticks)
thisCols = {colorSCC};
x1 = 15; x2 = x1+dx*0.1; y1 = Ylim(1):dy/7:Ylim(2); y1 = y1(end-2:end); y2 = y1;
legs = {'Spon'};
for ii = 1:length(legs)
    plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols{ii},'linewidth',2);
    text(x2+0.02*dx,y1(ii),legs{ii},'Color',thisCols{ii},'FontSize',legendFontSize);
end
x = 5; y = y1(end);
if ss==1; txtStr = 'Forelimb';
elseif ss==2; txtStr = 'Hindlimb';
elseif ss==3; txtStr = 'Visual';
end
text(x,y,txtStr,'Color','k','FontSize',legendFontSize+2,'HorizontalAlignment','center');
save2pdf(makeName('temp.pdf',getpdffolder),hf,600);


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



