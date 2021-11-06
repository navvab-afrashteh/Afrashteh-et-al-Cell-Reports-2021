function Artur_fig_bar_graph_amplitude_riseTime_1_VC(ss,cc)
if ~exist('ss','var')
    ss = 3;
end

if ~exist('cc','var')
    cc = 0;
end

pdfFileName = sprintf('%s_%d.pdf',mfilename,ss);
if cc
    pdfFileName = sprintf('onset_%s',pdfFileName);
end
pdfFileName = makeName(pdfFileName,getpdffolder)


if exist(pdfFileName,'file')
    delete(pdfFileName);
end
if exist(makeName('temp.pdf',getpdffolder),'file')
    delete(makeName('temp.pdf',getpdffolder));
end

% load data from base work space where all the data has to be pre loaded
mdataE = evalin('base','mdataE_VC');
mdataS = evalin('base','mdataS_VC');
meta = evalin('base','meta');

% data related settings
animals = [1 2 3 4];
selLMH_all = mdataE{ss}.selLMH.selLMH_all;
roi = 5;
dfth = 2;
match = {'FL' 'HL' 'VC'};
CCTh = 1;
cols = getColors;
selSCC = 3;
amplitude_threshold = meta.amplitude_threshold;
% mfilename

for an = 1:length(animals)
    valsE{an} = getDistributions_simpleAmps_artur_2(mdataE{ss},animals(an),roi,ss,amplitude_threshold);
    valsS{an} = getDistributionsSpon_simple_artur(mdataS{1},animals(an),roi,ss,dfth,amplitude_threshold);
    for LMH = 3%length(ssLMH)
        valsSCC{an,LMH} = getDistributionsSpon_motifs_simple_artur (mdataE,mdataS{1},animals(an),roi,match{ss},'Max',CCTh,LMH,amplitude_threshold,0);
    end
end

FD = load('Fractional Analysis_VC.mat');

if cc
    onset_error = [];
    for an = 1:length(animals)
        for aaa = 1
            thisOnset = valsE{an}.onsets;
            onset_error = [onset_error thisOnset-27];
        end
    end

    onset_error(isnan(onset_error)) = [];
    min_onset = min(onset_error);
    max_onset = max(onset_error);
    bins = min_onset:max_onset;
    ff = makeFigureWindow__one_axes_only(4,[1 1 1.75 2],[0.23 0.25 0.6 0.7]);
    axes(ff.ha);
    [counts,centers] = hist(onset_error,bins);
    counts = 100*counts/sum(counts)
    bar(centers,counts,'FaceColor','b','EdgeColor','b');
    xlim([min_onset max_onset]);
    ylim([0 max(counts)+5]);
    xlabel('Onset Difference (frames)');
    ylabel('Percentage');
    set(gca,'TickDir','out','FontSize',8,'FontWeight','bold');
    save2pdf(pdfFileName,ff.hf,600);
    return;
end

dur = plotAllCurves(12,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))

ff = makeFigureWindow__one_axes_only(2,[1 4 1 1.75],[0.38 0.2 0.6 0.69]);
plotBarGraphFDs (ff,animals,FD,cols,meta)
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))

ff = makeFigureWindow__one_axes_only(2,[1 4 1 1.75],[0.38 0.2 0.6 0.69]);
plotBarGraphRiseTime (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))
% 
ff = makeFigureWindow__one_axes_only(3,[5 4 1 1.75],[0.38 0.2 0.6 0.69]);
plotBarGraphDuration(ff,animals,mdataE,dur,ss,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))
% 
ff = makeFigureWindow__one_axes_only(4,[1 1 1 1.75],[0.39 0.2 0.6 0.68]);
plotBarGraphAmplitude (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))


function plotBarGraphFDs (ff,animals,FD,cols,meta)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;

for an = 1:length(animals)
    for iii = 1
        avgValsE(an,iii) = FD.FD.FD_E(an,iii);
    end
    avgValsSCC(an) = FD.FD.FD_S(an);
end
figure(1000);
annovaVar = [avgValsSCC' avgValsE(:,1)];% avgValsE(:,2) avgValsE(:,3)];
thisCols = {cols.sponColorsN{1},cols.stimLevelColors{3}};
[hrT,prT] = ttest(avgValsSCC',avgValsE);
disp('Fractal Dimension VC paired t-test')
fprintf('pvalue = %0.10f\n\n',prT)
axes(ff.ha);
[mVals semVals] = findMeanAndStandardError(annovaVar,1);
combs = nchoosek(1:size(annovaVar,2),2);

plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols,'ySpacingFactor',10);
xlim([0.4 0.6+size(annovaVar,2)]);
YLim = get(gca,'ylim'); YLim(1) = 1; ylim(YLim)
hy = ylabel('Fractal Dimension (A.U.)');
pos = get(hy,'Position');pos(1) = pos(1);set(hy,'Position',pos);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
set(ff.ha,'XTick',1:size(annovaVar,2));
% set(gca,'XTickLabel',{'Spont' 'SpontCC' 'Evok_1' 'Evok_2' 'Evok_3'},'FontSize',axesFontSize-2);
% legs = {'Spont-CC','Lo Evok','Med Evok','Hi Evok'};
% set(gca,'XTickLabel',legs,'FontSize',axesFontSize-2);
set(ff.ha,'XTickLabel',meta.stimNames);
xtickangle(ff.ha,45);
set(ff.ha,'Fontweight','bold');
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);



function plotBarGraphAmplitude (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;

for an = 1:length(animals)
    avgValsE(an,1) = mean(valsE{an}.alldf{1});
    avgValsSCC(an) = mean(valsSCC{an,selSCC}.alldf);
    avgValsS(an) = mean(valsS{an}.alldf);
end
figure(1000);
post_hoc = 'hsd';
annovaVar = [avgValsSCC' avgValsE(:,1)];
thisCols = {cols.sponColorsN{1},cols.stimLevelColors{3},cols.stimColors{4}};
[hrT,prT] = ttest(avgValsSCC',avgValsE);
disp('Amp VC paired t-test')
fprintf('pvalue = %0.10f\n\n',prT)
axes(ff.ha);
[mVals semVals] = findMeanAndStandardError(annovaVar,1);
combs = nchoosek(1:size(annovaVar,2),2);

plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols,'ySpacingFactor',10);
xlim([0.4 0.6+size(annovaVar,2)]);
ylabel('VSD Signal Amplitude (\DeltaF/F_0 %)','FontSize',ylabelFontSize);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
set(ff.ha,'XTick',1:size(annovaVar,2));
% set(gca,'XTickLabel',{'Spont' 'SpontCC' 'Evok_1' 'Evok_2' 'Evok_3'},'FontSize',axesFontSize-2);
% legs = {'Spont-CC','Lo Evok','Med Evok','Hi Evok'};
% set(gca,'XTickLabel',legs,'FontSize',axesFontSize-2);
set(ff.ha,'XTickLabel',meta.stimNamesVC,'FontSize',axesFontSize-4);
set(ff.ha,'Fontweight','bold','FontSize',axesFontSize);
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);

function plotBarGraphRiseTime (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;

for an = 1:length(animals)
    avgValsE(an,1) = mean(valsE{an}.allRiseTime{1});
    avgValsSCC(an) = mean(valsSCC{an,selSCC}.riseTime);
    avgValsS(an) = mean(valsS{an}.riseTime);
end
figure(1000);
annovaVar = [avgValsSCC' avgValsE(:,1)];
thisCols = {cols.sponColorsN{1},cols.stimLevelColors{3},cols.stimColors{4}};
[hrT,prT] = ttest(avgValsSCC',avgValsE);
disp('PeakTime VC paired t-test')
fprintf('pvalue = %0.10f\n\n',prT)
axes(ff.ha);
[mVals, semVals] = findMeanAndStandardError(annovaVar,1);
combs = nchoosek(1:size(annovaVar,2),2);

plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols,'ySpacingFactor',10);
xlim([0.4 0.6+size(annovaVar,2)]);
ylabel('Time to peak (msec)','FontSize',ylabelFontSize);
set(gca,'TickDir','out','FontSize',axesFontSize);
set(gca,'XTick',1:size(annovaVar,2));
% set(gca,'XTickLabel',{'Spont' 'SpontCC' 'Evok_1' 'Evok_2' 'Evok_3'},'FontSize',axesFontSize-2);
% legs = {'Spont-CC','Lo Evok','Med Evok','Hi Evok'};
% set(gca,'XTickLabel',legs,'FontSize',axesFontSize-2);
set(gca,'XTickLabel',meta.stimNamesVC,'FontSize',axesFontSize-4);
set(gca,'Fontweight','bold','FontSize',axesFontSize);
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);

function plotBarGraphDuration (ff,animals,mdataE,dur,ss,cols,meta)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;
avgValsE = dur.durE;
avgValsSCC = dur.durSCC';
avgValsS = dur.durS';

figure(1000);
annovaVar = [avgValsSCC' avgValsE(:,1)];
thisCols = {cols.sponColorsN{1},cols.stimLevelColors{3},cols.stimColors{4}};
[hrT,prT] = ttest(avgValsSCC',avgValsE);
axes(ff.ha);
[mVals semVals] = findMeanAndStandardError(annovaVar,1);
combs = nchoosek(1:size(annovaVar,2),2);

plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols,'ySpacingFactor',10);
xlim([0.4 0.6+size(annovaVar,2)]);
ylabel('Duration - FWHM (msec)','FontSize',ylabelFontSize);
set(gca,'TickDir','out','FontSize',axesFontSize);
set(gca,'XTick',1:size(annovaVar,2));
% set(gca,'XTickLabel',{'Spont' 'SpontCC' 'Evok_1' 'Evok_2' 'Evok_3'},'FontSize',axesFontSize-2);
% legs = {'Spont-CC','Lo Evok','Med Evok','Hi Evok'};
% set(gca,'XTickLabel',legs,'FontSize',axesFontSize-2);
set(gca,'XTickLabel',meta.stimNamesVC,'FontSize',axesFontSize-4);
set(gca,'Fontweight','bold','FontSize',axesFontSize);
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);


function dur = plotAllCurves(figN,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
% selLMH_all_majid {1} = [1 4 6;1 3 5;2 3 5;1 2 4;1 2 4];
% selLMH_all_majid {2} = [1 2 5;2 3 5;1 2 3;1 3 4;1 2 3];
for an = 1:length(animals)
    lE(an) = length(valsE{an}.allThisdf);
end
gCols = (max(lE)+2);
gRows = length(animals);
totalPlots = gRows * gCols;
grcs = reshape(1:totalPlots,gCols,gRows)';
yl = 1;
amplitude_threshold = 0;%meta.amplitude_threshold;
figure(figN);clf
eaaa = 0;
for an = 1:length(animals)
    if ~isempty(valsSCC{an,selSCC}.allThisdf)
        inds = valsSCC{an,selSCC}.alldf > amplitude_threshold;
        thisPlot = valsSCC{an,selSCC}.allThisdf(inds,:);% valsE{an}.allThisdf{aaa};
        if ~isempty(thisPlot)
            subplot(gRows,gCols,grcs(an,2));hold on;
            avgPlot = mean(thisPlot);
            allAvgPlotsSCC(an,:) = avgPlot;
            plot(thisPlot');%(:,27:(27+25))');
            plot(avgPlot,'k','linewidth',1.25);
            box off;
            ylim([-0.1 yl]);
%             xlim([0 50]);
            title('Spont-CC','FontSize',10);
            plot([16 16],ylim,'-b');
        end
    end

    if ~isempty(valsS{an}.allThisdf)
        inds = valsS{an}.alldf > amplitude_threshold;
        thisPlot = valsS{an}.allThisdf(inds,:);
        if ~isempty(thisPlot)
            subplot(length(animals),gCols,grcs(an,1));hold on;
            avgPlot = mean(thisPlot);
            allAvgPlotsS(an,:) = avgPlot;
            plot(thisPlot');
            plot(avgPlot,'k','linewidth',1.25);
            box off;
            ylim([-0.1 yl]);
%             xlim([0 50]);
            title('Spont','FontSize',10);
            plot([16 16],ylim,'-b');
        end
    end

    for aaa = 1:length(valsE{an}.allThisdf)
        if ~isempty(valsE{an}.allThisdf{aaa})
            inds = valsE{an}.alldf{aaa} > amplitude_threshold;
            thisPlot = valsE{an}.allThisdf{aaa}(inds,:);
            if ~isempty(thisPlot)
                subplot(length(animals),gCols,grcs(an,aaa+2));hold on;
                avgPlot = mean(thisPlot);
                allAvgPlotsE{an,aaa} = avgPlot;
                plot(thisPlot');
                plot(avgPlot,'k','linewidth',1.25);
                box off;
                titleT = sprintf('A = %.3f',mdataE{ss}.amps(an,aaa));
%                 ylim([-0.1 yl]);
%                 xlim([0 50]);
                plot([16 16],ylim,'-b');
            end
        end
    end
end
for an = 1:length(animals)
    thisCurve = allAvgPlotsSCC(an,:);
    durSCC(an,1) = findDuration(thisCurve,16,6.67);
    thisCurve = allAvgPlotsS(an,:);
    durS(an,1) = findDuration(thisCurve,16,6.67);
    for aaa = 1:size(allAvgPlotsE,2)
        thisCurve = allAvgPlotsE{an,aaa};
        durE(an,1) = findDuration(thisCurve,16,6.67);
    end
end
% figure(100);clf;
x_axis = (0:(size(allAvgPlotsSCC,2)-1)) * 1000/150 - (15000/150);
% x_axis = 1:size(allAvgPlotsSCC,2);
ff = makeFigureWindow__one_axes_only(100,[4 1 1.5 1.75],[0.25 0.2 0.7 0.75]);
thisCols = {cols.sponColorsN{1},cols.stimLevelColors{3},cols.stimColors{4}};
hold on;
avgAvgPlotSCC = mean(allAvgPlotsSCC);
plot(x_axis,avgAvgPlotSCC,'color',thisCols{1});
semAvgPlotSCC = std(allAvgPlotsSCC)/sqrt(length(allAvgPlotsSCC));
shadedErrorBar(x_axis,avgAvgPlotSCC,semAvgPlotSCC,{'color',thisCols{1}},0.75);
% avgAvgPlotS = mean(allAvgPlotsS);
% plot(x_axis,avgAvgPlotS,'color',thisCols{1});
for aaa = 1
    allAvgPlotsE_LMH = [];
    for an = 1:length(animals)
        allAvgPlotsE_LMH(an,:) = allAvgPlotsE{an};
    end
    avgAvgPlotE{aaa} = mean(allAvgPlotsE_LMH);
    plot(x_axis,avgAvgPlotE{aaa},'color',thisCols{aaa+1});
    semAvgPlotE{aaa} = std(allAvgPlotsE_LMH)/sqrt(length(allAvgPlotsE_LMH));
    shadedErrorBar(x_axis,avgAvgPlotE{aaa},semAvgPlotE{aaa},{'color',thisCols{aaa+1}},0.75);
end
xlim([-40 max(x_axis)])
% plot([0 0],ylim,'--m');
ylim([-0.05 0.6]);
set(gca,'TickDir','out','FontSize',6);
set(gca,'linewidth',1.25,'Fontweight','bold','FontSize',6);
xlabel('Time (ms)');
ylabel('Avg VSD Signal \DeltaF/F_0 %');
% thisCols = {cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3}};
x1 = 280; x2 = x1+50; y1 = (0.4:0.08:7); y1 = y1(1:4); y2 = y1;
legendFontSize = 6-1;
% legs = {'Spont-CC','Lo Evok','Med Evok','Hi Evok'};
legs = meta.stimNamesVC;
for ii = 1:2
    plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols{ii},'linewidth',2);
    text(x2+10,y1(ii),legs{ii},'Color',thisCols{ii},'FontSize',legendFontSize+2);
end
n = 0;
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);
dur.durSCC = durSCC;
dur.durS = durS;
dur.durE = durE;