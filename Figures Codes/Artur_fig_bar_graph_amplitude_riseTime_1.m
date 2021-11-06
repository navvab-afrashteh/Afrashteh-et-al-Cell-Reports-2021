function Artur_fig_bar_graph_amplitude_riseTime_1(ss,cc)
match = {'FL' 'HL'};
if ~exist('ss','var')
    ss = 1;
end

if ~exist('cc','var')
    cc = 0;
end

pdfFileName = sprintf('%s_%s.pdf',mfilename,match{ss});
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
mdataE = evalin('base','mdataE');
mdataS = evalin('base','mdataSDFth');
meta = evalin('base','meta');

% data related settings
animals = [1 2 3 4 5];
selLMH_all = mdataE{ss}.selLMH.selLMH_all;
roi = [2 3];
dfth = 2;

CCTh = 1;
cols = getColors;
selSCC = 3;
amplitude_threshold = meta.amplitude_threshold;
% mfilename
reWrite = 0;
fileName = sprintf('%s_%d.mat',mfilename,ss);
if reWrite
    for an = 1:length(animals)
        valsE{an} = getDistributions_simpleAmps_artur_2(mdataE{ss},animals(an),roi(ss),ss,amplitude_threshold);
        valsS{an} = getDistributionsSpon_simple_artur(mdataS{1},animals(an),roi(ss),ss,dfth,amplitude_threshold);
        for LMH = 3%length(ssLMH)
            valsSCC{an,LMH} = getDistributionsSpon_motifs_simple_artur (mdataE,mdataS{1},animals(an),roi(ss),match{ss},'Max',CCTh,LMH,amplitude_threshold);
        end
    end
    save(fileName,'valsE','valsS','valsSCC');
else
    load(fileName);
end

if cc
    onset_error = [];
    for an = 1:length(animals)
        for aaa = 1:3
            thisOnset = valsE{an}.onsets(selLMH_all{ss}(an,aaa),:);
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
plotBarGraphRiseTime (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))

% % 
ff = makeFigureWindow__one_axes_only(4,[1 1 1 1.75],[0.39 0.2 0.6 0.68]);
plotBarGraphAmplitude (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))


% ff = makeFigureWindow__one_axes_only(3,[5 4 1 1.75],[0.38 0.2 0.6 0.69]);
% plotBarGraphDuration(ff,animals,mdataE,dur,ss,cols,meta);
% append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))

if exist(makeName('temp.pdf',getpdffolder),'file')
    delete(makeName('temp.pdf',getpdffolder));
end

function plotBarGraphAmplitude (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;

for an = 1:length(animals)
    for iii = 1:length(selLMH_all{ss}(an,:))
        ii = selLMH_all{ss}(an,iii);
        avgValsE(an,iii) = mean(valsE{an}.alldf{ii});
    end
    avgValsSCC(an) = mean(valsSCC{an,selSCC}.alldf);
    avgValsS(an) = mean(valsS{an}.alldf);
end
figure(1000);
post_hoc = 'hsd';
%             annovaVar = [avgValsS' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3) avgValsE(:,4)];
%             [p,tbl,stats] = anova1(annovaVar,{'Spon','Evok_1','Evok_2','Evok_3','Evok_4'},'off');

% annovaVar = [avgValsS' avgValsSCC' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3)];
% [p,tbl,stats] = anova1(annovaVar,{'Spont','SpontCC','Evok_1','Evok_2','Evok_3'},'off');
% thisCols = {cols.sponColorsN{1},cols.eSponColors{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3},cols.stimColors{4}};

annovaVar = [avgValsSCC' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3)];
disp('Amp')
prT = ranovaFLHL(annovaVar);
thisCols = {cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3},cols.stimColors{4}};

% [p,tbl,stats] = anova1(annovaVar,{'SpontCC','Evok_1','Evok_2','Evok_3'},'off');
% [c,~,~,gnames] = multcompare(stats,'CType',post_hoc);
% prT = c(:,6);

hrT = prT<0.05;
annovaVar = [avgValsSCC' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3)];
axes(ff.ha);
[mVals semVals] = findMeanAndStandardError(annovaVar,1);
combs = nchoosek(1:size(annovaVar,2),2);

plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols,'ySpacingFactor',10);
xlim([0.4 0.6+size(annovaVar,2)]);
hy = ylabel('VSD Signal Amplitude (\DeltaF/F_0 %)');
pos = get(hy,'Position');pos(1) = pos(1) + 0.4;set(hy,'Position',pos);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
set(ff.ha,'XTick',1:size(annovaVar,2));
% set(gca,'XTickLabel',{'Spont' 'SpontCC' 'Evok_1' 'Evok_2' 'Evok_3'},'FontSize',axesFontSize-2);
% legs = {'Spont-CC','Lo Evok','Med Evok','Hi Evok'};
% set(gca,'XTickLabel',legs,'FontSize',axesFontSize-2);
set(ff.ha,'XTickLabel',meta.stimNames);
xtickangle(ff.ha,45);
set(ff.ha,'Fontweight','bold');
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);

function plotBarGraphRiseTime (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;

for an = 1:length(animals)
    for iii = 1:length(selLMH_all{ss}(an,:))
        ii = selLMH_all{ss}(an,iii);
        avgValsE(an,iii) = mean(valsE{an}.allRiseTime{ii});
    end
    avgValsSCC(an) = mean(valsSCC{an,selSCC}.riseTime);
    avgValsS(an) = mean(valsS{an}.riseTime);
end
figure(1000);
post_hoc = 'hsd';
%             annovaVar = [avgValsS' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3) avgValsE(:,4)];
%             [p,tbl,stats] = anova1(annovaVar,{'Spon','Evok_1','Evok_2','Evok_3','Evok_4'},'off');

% annovaVar = [avgValsS' avgValsSCC' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3)];
% [p,tbl,stats] = anova1(annovaVar,{'Spont','SpontCC','Evok_1','Evok_2','Evok_3'},'off');
% thisCols = {cols.sponColorsN{1},cols.eSponColors{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3},cols.stimColors{4}};

annovaVar = [avgValsSCC' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3)];
disp('RiseTime')
prT = ranovaFLHL(annovaVar);
thisCols = {cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3},cols.stimColors{4}};

% [p,tbl,stats] = anova1(annovaVar,{'SpontCC','Evok_1','Evok_2','Evok_3'},'off');
% [c,~,~,gnames] = multcompare(stats,'CType',post_hoc);
% prT = c(:,6);

hrT = prT<0.05;
annovaVar = [avgValsSCC' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3)];
axes(ff.ha);
[mVals semVals] = findMeanAndStandardError(annovaVar,1);
combs = nchoosek(1:size(annovaVar,2),2);

plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols,'ySpacingFactor',10);
xlim([0.4 0.6+size(annovaVar,2)]);
hy = ylabel('Peak Time (msec)');
pos = get(hy,'Position');pos(1) = pos(1) + 0.5;set(hy,'Position',pos);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
set(ff.ha,'XTick',1:size(annovaVar,2));
set(ff.ha,'yTick',0:40:80);
% set(gca,'XTickLabel',{'Spont' 'SpontCC' 'Evok_1' 'Evok_2' 'Evok_3'},'FontSize',axesFontSize-2);
% legs = {'Spont-CC','Lo Evok','Med Evok','Hi Evok'};
% set(gca,'XTickLabel',legs,'FontSize',axesFontSize-2);
set(ff.ha,'XTickLabel',meta.stimNames);
xtickangle(ff.ha,45);
set(ff.ha,'Fontweight','bold');
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);

function plotBarGraphDuration (ff,animals,mdataE,dur,ss,cols,meta)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;
avgValsE = dur.durE;
avgValsSCC = dur.durSCC';
avgValsS = dur.durS';

figure(1000);
post_hoc = 'hsd';
%             annovaVar = [avgValsS' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3) avgValsE(:,4)];
%             [p,tbl,stats] = anova1(annovaVar,{'Spon','Evok_1','Evok_2','Evok_3','Evok_4'},'off');

% annovaVar = [avgValsS' avgValsSCC' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3)];
% [p,tbl,stats] = anova1(annovaVar,{'Spont','SpontCC','Evok_1','Evok_2','Evok_3'},'off');
% thisCols = {cols.sponColorsN{1},cols.eSponColors{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3},cols.stimColors{4}};

annovaVar = [avgValsSCC' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3)];
[p,tbl,stats] = anova1(annovaVar,{'SpontCC','Evok_1','Evok_2','Evok_3'},'off');
thisCols = {cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3},cols.stimColors{4}};

[c,~,~,gnames] = multcompare(stats,'CType',post_hoc);
prT = c(:,6);
hrT = prT<0.05;
axes(ff.ha);
[mVals semVals] = findMeanAndStandardError(annovaVar,1);
combs = nchoosek(1:size(annovaVar,2),2);

plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols,'ySpacingFactor',10);
xlim([0.4 0.6+size(annovaVar,2)]);
hy = ylabel('Duration - FWHM (msec)');
pos = get(hy,'Position');pos(1) = pos(1) + 0.9;set(hy,'Position',pos);
set(gca,'TickDir','out','FontSize',axesFontSize);
set(gca,'XTick',1:size(annovaVar,2));
% set(gca,'XTickLabel',{'Spont' 'SpontCC' 'Evok_1' 'Evok_2' 'Evok_3'},'FontSize',axesFontSize-2);
% legs = {'Spont-CC','Lo Evok','Med Evok','Hi Evok'};
% set(gca,'XTickLabel',legs,'FontSize',axesFontSize-2);
set(gca,'XTickLabel',meta.stimNames);
xtickangle(ff.ha,45);
set(gca,'Fontweight','bold');
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);


function dur = plotAllCurves(figN,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
selLMH_all_majid {1} = [1 4 6;1 3 5;2 3 5;1 2 4;1 2 4];
selLMH_all_majid {2} = [1 2 5;2 3 5;1 2 3;1 3 4;1 2 3];
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
                if aaa == selLMH_all{ss}(an,1)
                    titleT = sprintf('A = %.3f - Lo',mdataE{ss}.amps(an,aaa));
                end
                if aaa == selLMH_all{ss}(an,2)
                    titleT = sprintf('A = %.3f - Med',mdataE{ss}.amps(an,aaa));
                end
                if aaa == selLMH_all{ss}(an,3)
                    titleT = sprintf('A = %.3f - Hi',mdataE{ss}.amps(an,aaa));
                end
                title(titleT,'FontSize',10);
                if aaa == selLMH_all_majid{ss}(an,1)
                    text(30,0.75,'LoM');
                end
                if aaa == selLMH_all_majid{ss}(an,2)
                    text(30,0.75,'MedM');
                end
                if aaa == selLMH_all_majid{ss}(an,3)
                    text(30,0.75,'HiM');
                end
                ylim([-0.1 yl]);
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
        if aaa == selLMH_all{ss}(an,1)
            thisCurve = allAvgPlotsE{an,aaa};
            durE(an,1) = findDuration(thisCurve,16,6.67);
        end
        if aaa == selLMH_all{ss}(an,2)
            thisCurve = allAvgPlotsE{an,aaa};
            durE(an,2) = findDuration(thisCurve,16,6.67);
        end
        if aaa == selLMH_all{ss}(an,3)
            thisCurve = allAvgPlotsE{an,aaa};
            durE(an,3) = findDuration(thisCurve,16,6.67);
        end
    end
end
% figure(100);clf;
x_axis = (0:(size(allAvgPlotsSCC,2)-1)) * 1000/150 - (15000/150);
% x_axis = 1:size(allAvgPlotsSCC,2);
ff = makeFigureWindow__one_axes_only(100,[4 1 2 1.75],[0.2 0.2 0.75 0.75]);
thisCols = {cols.sponColorsN{1},cols.eSponColors{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3},cols.stimColors{4}};
hold on;
avgAvgPlotSCC = mean(allAvgPlotsSCC);
plot(x_axis,avgAvgPlotSCC,'color',thisCols{1});
semAvgPlotSCC = std(allAvgPlotsSCC)/sqrt(length(allAvgPlotsSCC));
shadedErrorBar(x_axis,avgAvgPlotSCC,semAvgPlotSCC,{'color',thisCols{1}},0.75);
% avgAvgPlotS = mean(allAvgPlotsS);
% plot(x_axis,avgAvgPlotS,'color',thisCols{1});
for aaa = 1:3
    allAvgPlotsE_LMH = [];
    for an = 1:length(animals)
        allAvgPlotsE_LMH(an,:) = allAvgPlotsE{an,selLMH_all{ss}(an,aaa)};
    end
    avgAvgPlotE{aaa} = mean(allAvgPlotsE_LMH);
    plot(x_axis,avgAvgPlotE{aaa},'color',thisCols{aaa+2});
    semAvgPlotE{aaa} = std(allAvgPlotsE_LMH)/sqrt(length(allAvgPlotsE_LMH));
    shadedErrorBar(x_axis,avgAvgPlotE{aaa},semAvgPlotE{aaa},{'color',thisCols{aaa+2}});
end
xlim([min(x_axis) max(x_axis)])
% plot([0 0],ylim,'--m');
ylim([-0.05 0.6]);
set(gca,'TickDir','out','FontSize',6);
set(gca,'linewidth',1.25,'Fontweight','bold','FontSize',6);
xlabel('Time (ms)');
hy = ylabel('Avg VSD Signal \DeltaF/F_0 %');
pos = get(hy,'Position');pos(1) = pos(1) + 10;set(hy,'Position',pos);
thisCols = {cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3}};
x1 = 240; x2 = x1+40; y1 = (0.35:0.06:7); y1 = y1(1:4); y2 = y1;
legendFontSize = 6-1;
% legs = {'Spont-CC','Lo Evok','Med Evok','Hi Evok'};
legs = meta.stimNames;
for ii = 1:4
    plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols{ii},'linewidth',2);
    text(x2+10,y1(ii),legs{ii},'Color',thisCols{ii},'FontSize',legendFontSize+2);
end
if ss == 1
    ahaw = 15;
    ah = 0.03;
    alw = 0.5;
    maxY = max(avgAvgPlotE{3}) + max(semAvgPlotE{3})+0.001;
    plot([-70 70],[maxY maxY],'color','k','linewidth',alw);
    plot([-30 -30],[0.01 maxY-0.01],'color','k','linewidth',alw);
    plot([-30 -30-ahaw],[maxY-0.01 maxY-0.01-ah],'color','k','linewidth',alw); 
    plot([-30 -30+ahaw],[maxY-0.01 maxY-0.01-ah],'color','k','linewidth',alw);
    plot([-30 -30-ahaw],[0.01 0.01+ah],'color','k','linewidth',alw); 
    plot([-30 -30+ahaw],[0.01 0.01+ah],'color','k','linewidth',alw);
    text(-55,maxY/3,'Amplitude','FontSize',7,'color','k','rotation',90);
end
n = 0;
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);
dur.durSCC = durSCC;
dur.durS = durS;
dur.durE = durE;