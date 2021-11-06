function R1_Artur_fig_bar_graph_amplitude_riseTime_1(ss,cc)
clc;close all;
match = {'FL' 'HL'};
if ~exist('ss','var')
    ss = 2;
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

dur = plotAllCurves(12,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))

ff = makeFigureWindow__one_axes_only(2,[1 4 1 1.75],[0.38 0.2 0.6 0.69]);
plotBarGraphRiseTime (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))

ff = makeFigureWindow__one_axes_only(4,[1 1 1 1.75],[0.39 0.2 0.6 0.68]);
plotBarGraphAmplitude (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))

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
    % spon events excluding PCC events
    sponFrames = valsS{an}.allFrames(:,2:3);
    PCCFrames = valsSCC{an,selSCC}.allFrames(:,2:3);
    Lia = ~ismember(sponFrames,PCCFrames,'rows').';
    amplitude_threshold = 0;%meta.amplitude_threshold;
    inds = valsS{an}.alldf > amplitude_threshold;
    avgValsS(an) = mean(valsS{an}.alldf(and(inds,Lia)));
end
figure(1000);
annovaVar = [avgValsS' avgValsSCC'];
disp('Amp')
thisCols = {'b',cols.sponColorsN{1}};

[prT,tbl,stats] = anova1(annovaVar,{'Spon','SpontCC'},'off');
hrT = prT<0.05;
annovaVar = [avgValsS' avgValsSCC'];% avgValsE(:,1) avgValsE(:,2) avgValsE(:,3)];
[mVals semVals] = findMeanAndStandardError(annovaVar,1);
combs = nchoosek(1:size(annovaVar,2),2);
disp('Amp paired anova')
fprintf('pvalue = %0.10f\n\n',prT)

axes(ff.ha);
plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols,'ySpacingFactor',10);
xlim([0.4 0.6+size(annovaVar,2)]);
hy = ylabel('VSD Signal Amplitude (\DeltaF/F_0 %)');
pos = get(hy,'Position');pos(1) = pos(1) + 0.1;set(hy,'Position',pos);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
set(ff.ha,'XTick',1:size(annovaVar,2));
set(ff.ha,'XTickLabel',{'ex PCC','PCC'});
xtickangle(ff.ha,0);
set(ff.ha,'Fontweight','bold');
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);

function plotBarGraphRiseTime (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;

for an = 1:length(animals)
    avgValsSCC(an) = mean(valsSCC{an,selSCC}.riseTime);    
    % spon events excluding PCC events
    sponFrames = valsS{an}.allFrames(:,2:3);
    PCCFrames = valsSCC{an,selSCC}.allFrames(:,2:3);
    Lia = ~ismember(sponFrames,PCCFrames,'rows').';
    amplitude_threshold = 0;%meta.amplitude_threshold;
    inds = valsS{an}.alldf > amplitude_threshold;
    avgValsS(an) = mean(valsS{an}.riseTime(and(inds,Lia)));
end
figure(1000);
annovaVar = [avgValsS' avgValsSCC'];
disp('RiseTime')
thisCols = {'b',cols.sponColorsN{1}};

[prT,tbl,stats] = anova1(annovaVar,{'Spon','SpontCC'},'off');
hrT = prT<0.05;
[mVals semVals] = findMeanAndStandardError(annovaVar,1);
combs = nchoosek(1:size(annovaVar,2),2);
disp('PeakTime paired anova')
fprintf('pvalue = %0.10f\n\n',prT)

axes(ff.ha);
plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols,'ySpacingFactor',10);
xlim([0.4 0.6+size(annovaVar,2)]);
hy = ylabel('Peak Time (msec)');
pos = get(hy,'Position');pos(1) = pos(1) + 0.0;set(hy,'Position',pos);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
set(ff.ha,'XTick',1:size(annovaVar,2));
set(ff.ha,'yTick',0:40:80);
set(ff.ha,'XTickLabel',{'ex PCC', 'PCC'});
xtickangle(ff.ha,0);
set(ff.ha,'Fontweight','bold');
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);

function plotAllCurves(figN,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
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
        % spon events excluding PCC events
        sponFrames = valsS{an}.allFrames(:,2:3);
        PCCFrames = valsSCC{an,selSCC}.allFrames(:,2:3);
        Lia = ~ismember(sponFrames,PCCFrames,'rows').';
        inds = valsS{an}.alldf > amplitude_threshold;
        thisPlot = valsS{an}.allThisdf(and(inds,Lia),:);
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
thisCols = {cols.eSponColors{1},cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3},cols.stimColors{4}};
thisCols = {'b',cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3},cols.stimColors{4}};
hold on;
avgAvgPlotSCC = mean(allAvgPlotsSCC);
plot(x_axis,avgAvgPlotSCC,'color',thisCols{2});
semAvgPlotSCC = std(allAvgPlotsSCC)/sqrt(length(allAvgPlotsSCC));
shadedErrorBar(x_axis,avgAvgPlotSCC,semAvgPlotSCC,{'color',thisCols{2}},0.75);
avgAvgPlotS = mean(allAvgPlotsS);
plot(x_axis,avgAvgPlotS,'color',thisCols{1});
semAvgPlotS = std(allAvgPlotsS)/sqrt(length(allAvgPlotsS));
shadedErrorBar(x_axis,avgAvgPlotS,semAvgPlotS,{'color',thisCols{1}},0.75);

xlim([min(x_axis) max(x_axis)])
ylim([-0.05 0.6]); ylim([-0.02 0.12]); set(gca,'YTick',[0,0.04,0.08,0.12])
set(gca,'TickDir','out','FontSize',6);
set(gca,'linewidth',1.25,'Fontweight','bold','FontSize',6);
xlabel('Time (ms)');
hy = ylabel('Avg VSD Signal \DeltaF/F_0 %');
pos = get(hy,'Position');pos(1) = pos(1) + 5;set(hy,'Position',pos);
x1 = 240; x2 = x1+40; y1 = (0.35:0.08:7)/5; y1 = y1(1:5); y2 = y1;
legendFontSize = 6-1;
legs = {'ex PCC','PCC','Lo-Ev','Med-Ev','Hi-Ev'};
for ii = 1:2%length(legs)
    plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols{ii},'linewidth',2);
    text(x2+10,y1(ii),legs{ii},'Color',thisCols{ii},'FontSize',legendFontSize+2);
end
if ss == 10
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
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);
