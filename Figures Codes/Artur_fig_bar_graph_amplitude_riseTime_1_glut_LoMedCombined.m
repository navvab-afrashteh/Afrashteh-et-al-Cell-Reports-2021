function Artur_fig_bar_graph_amplitude_riseTime_1_glut_LoMedCombined(ss)
close all
if ~exist('ss','var')
    ss = 1;
end

pdfFileName = sprintf('%s_%d.pdf',mfilename,ss);
pdfFileName = makeName(pdfFileName,getpdffolder)


if exist(pdfFileName,'file')
    delete(pdfFileName);
end
if exist(makeName('temp.pdf',getpdffolder),'file')
    delete(makeName('temp.pdf',getpdffolder));
end

% load data from base work space where all the data has to be pre loaded
mdataE = evalin('base','mdataE_glut_BPF_050Hz_6Hz');
mdataS = evalin('base','mdataS_glut_BPF_050Hz_6Hz');
meta = evalin('base','meta');

% data related settings
animals = [1 2 3 4 5];
selLMH_all = mdataE{ss}.selLMH.selLMH_all;
roi = [2];
dfth = 0;
% match = {'FL' 'HL'};
CCTh = 1;
cols = getColors;
selSCC = 3;
amplitude_threshold = meta.amplitude_threshold;

% for an = 1:length(animals)
%     valsE{an} = getDistributions_simpleAmps_artur_2_glut(mdataE{ss},animals(an),roi(ss),ss,amplitude_threshold);
%     valsS{an} = getDistributionsSpon_simple_artur_glut(mdataS{1},animals(an),roi(ss),ss,dfth,amplitude_threshold);
%     for LMH = 3%length(ssLMH)
%         valsSCC{an,LMH} = getDistributionsSpon_motifs_simple_artur_glut_2 (mdataE,mdataS{1},animals(an),roi(ss),'AC','Max',CCTh,LMH,amplitude_threshold);
%     end
% end
% save('ac_glut_motifs_averaged_data.mat','valsE','valsSCC','valsS');
load('ac_glut_motifs_averaged_data.mat');
load('Y:\homes\navvab.afrashteh\CloudStation\Navvab-Sam\NeuroPhotonic Paper\Test\CorrMaxDF_BPF_05_6Hz_frames_AC_glut_all_04corr.mat');

% dur = plotAllCurves(12,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta,[]);
dur = plotAllCurves(12,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta,CorrMaxDF_frames);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))

ff = makeFigureWindow__one_axes_only(2,[1 4 1 1.75],[0.38 0.2 0.6 0.69]);
plotBarGraphRiseTime (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))

ff = makeFigureWindow__one_axes_only(3,[5 4 1.45 2],[0.23 0.18 0.75 0.75]);
plotBarGraphDuration(ff,animals,mdataE,dur,ss,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))

% ff = makeFigureWindow__one_axes_only(1,[1 1 1.45 2],[0.25 0.18 0.75 0.75]);
ff = makeFigureWindow__one_axes_only(4,[1 1 1 1.75],[0.39 0.2 0.6 0.68]);
plotBarGraphAmplitude (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))

function plotBarGraphAmplitude (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;

for an = 1:length(animals)
    for ii = [1,3]%1:length(selLMH_all{ss}(an,:))
        if ii == 1
            avgValsE(an,ii) = mean([valsE{an}.avg.Amp(ii),valsE{an}.avg.Amp(ii+1)]);
        else
            avgValsE(an,ii) = valsE{an}.avg.Amp(ii);
        end
    end
    avgValsSCC(an) = (valsSCC{an,selSCC}.avg.Amp);
    avgValsS(an) = mean(valsS{an}.alldf);
end
figure(1000);
post_hoc = 'hsd';
%             annovaVar = [avgValsS' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3) avgValsE(:,4)];
%             [p,tbl,stats] = anova1(annovaVar,{'Spon','Evok_1','Evok_2','Evok_3','Evok_4'},'off');

% annovaVar = [avgValsS' avgValsSCC' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3)];
% [p,tbl,stats] = anova1(annovaVar,{'Spont','SpontCC','Evok_1','Evok_2','Evok_3'},'off');
% thisCols = {cols.sponColorsN{1},cols.eSponColors{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3},cols.stimColors{4}};

annovaVar = [avgValsSCC' avgValsE(:,1) avgValsE(:,3)];
% annovaVar = annovaVar./repmat(mean(annovaVar,2),1,size(annovaVar,2));
disp('Amp')
prT = ranovaGlut(annovaVar);
thisCols = {cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{3},cols.stimColors{4}};

% % annovaVar = annovaVar./repmat(mean(annovaVar,2),1,4);
% [p,tbl,stats] = anova1(annovaVar,{'SpontCC','Evok_1','Evok_2'},'off');
% [c,~,~,gnames] = multcompare(stats,'CType',post_hoc);
% prT = c(:,6);

hrT = prT<0.05;
annovaVar = [avgValsSCC' avgValsE(:,1) avgValsE(:,3)];
axes(ff.ha);
[mVals semVals] = findMeanAndStandardError(annovaVar,1);
combs = nchoosek(1:size(annovaVar,2),2);

plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols,'ySpacingFactor',10);
xlim([0.4 0.6+size(annovaVar,2)]);
hy = ylabel('Glu Signal Amplitude (\DeltaF/F_0 %)','FontSize',ylabelFontSize);
pos = get(hy,'Position');pos(1) = pos(1) + 0.2;set(hy,'Position',pos);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
set(ff.ha,'XTick',1:size(annovaVar,2));
% set(gca,'XTickLabel',{'Spont' 'SpontCC' 'Evok_1' 'Evok_2' 'Evok_3'},'FontSize',axesFontSize-2);
% legs = {'Spont-CC','Lo Evok','Med Evok','Hi Evok'};
% set(gca,'XTickLabel',legs,'FontSize',axesFontSize-2);
set(ff.ha,'XTickLabel',meta.stimNames([1,2,4]));
set(ff.ha,'yTick',[0 0.4 0.8 1.2]);
xtickangle(ff.ha,45);
set(ff.ha,'Fontweight','bold','FontSize',axesFontSize);
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);

function plotBarGraphRiseTime (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;

for an = 1:length(animals)
    for ii = [1,3]%1:length(selLMH_all{ss}(an,:))
        avgValsE(an,ii) = mean(valsE{an}.allRiseTime{ii});
        avgValsE(an,ii) = mean(valsE{an}.avg.riseTime(ii));
        if ii == 1
            avgValsE(an,ii) = mean([valsE{an}.allRiseTime{ii},...
                                    valsE{an}.allRiseTime{ii+1}]);
            avgValsE(an,ii) = mean([valsE{an}.avg.riseTime(ii),...
                                    valsE{an}.avg.riseTime(ii+1)]);
        end
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

annovaVar = [avgValsSCC' avgValsE(:,1) avgValsE(:,3)];
disp('RiseTime')
prT = ranovaGlut(annovaVar);
thisCols = {cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{3},cols.stimColors{4}};

% [p,tbl,stats] = anova1(annovaVar,{'SpontCC','Evok_1','Evok_2'},'off');
% [c,~,~,gnames] = multcompare(stats,'CType',post_hoc);
% prT = c(:,6);

hrT = prT<0.05;
annovaVar = [avgValsSCC' avgValsE(:,1) avgValsE(:,3)];
axes(ff.ha);cla
[mVals semVals] = findMeanAndStandardError(annovaVar,1);
combs = nchoosek(1:size(annovaVar,2),2);

plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols,'ySpacingFactor',10);
xlim([0.4 0.6+size(annovaVar,2)]);
hy = ylabel('Peak Time (msec)');
pos = get(hy,'Position');pos(1) = pos(1) + 0.1;set(hy,'Position',pos);
set(gca,'TickDir','out','FontSize',axesFontSize);
set(gca,'XTick',1:size(annovaVar,2));
% set(gca,'XTickLabel',{'Spont' 'SpontCC' 'Evok_1' 'Evok_2' 'Evok_3'},'FontSize',axesFontSize-2);
% legs = {'Spont-CC','Lo Evok','Med Evok','Hi Evok'};
% set(gca,'XTickLabel',legs,'FontSize',axesFontSize-2);
set(gca,'XTickLabel',meta.stimNames([1,2,4]),'FontSize',axesFontSize-4);
xtickangle(ff.ha,45);
set(gca,'Fontweight','bold','FontSize',axesFontSize);
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);

function plotBarGraphDuration (ff,animals,mdataE,dur,ss,cols,meta)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;
avgValsE = dur.durE;
avgValsSCC = dur.durSCC';
% avgValsS = dur.durS';

figure(1000);
post_hoc = 'hsd';
%             annovaVar = [avgValsS' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3) avgValsE(:,4)];
%             [p,tbl,stats] = anova1(annovaVar,{'Spon','Evok_1','Evok_2','Evok_3','Evok_4'},'off');

% annovaVar = [avgValsS' avgValsSCC' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3)];
% [p,tbl,stats] = anova1(annovaVar,{'Spont','SpontCC','Evok_1','Evok_2','Evok_3'},'off');
% thisCols = {cols.sponColorsN{1},cols.eSponColors{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3},cols.stimColors{4}};

annovaVar = [avgValsSCC' avgValsE(:,1) avgValsE(:,2)];
[p,tbl,stats] = anova1(annovaVar,{'SpontCC','Evok_1','Evok_2'},'off');
thisCols = {cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{3},cols.stimColors{4}};

[c,~,~,gnames] = multcompare(stats,'CType',post_hoc);
prT = c(:,6);
hrT = prT<0.05;
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
set(gca,'XTickLabel',meta.stimNames([1,2,4]),'FontSize',axesFontSize-4);
xtickangle(ff.ha,45);
set(gca,'Fontweight','bold','FontSize',axesFontSize);
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);


function dur = plotAllCurves(figN,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta,CorrMaxDF_frames)
selLMH_all_majid {1} = [1 4 6;1 3 5;2 3 5;1 2 4;1 2 4];
selLMH_all_majid {2} = [1 2 5;2 3 5;1 2 3;1 3 4;1 2 3];
for an = 1:length(animals)
    lE(an) = length(valsE{an}.allThisdf);
end
gCols = max(lE);
gRows = length(animals);
totalPlots = gRows * gCols;
grcs = reshape(1:totalPlots,gCols,gRows)';
yl = 2.5;
amplitude_threshold = meta.amplitude_threshold;
hf = figure(figN);clf
eaaa = 0;
for an = 1:length(animals)
    if ~isempty(valsSCC{an,selSCC}.allThisdf)
        inds = valsSCC{an,selSCC}.alldf > amplitude_threshold;
        if ~isempty(CorrMaxDF_frames)
            TrjFrames = CorrMaxDF_frames{an}';
            CorrFrames = valsSCC{an, 3}.allFrames(:,3);
            idx = ismember(CorrFrames,TrjFrames);
            inds = inds & idx';
        end
        thisPlot = valsSCC{an,selSCC}.allThisdf(inds,:);% valsE{an}.allThisdf{aaa};
        if ~isempty(thisPlot)
            subplot(gRows,gCols,grcs(an,1));hold on;
            plot((-10:size(thisPlot,2)-11)/100,thisPlot');
            avgPlot = mean(thisPlot);
            errbar = std(thisPlot)/sqrt(size(thisPlot,1));
            errbar = std(thisPlot);
            errorbar((-10:4:size(avgPlot,2)-11)/100,avgPlot(1:4:end),...
                errbar(1:4:end),'y','linewidth',1);
            plot((-10:size(avgPlot,2)-11)/100,avgPlot,'y','linewidth',5);
            avgPlot = valsSCC{an,selSCC}.avg.df_chopped';
            allAvgPlotsSCC(an,:) = avgPlot;
            plot((-10:size(avgPlot,2)-11)/100,avgPlot,'k','linewidth',2);
            box off;
            ylim([-0.2 yl]);
            xlim([-10 40]/100);
            title('Spont-CC','FontSize',10);
            plot([0 0],ylim,'-b');
        end
    end

    for aaa = [1,3]%1:length(valsE{an}.allThisdf)
        if ~isempty(valsE{an}.allThisdf{aaa})
            inds = valsE{an}.alldf{aaa} > amplitude_threshold;
            thisPlot = valsE{an}.allThisdf{aaa}(inds,:);
            if aaa == 1
                inds = valsE{an}.alldf{aaa+1} > amplitude_threshold;
                thisPlot = [thisPlot; valsE{an}.allThisdf{aaa+1}(inds,:)];
            end
            if ~isempty(thisPlot)
                if aaa == 1
                    subplot(length(animals),gCols,grcs(an,aaa+1));hold on;
                else
                    subplot(length(animals),gCols,grcs(an,aaa));hold on;
                end
                plot((-17:size(thisPlot,2)-18)/150,thisPlot');
                avgPlot = mean(thisPlot);
                plot((-17:size(avgPlot,2)-18)/150,avgPlot,'y','linewidth',5);
                if aaa == 1
                    avgPlot = mean([valsE{an}.avg.df_chopped{aaa:aaa+1}],2)';
                else
                    avgPlot = valsE{an}.avg.df_chopped{aaa}';
                end
                allAvgPlotsE{an,aaa} = avgPlot;
                plot((-17:size(avgPlot,2)-18)/150,avgPlot,'k','linewidth',2);
                box off;
                titleT = sprintf('A = %.3f',mdataE{ss}.amps(an,aaa));
                if aaa == selLMH_all{ss}(an,1)
                    titleT = sprintf('A = %.3f - Lo',mdataE{ss}.amps(an,aaa));
                end
                if aaa == selLMH_all{ss}(an,3)
                    titleT = sprintf('A = %.3f - Hi',mdataE{ss}.amps(an,aaa));
                end
                title(titleT,'FontSize',10);
                if aaa == selLMH_all_majid{ss}(an,1)
                    text(30,0.75,'LoM');
                end
                if aaa == selLMH_all_majid{ss}(an,3)
                    text(30,0.75,'HiM');
                end
                ylim([-0.2 yl]);
                xlim([-17 58]/150);
                plot([0 0],ylim,'-b');
            end
        end
    end
end
save2pdf(makeName('allCurves_df_glut_LoMedCombined.pdf',getpdffolder),hf,600);

for an = 1:length(animals)
    thisCurve = allAvgPlotsSCC(an,:);
    durSCC(an,1) = findDuration(thisCurve,16,6.67);
%     thisCurve = allAvgPlotsS(an,:);
%     durS(an,1) = findDuration(thisCurve,16,6.67);
    for aaa = [1,3]%1:size(allAvgPlotsE,2)
        thisCurve = allAvgPlotsE{an,aaa};
        if aaa == 1
            durE(an,1) = findDuration(thisCurve,16,6.67);
        else
            durE(an,2) = findDuration(thisCurve,16,6.67);
        end
    end
end
% figure(100);clf;
x_axisS = (0:(size(allAvgPlotsSCC,2)-1)) * 1000/100 - (10000/100);
x_axisE = (0:(size(allAvgPlotsE{1,1},2)-1)) * 1000/150 - (15000/150);
ff = makeFigureWindow__one_axes_only(100,[4 1 2 1.75],[0.2 0.2 0.75 0.75]);
thisCols = {cols.sponColorsN{1},cols.eSponColors{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3},cols.stimColors{4}};
hold on;
avgAvgPlotSCC = mean(allAvgPlotsSCC);
plot(x_axisS,avgAvgPlotSCC,'color',thisCols{1});
semAvgPlotSCC = std(allAvgPlotsSCC)/sqrt(length(allAvgPlotsSCC));
shadedErrorBar(x_axisS,avgAvgPlotSCC,semAvgPlotSCC,{'color',thisCols{1}});
% avgAvgPlotS = mean(allAvgPlotsS);
% plot(x_axis,avgAvgPlotS,'color',thisCols{1});
for aaa = [1,3]
    allAvgPlotsE_LMH = [];
    for an = 1:length(animals)
        allAvgPlotsE_LMH(an,:) = allAvgPlotsE{an,selLMH_all{ss}(an,aaa)};
    end
    avgAvgPlotE{aaa} = mean(allAvgPlotsE_LMH);
    plot(x_axisE,avgAvgPlotE{aaa},'color',thisCols{aaa+2});
    semAvgPlotE{aaa} = std(allAvgPlotsE_LMH)/sqrt(length(allAvgPlotsE_LMH));
    shadedErrorBar(x_axisE,avgAvgPlotE{aaa},semAvgPlotE{aaa},{'color',thisCols{aaa+2}});
end
xlim([-40 max([x_axisS x_axisE])])
% plot([0 0],ylim,'--m');
ylim([-0.05 1]);
set(gca,'TickDir','out','FontSize',6);
set(gca,'linewidth',1.25,'Fontweight','bold','FontSize',6);
xlabel('Time (ms)');
hy = ylabel('Avg Glut Signal \DeltaF/F_0 %');
pos = get(hy,'Position');pos(1) = pos(1) + 5;set(hy,'Position',pos);
thisCols = {cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3}};
x1 = 240; x2 = x1+40; y1 = (0.79:0.1:7); y1 = y1(1:4); y2 = y1;
legendFontSize = 6-1;
% legs = {'Spont-CC','Lo Evok','Med Evok','Hi Evok'};
legs = meta.stimNames;
iii = 0;
for ii = [1,2,4]
    iii = iii+1;
    plot([x1 x2],[y1(iii) y2(iii)],'color',thisCols{ii},'linewidth',2);
    text(x2+10,y1(iii),legs{ii},'Color',thisCols{ii},'FontSize',legendFontSize+2);
end
n = 0;
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);
dur.durSCC = durSCC;
dur.durE = durE;

