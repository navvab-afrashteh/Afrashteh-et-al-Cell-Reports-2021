function Artur_fig_bar_graph_speed_directions_AC(ss,roi)

if ~exist('ss','var')
    ss = 1; roi = 2;
end
pdfFileName = sprintf('%s_ss_%d_roi_%d.pdf',mfilename,ss,roi);
pdfFileName = makeName(pdfFileName,getpdffolder)

if exist(pdfFileName,'file')
    delete(pdfFileName);
end
if exist(makeName('temp.pdf',getpdffolder),'file')
    delete(makeName('temp.pdf',getpdffolder));
end

% load data from base work space where all the data has to be pre loaded
mdataE = evalin('base','mdataE_AC');
mdataS = evalin('base','mdataS_AC');
meta = evalin('base','meta');

% data related settings
animals = [1 2 3 4 5 6 7];
selLMH_all = mdataE{ss}.selLMH.selLMH_all;
selLMH = [1 2 3];
% roi = [2];
dfth = 0;
% match = {'FL' 'HL'};
CCTh = 1;
cols = getColors;
selSCC = 3;
amplitude_threshold = 0;%meta.amplitude_threshold;
% mfilename
multLevels = 1;
dataType = 'AC';
for an = 1:length(animals)
    valsE{an} = getDistributions_simpleAmps_artur_2_AC(mdataE{ss},animals(an),roi(ss),ss,amplitude_threshold);
    %     valsS{an} = getDistributionsSpon_simple_artur_glut(mdataS{1},animals(an),roi(ss),ss,dfth,amplitude_threshold);
    for LMH = 3%length(ssLMH)
        valsSCC{an,LMH} = getDistributionsSpon_motifs_simple_artur_2_AC (mdataE,mdataS{1},animals(an),roi(ss),'AC','Max',CCTh,LMH,amplitude_threshold,multLevels,dataType);
    end
    valsS{an} = valsSCC{an,LMH};
end

plotAllCurves(12,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder));

ff = makeFigureWindow__one_axes_only(2,[1 4 1 1.75],[0.38 0.2 0.6 0.69]);
plotBarGraphSpeed (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder));

ff = makeFigureWindow__one_axes_only(3,[3 4 1 1.75],[0.38 0.2 0.6 0.69]);
plotBarGraphDirection (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder));

% save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);

function plotBarGraphSpeed (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;
% 
% for an = 1:length(animals)
%     for iii = 1:length(selLMH_all{an})
%         ii = selLMH_all{an}{iii};
%         avgValsE(an,iii) = mean(valsE{an}.absalluv{ii});
%     end
%     avgValsSCC(an) = mean(valsSCC{an,selSCC}.absalluv);
%     avgValsS(an) = mean(valsS{an}.absalluv);
% end

for an = 1:length(animals)
    for iii = 1:3%length(selLMH_all{ss}(an,:))
%         ii = iii;%selLMH_all{ss}(iii);
        avgValsE(an,iii) = valsE{an}.avgSpeed(iii);
%             avgValsE(an,iii) = abs(valsE{an}.Avgalluv{ii});
    end
%     avgValsSCC(an) = valsSCC{an,selSCC}.avgSpeed;
    avgValsSCC(an) = mean(valsSCC{an,selSCC}.absalluv);
    avgValsS(an) = mean(valsS{an}.absalluv);
end

figure(1000);
post_hoc = 'hsd';
%             annovaVar = [avgValsS' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3) avgValsE(:,4)];
%             [p,tbl,stats] = anova1(annovaVar,{'Spon','Evok_1','Evok_2','Evok_3','Evok_4'},'off');

% annovaVar = [avgValsS' avgValsSCC' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3)];
% [p,tbl,stats] = anova1(annovaVar,{'Spont','SpontCC','Evok_1','Evok_2','Evok_3'},'off');
% thisCols = {cols.sponColorsN{1},cols.eSponColors{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3},cols.stimColors{4}};

annovaVar = [avgValsSCC' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3)];
annovaVar = annovaVar./repmat(mean(annovaVar,2),1,4);
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
ylabel('Speed (mm/sec)','FontSize',ylabelFontSize);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
set(ff.ha,'XTick',1:size(annovaVar,2));
% set(gca,'XTickLabel',{'Spont' 'SpontCC' 'Evok_1' 'Evok_2' 'Evok_3'},'FontSize',axesFontSize-2);
% legs = {'Spont-CC','Lo Evok','Med Evok','Hi Evok'};
% set(gca,'XTickLabel',legs,'FontSize',axesFontSize-2);
set(ff.ha,'XTickLabel',meta.stimNames,'FontSize',axesFontSize-2);
xtickangle(ff.ha,45);
set(ff.ha,'Fontweight','bold','FontSize',axesFontSize);
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);

function plotBarGraphDirection (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;
for an = 1:length(animals)
    for iii = 1:length(selLMH_all{an})
        ii = selLMH_all{an}{iii};
        if animals(an) == 6
            theAng = angle(-ctranspose(mean(valsE{an}.allang{ii}))); theMag = abs(mean(valsE{an}.allang{ii}));
            avgValsE(an,iii) = conj(theMag * exp(1i * theAng));
        else
            avgValsE(an,iii) = conj(mean(valsE{an}.allang{ii}));
        end
    end
    if animals(an) == 6
        theAng = angle(-ctranspose(mean(valsSCC{an,selSCC}.allang))); theMag = abs(mean(valsSCC{an,selSCC}.allang));
        avgValsSCC(an) = conj(theMag * exp(1i * theAng));
        theAng = angle(-ctranspose(mean(valsS{an}.allang))); theMag = abs(mean(valsS{an}.allang));
        avgValsS(an) = conj(theMag * exp(1i * theAng));
    else
        avgValsSCC(an) = conj(mean(valsSCC{an,selSCC}.allang));
        avgValsS(an) = conj(mean(valsS{an}.allang));
    end
end
avgValsE = avgValsE./abs(avgValsE);
avgValsSCC = avgValsSCC./abs(avgValsSCC);
avgValsS = avgValsS./abs(avgValsS);

avgValsE = abs(mean(avgValsE));
avgValsSCC = abs(mean(avgValsSCC));

figure(1000);
post_hoc = 'hsd';

annovaVar = [avgValsSCC' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3)];
[p,tbl,stats] = anova1(annovaVar,{'SpontCC','Evok_1','Evok_2','Evok_3'},'off');
thisCols = {cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3},cols.stimColors{4}};

% [c,~,~,gnames] = multcompare(stats,'CType',post_hoc);
prT = 0;%c(:,6);
hrT = 0;%prT<0.05;
axes(ff.ha);
mVals = annovaVar;
semVals = zeros(size(mVals));
combs = nchoosek(1:size(annovaVar,2),2);

plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols,'ySpacingFactor',10);
xlim([0.4 0.6+size(annovaVar,2)]);
ylim([0 0.7]);
ylabel('A.U.','FontSize',ylabelFontSize);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
set(ff.ha,'XTick',1:size(annovaVar,2));
% set(gca,'XTickLabel',{'Spont' 'SpontCC' 'Evok_1' 'Evok_2' 'Evok_3'},'FontSize',axesFontSize-2);
% legs = {'Spont-CC','Lo Evok','Med Evok','Hi Evok'};
% set(gca,'XTickLabel',legs,'FontSize',axesFontSize-2);
set(ff.ha,'XTickLabel',meta.stimNames,'FontSize',axesFontSize-2);
xtickangle(ff.ha,45);
set(ff.ha,'Fontweight','bold','FontSize',axesFontSize);
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);

function plotAllCurves(figN,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
for an = 1:length(animals)
    lE(an) = length(valsE{an}.allThisdf);
end
gCols = (max(lE)+2);
gRows = length(animals);
totalPlots = gRows * gCols;
grcs = reshape(1:totalPlots,gCols,gRows)';
yl = 10;
amplitude_threshold = 0;%meta.amplitude_threshold;
figure(figN);clf
eaaa = 0;
for an = 1:length(animals)
    if ~isempty(valsSCC{an,selSCC}.allThisdf)
        inds = valsSCC{an,selSCC}.alldf > amplitude_threshold;
        thisPlot = valsSCC{an,selSCC}.allThisuv(inds,:);% valsE{an}.allThisdf{aaa};
        if ~isempty(thisPlot)
            if valsSCC{an,selSCC}.Fs == 220
                [K,N] = size(thisPlot);
                N1 = ceil(N*150/220);
                t = 0:N-1; t = t/220;
                t1 = 0:N1-1; t1 = t1/150;
                temp = [];
                for k=1:K
                    temp(k,:) = interp1(t,thisPlot(k,:),t1);
                end
                thisPlot = temp;
            end
            subplot(gRows,gCols,grcs(an,2));hold on;
            avgPlot = mean(thisPlot);
            allAvgPlotsSCC(an,:) = avgPlot;
            plot(thisPlot');%(:,27:(27+25))');
            plot(avgPlot,'k','linewidth',1.25);
            box off;
            ylim([-0.1 yl]);
            xlim([0 50]);
            title('Spont-CC','FontSize',10);
            plot([11 11],ylim,'-b');
        end
    end
    
    if ~isempty(valsS{an}.allThisdf)
        inds = valsS{an}.alldf > amplitude_threshold;
        thisPlot = valsS{an}.allThisuv(inds,:);
        if ~isempty(thisPlot)
            if valsSCC{an,selSCC}.Fs == 220
                [K,N] = size(thisPlot);
                N1 = ceil(N*150/220);
                t = 0:N-1; t = t/220;
                t1 = 0:N1-1; t1 = t1/150;
                temp = [];
                for k=1:K
                    temp(k,:) = interp1(t,thisPlot(k,:),t1);
                end
                thisPlot = temp;
            end
            subplot(length(animals),gCols,grcs(an,1));hold on;
            avgPlot = mean(thisPlot);
            allAvgPlotsS(an,:) = avgPlot;
            plot(thisPlot');
            plot(avgPlot,'k','linewidth',1.25);
            box off;
            ylim([-0.1 yl]);
            xlim([0 50]);
            title('Spont','FontSize',10);
            plot([11 11],ylim,'-b');
        end
    end
    
    for aaa = 1:length(valsE{an}.allThisdf)
        if ~isempty(valsE{an}.allThisdf{aaa})
            inds = valsE{an}.alldf{aaa} > amplitude_threshold;
            thisPlot = valsE{an}.allThisuv{aaa}(inds,:);
            if ~isempty(thisPlot)
                if valsE{an}.Fs == 220
                    [K,N] = size(thisPlot);
                    N1 = ceil(N*150/220);
                    t = 0:N-1; t = t/220;
                    t1 = 0:N1-1; t1 = t1/150;
                    temp = [];
                    for k=1:K
                        temp(k,:) = interp1(t,thisPlot(k,:),t1);
                    end
                    thisPlot = temp;
                end
                subplot(length(animals),gCols,grcs(an,aaa+2));hold on;
                avgPlot = mean(thisPlot);
                allAvgPlotsE{an,aaa} = avgPlot;
                plot(thisPlot');
                plot(avgPlot,'k','linewidth',1.25);
                box off;
                titleT = sprintf('A = %.3f',mdataE{ss}.amps(an,aaa));
                if aaa == selLMH_all{an}{1}
                    titleT = sprintf('A = %.3f - Lo',mdataE{ss}.amps(an,aaa));
                end
                if aaa == selLMH_all{an}{2}
                    titleT = sprintf('A = %.3f - Med',mdataE{ss}.amps(an,aaa));
                end
                if aaa == selLMH_all{an}{3}
                    titleT = sprintf('A = %.3f - Hi',mdataE{ss}.amps(an,aaa));
                end
                title(titleT,'FontSize',10);
                ylim([-0.1 yl]);
                xlim([0 50]);
                plot([16 16],ylim,'-b');
            end
        end
    end
end

% figure(100);clf;
x_axisS = (0:(size(allAvgPlotsSCC,2)-1)) * 1000/100 - (10000/100);
x_axisE = (0:(size(allAvgPlotsE{1,1},2)-1)) * 1000/150 - (15000/150);
ff = makeFigureWindow__one_axes_only(100,[4 1 1.5 1.75],[0.25 0.2 0.7 0.75]);
thisCols = {cols.sponColorsN{1},cols.eSponColors{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3},cols.stimColors{4}};
hold on;
avgAvgPlotSCC = mean(allAvgPlotsSCC);
plot(x_axisS,avgAvgPlotSCC,'color',thisCols{1});
semAvgPlotSCC = std(allAvgPlotsSCC)/sqrt(length(allAvgPlotsSCC));
shadedErrorBar(x_axisS,avgAvgPlotSCC,semAvgPlotSCC,{'color',thisCols{1}},0.75);
% avgAvgPlotS = mean(allAvgPlotsS);
% plot(x_axis,avgAvgPlotS,'color',thisCols{1});
for aaa = 1:3
    allAvgPlotsE_LMH = [];
    for an = 1:length(animals)
        allAvgPlotsE_LMH(an,:) = allAvgPlotsE{an,selLMH_all{an}{aaa}};
    end
    avgAvgPlotE{aaa} = mean(allAvgPlotsE_LMH);
    plot(x_axisE,avgAvgPlotE{aaa},'color',thisCols{aaa+2});
    semAvgPlotE{aaa} = std(allAvgPlotsE_LMH)/sqrt(length(allAvgPlotsE_LMH));
    shadedErrorBar(x_axisE,avgAvgPlotE{aaa},semAvgPlotE{aaa},{'color',thisCols{aaa+2}},0.75);
end
xlim([-40 max(x_axisE)])
% plot([0 0],ylim,'--m');
set(gca,'TickDir','out','FontSize',6);
set(gca,'linewidth',1.25,'Fontweight','bold','FontSize',6);
xlabel('Time (ms)');
ylabel('Avg Speed (mm/sec)');
thisCols = {cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3}};
x1 = 60; x2 = x1+20; y1 = (3:1:7); y1 = y1(1:4); y2 = y1;
legendFontSize = 6-1;
% legs = {'Spont-CC','Lo Evok','Med Evok','Hi Evok'};
legs = meta.stimNames;
for ii = 1:4
    plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols{ii});
    text(x2+10,y1(ii),legs{ii},'Color',thisCols{ii},'FontSize',legendFontSize);
end
% set(gca,'YScale','log');
ylims = ylim;
ylim([-0.05 15]);
n = 0;
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);
