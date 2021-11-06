function Artur_fig_bar_graph_speed_directions_2_AC_average(ss,roi)

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
animals = [1 2 3 4 5];
selLMH_all = mdataE{ss}.selLMH.selLMH_all;
selLMH = [1 2 3];
% roi = [2];
dfth = 0;
% match = {'FL' 'HL'};
CCTh = 1;
cols = getColors;
selSCC = 3;
amplitude_threshold = meta.amplitude_threshold;
% mfilename

multLevels = 1;
dataType = 'AC';

% for an = 1:length(animals)
%     valsE{an} = getDistributions_simpleAmps_artur_2_AC(mdataE{ss},animals(an),roi(ss),ss,amplitude_threshold);
% %     valsS{an} = getDistributionsSpon_simple_artur_glut(mdataS{1},animals(an),roi(ss),ss,dfth,amplitude_threshold);
%     for LMH = 3%length(ssLMH)
%         valsSCC{an,LMH} = getDistributionsSpon_motifs_simple_artur_2_AC (mdataE,mdataS{1},animals(an),roi(ss),'AC','Max',CCTh,LMH,amplitude_threshold,multLevels,dataType);
%     end
%     valsS{an} = valsSCC{an,LMH};
% end
% save('ac_motifs_averaged_data.mat','valsE','valsS','valsSCC');
load('ac_motifs_averaged_data.mat');


for ann = 1:length(animals)
    an = animals(ann);
    for iii = 1:3
        ii = iii;
        avgValsE(ann,iii) = abs(valsE{an}.avg.Speed(ii));

    end
    avgValsSCC(ann) = valsSCC{an,selSCC}.avg.Speed;
end
annovaVar = [avgValsSCC' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3)];
meta.mES = mean(annovaVar,2);

% valS(6) = [];
plotAllCurves(12,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder));

ff = makeFigureWindow__one_axes_only(2,[1 4 1 1.75],[0.38 0.2 0.6 0.69]);
plotBarGraphSpeed (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder));
% 
% ff = makeFigureWindow__one_axes_only(3,[3 4 1 1.75],[0.38 0.2 0.6 0.69]);
% plotBarGraphDirection (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
% append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder));

% save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);

function plotBarGraphSpeed (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;

for ann = 1:length(animals)
    an = animals(ann);
    for iii = 1:3
        ii = iii;
        avgValsE(ann,iii) = abs(valsE{an}.avg.Speed(ii));
%         avgValsE(ann,iii) = mean(valsE{an}.absalluv{ii});
%         try
%             avgValsE(ann,iii) = abs(valsE{an}.Avgalluv{ii});
%         catch
%             avgValsE(an,iii) = 0;
%         end
    end
    avgValsSCC(ann) = valsSCC{an,selSCC}.avg.Speed;%valsSCC{an,selSCC}.avgSpeed;
%     avgValsSCC(an) = valsSCC{an,selSCC}.avgSpeed;
%     avgValsSCC(ann) = mean(valsSCC{an,selSCC}.absalluv);
    avgValsS(ann) = mean(valsS{an}.absalluv);
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
    for iii = 1:length(selLMH_all{ss}(an,:))
        ii = selLMH_all{ss}(an,iii);
        if animals(an) == 4
            theAng = angle(mean(valsE{an}.allang{ii}))+pi; theMag = abs(mean(valsE{an}.allang{ii}));
            avgValsE(an,iii) = conj(theMag * exp(1i * theAng));
        else
            avgValsE(an,iii) = conj(mean(valsE{an}.allang{ii}));
        end
    end
    if animals(an) == 4
        theAng = angle(mean(valsSCC{an,selSCC}.allang))+pi; theMag = abs(mean(valsSCC{an,selSCC}.allang));
        avgValsSCC(an) = conj(theMag * exp(1i * theAng));
        theAng = angle(mean(valsS{an}.allang))+pi; theMag = abs(mean(valsS{an}.allang));
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



function dur = plotAllCurves(figN,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
for an = 1:length(animals)
    lE(an) = length(valsE{an}.allThisdf);
end
gCols = (max(lE)+2);
gRows = length(animals);
totalPlots = gRows * gCols;
grcs = reshape(1:totalPlots,gCols,gRows)';
yl = 2;
amplitude_threshold = 0*meta.amplitude_threshold;
figure(figN);clf
eaaa = 0;
% animals = 3;
for ann = 1:length(animals)
    an = animals(ann);
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
            subplot(gRows,gCols,grcs(ann,2));hold on;
            avgPlot = valsSCC{an,selSCC}.avg.uv_chopped';
            if valsSCC{an,selSCC}.Fs == 220
                N = length(avgPlot);
                N1 = ceil(N*150/220);
                t = 0:N-1; t = t/220;
                t1 = 0:N1-1; t1 = t1/150;
                avgPlot = interp1(t,avgPlot,t1);
            end
            avgPlotA = avgPlot;
            avgPlot = mean(thisPlot);
            allAvgPlotsSCC(an,:) = avgPlotA/meta.mES(ann);
%             plot(thisPlot'/meta.mES(ann));%(:,27:(27+25))');
%             plot(avgPlot/meta.mES(ann),'k','linewidth',1.25);
            plot(avgPlotA/meta.mES(ann),'k','linewidth',1.25);
            box off;
            ylim([-0.1 yl]);
%             xlim([0 50]);
            title('Spont-CC','FontSize',10);
            plot([16 16],ylim,'-b');
        end
    end

    if ~isempty(valsS{an}.allThisdf)
        inds = valsS{an}.alldf > amplitude_threshold;
        thisPlot = valsS{an}.allThisuv(inds,:);
        if ~isempty(thisPlot)
            if valsS{an}.Fs == 220
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
            subplot(length(animals),gCols,grcs(ann,1));hold on;
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
            thisPlot = valsE{an}.allThisuv{aaa}(inds,:);
            thisPlotA = abs(valsE{an}.AvgallThisuv{1}{aaa})';
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
                    thisPlotA = interp1(t,thisPlotA,t1);
                end
                
                subplot(length(animals),gCols,grcs(ann,aaa+2));hold on;
                avgPlot = mean(thisPlot);
%                 avgPlot = abs(valsE{an}.AvgallThisuv{1}{aaa});
                allAvgPlotsE{an,aaa} = avgPlot/meta.mES(ann);
%                 plot(thisPlot'/meta.mES(ann));
%                 plot(avgPlot/meta.mES(ann),'k','linewidth',1.25);
                plot(thisPlotA/meta.mES(ann),'k','linewidth',1.25);
                box off;
                titleT = sprintf('A = %.3f',mdataE{ss}.amps(an,aaa));
                if aaa == 1
                    titleT = sprintf('A = %.3f - Lo',mdataE{ss}.amps(an,aaa));
                end
                if aaa == 2
                    titleT = sprintf('A = %.3f - Med',mdataE{ss}.amps(an,aaa));
                end
                if aaa == 3
                    titleT = sprintf('A = %.3f - Hi',mdataE{ss}.amps(an,aaa));
                end
                title(titleT,'FontSize',10);
                ylim([-0.1 yl]);
%                 xlim([0 50]);
                plot([16 16],ylim,'-b');
            end
        end
    end
end
save2pdf(makeName('allCurves.pdf',getpdffolder),gcf,600);
for an = 1:length(animals)
    thisCurve = allAvgPlotsSCC(an,:);
    durSCC(an,1) = findDuration(thisCurve,16,6.67);
%     thisCurve = allAvgPlotsS(an,:);
%     durS(an,1) = findDuration(thisCurve,16,6.67);
    for aaa = 1:size(allAvgPlotsE,2)
        if aaa == 1
            thisCurve = allAvgPlotsE{an,aaa};
            durE(an,1) = findDuration(thisCurve,16,6.67);
        end
        if aaa == 2
            thisCurve = allAvgPlotsE{an,aaa};
            durE(an,2) = findDuration(thisCurve,16,6.67);
        end
        if aaa == 3
            thisCurve = allAvgPlotsE{an,aaa};
            durE(an,3) = findDuration(thisCurve,16,6.67);
        end
    end
end
% figure(100);clf;
x_axis = (0:(size(allAvgPlotsSCC,2)-1)) * 1000/150 - (15000/150);
% x_axis = 1:size(allAvgPlotsSCC,2);
ff = makeFigureWindow__one_axes_only(100,[4 1 2.5 2],[0.15 0.2 0.84 0.75]);
thisCols = {cols.sponColorsN{1},cols.eSponColors{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3},cols.stimColors{4}};
avgAvgPlotSCC = mean(allAvgPlotsSCC);
avgAvgPlotSCC = avgAvgPlotSCC-min(avgAvgPlotSCC);%%%%%%
plot(x_axis,avgAvgPlotSCC,'color',thisCols{1});
semAvgPlotSCC = std(allAvgPlotsSCC)/sqrt(length(allAvgPlotsSCC));
shadedErrorBar(x_axis,avgAvgPlotSCC,semAvgPlotSCC,{'color',thisCols{1}},0.75);
% avgAvgPlotS = mean(allAvgPlotsS);
% plot(x_axis,avgAvgPlotS,'color',thisCols{1});
for aaa = 1:3
    allAvgPlotsE_LMH = [];
    for ann = 1:length(animals)
        an = animals(ann);
        allAvgPlotsE_LMH(an,:) = allAvgPlotsE{an,aaa};
    end
    avgAvgPlotE{aaa} = mean(allAvgPlotsE_LMH);
    avgAvgPlotE{aaa} = avgAvgPlotE{aaa} - min(avgAvgPlotE{aaa});%%%%%%%%%
    plot(x_axis,avgAvgPlotE{aaa},'color',thisCols{aaa+2});
    semAvgPlotE{aaa} = std(allAvgPlotsE_LMH)/sqrt(length(allAvgPlotsE_LMH));
    shadedErrorBar(x_axis,avgAvgPlotE{aaa},semAvgPlotE{aaa},{'color',thisCols{aaa+2}},0.75);
end
xlim([-40 max(x_axis)])
% plot([0 0],ylim,'--m');
ylim([-0.05 2]);
set(gca,'TickDir','out','FontSize',6);
set(gca,'linewidth',1.25,'Fontweight','bold','FontSize',6);
xlabel('Time (ms)');
ylabel('Avg Speed (mm/sec)');
thisCols = {cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3}};
x1 = 300; x2 = x1+20; y1 = (0.4:0.04:7); y1 = y1(1:4); y2 = y1;
legendFontSize = 6-1;
% legs = {'Spont-CC','Lo Evok','Med Evok','Hi Evok'};
legs = meta.stimNames;
for ii = 1:length(legs)
    plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols{ii});
    text(x2+10,y1(ii),legs{ii},'Color',thisCols{ii},'FontSize',legendFontSize);
end
n = 0;
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);
dur.durSCC = durSCC;
% dur.durS = durS;
dur.durE = durE;
% figure(10000);clf;
% for an = 1:length(animals)
%     subplot(1,4,1)
%     avgPlot = valsSCC{an,selSCC}.avg.uv_chopped';
%     plot(avgPlot,'k','linewidth',1.25);
%     box off;
%     ylim([-0.1 yl]);
% %             xlim([0 50]);
%     title('Spont-CC','FontSize',10);
%     plot([16 16],ylim,'-b');
%     for aaa = 1:length(valsE{an}.allThisdf)
%         if ~isempty(valsE{an}.allThisdf{aaa})
%             inds = valsE{an}.alldf{aaa} > amplitude_threshold;
%             thisPlot = valsE{an}.allThisdf{aaa}(inds,:);
%             thisPlot = abs(valsE{an}.AvgallThisuv{1}{aaa})';
%             if ~isempty(thisPlot)
%                 if valsE{an}.Fs == 220
%                     [K,N] = size(thisPlot);
%                     N1 = ceil(N*150/220);
%                     t = 0:N-1; t = t/220;
%                     t1 = 0:N1-1; t1 = t1/150;
%                     temp = [];
%                     for k=1:K
%                         temp(k,:) = interp1(t,thisPlot(k,:),t1);
%                     end
%                     thisPlot = temp;
%                 end
%                 subplot(length(animals),gCols,grcs(an,aaa+2));hold on;
% %                 avgPlot = mean(thisPlot);
%                 avgPlot = thisPlot;
% %                 avgPlot = abs(valsE{an}.AvgallThisuv{1}{aaa});
%                 allAvgPlotsE{an,aaa} = avgPlot;
% %                 plot(thisPlot');
%                 plot(avgPlot,'k','linewidth',1.25);
%                 box off;
%                 titleT = sprintf('A = %.3f',mdataE{ss}.amps(an,aaa));
%                 if aaa == 1
%                     titleT = sprintf('A = %.3f - Lo',mdataE{ss}.amps(an,aaa));
%                 end
%                 if aaa == 2
%                     titleT = sprintf('A = %.3f - Med',mdataE{ss}.amps(an,aaa));
%                 end
%                 if aaa == 3
%                     titleT = sprintf('A = %.3f - Hi',mdataE{ss}.amps(an,aaa));
%                 end
%                 title(titleT,'FontSize',10);
%                 ylim([-0.1 yl]);
% %                 xlim([0 50]);
%                 plot([16 16],ylim,'-b');
%             end
%         end
%     end
% end