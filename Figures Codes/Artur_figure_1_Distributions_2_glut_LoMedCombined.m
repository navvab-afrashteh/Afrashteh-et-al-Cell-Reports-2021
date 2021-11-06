function Artur_figure_1_Distributions_2_glut_LoMedCombined(ss)
ss = 1;
close all
clc
% all units are in percentage
ff.numberOfRows = 1; % specify how many rows do you want in the figure
ff.spaceBetweenRows = [0.05]; % specify based on number of rows what would be the spacing between rows
ff.rowHeights = [0.75]; % specify row heights ... should have the same number as number of rows

ff.numberOfCols = [1]; % for each row specify the number of columns
ff.colWidths = {[0.85]}; % for each row specify the column widths
ff.spaceBetweenCols = {[0.001]}; % specify space between columns

ff.leftOffset = 0.13; % how much the whole figure is offset to left from figure's left border
ff.bottomOffset = 0.21; % how much the whole figure is offset upwards from the bottom border of the figure

ff = getPanelsProps(ff); % will give a structure output which has all the lefts and bottoms of panels

figNum = 1; figPosition = [1 5 3 1.75];

ff.hf = makeFigureWindow(figNum,figPosition,1);

% ff = makeAxes(ff.hf,ff,1,1,[0 0 0 0]); % for generating axes

% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;

% load data from base work space where all the data has to be pre loaded
mdataE = evalin('base','mdataE_glut_BPF_050Hz_6Hz');
mdataS = evalin('base','mdataS_glut_BPF_050Hz_6Hz');
meta = evalin('base','meta');

% data related settings
animals = [1 2 3 4 5];
selLMH_all = mdataE{ss}.selLMH.selLMH_all;
selLMH = [1 2 3];
roi = [2];
dfth = 0;
% match = {'FL' 'HL'};
CCTh = 1;
cols = getColors;
selSCC = 3;
amplitude_threshold = meta.amplitude_threshold;
% mfilename

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

minY = 0;
maxY = 3;
nbins = 20;
ninterp = 45*2;
minX = 0;
maxX = 2;
xbins = linspace(minX,maxX,nbins);
cInterp = linspace(minX,maxX,ninterp);
cE_all = zeros(length(animals),ninterp);
cSCC_all = zeros(length(animals),ninterp);
cS_all = zeros(length(animals),ninterp);
for an = 1:length(animals)
    ann = animals(an);
    bVarEC = [];
    for ii = [1,3]%1:length(selLMH)
        if ii == 1
            bVarEA{ii} = [valsE{ann}.alldf{ii}, valsE{ann}.alldf{ii+1}];
        else
            bVarEA{ii} = valsE{ann}.alldf{ii};
        end
        bVarEC = [bVarEC bVarEA{ii}];
        avgValsE(an,ii) = mean(bVarEA{ii});
    end
    amplitude_threshold = 0;
    inds = valsSCC{ann,selSCC}.alldf > amplitude_threshold;
    if ~isempty(CorrMaxDF_frames)
        TrjFrames = CorrMaxDF_frames{an}';
        CorrFrames = valsSCC{an, 3}.allFrames(:,3);
        idx = ismember(CorrFrames,TrjFrames);
        inds = inds & idx';
    end
    bVarSCC = valsSCC{ann,selSCC}.alldf(inds);
    avgValsSCC(an) = mean(bVarSCC);
    bVarS = valsS{ann}.alldf;
    avgValsS(an) = mean(valsS{ann}.alldf);
    for ii = [1,3]%1:length(selLMH)
        bVarE = bVarEA{ii};
        [cEi, cSCCi, cSi] = calcCounts(bVarE, bVarSCC, bVarS, xbins, cInterp);
        cEAi{ii} = cEi;
    end
    for ii = [1,3]%1:length(selLMH)
        cEA(ii).cE_all(an,:) = cEAi{ii};
    end
    [cECi, cSCCi, cSi] = calcCounts(bVarEC, bVarSCC, bVarS, xbins, cInterp);
    cEC_all(an,:) = cECi;
    cSCC_all(an,:) = cSCCi;
    cS_all(an,:) = cSi;
%     [pan(an),han(an),stat] = ranksum(bVarEC,bVarSCC,'tail','right');
end
for ii = [1,3]%1:length(selLMH)
    cEmean{ii} = mean(cEA(ii).cE_all,1);
end
cECmean = mean(cEC_all,1);
cSCCmean = mean(cSCC_all,1);
cSmean = mean(cS_all,1);
for ii = [1,3]%1:length(selLMH)
    cEsem{ii} = std(cEA(ii).cE_all,[],1)/sqrt(length(animals));
end
cECsem = std(cEC_all,[],1)/sqrt(length(animals));
cSCCsem = std(cSCC_all,[],1)/sqrt(length(animals));
cSsem = std(cS_all,[],1)/sqrt(length(animals));
for ii = 1:5
    pLp3s(ii) = sum(cSCC_all(ii,find(cInterp<0.3)));
    pLap3e(ii) = sum(cEC_all(ii,find(cInterp>0.3)));
end
ff = makeAxes(ff.hf,ff,1,1,[0 0 0 0]);
hold on;
faceAlpha = 0.75;
for ii = [1,3]%1:length(selLMH)
    axes(ff.ha); shadedErrorBar(cInterp,cEmean{ii},cEsem{ii},{'color',cols.stimLevelColors{selLMH(ii)}/1.1});
end
axes(ff.ha); shadedErrorBar(cInterp,cSCCmean,cSCCsem,{'color',cols.sponColorsN{1}/1.1});
plot(cInterp,cSCCmean,'color',cols.sponColorsN{1},'linewidth',1.25)
for ii = [1,3]%1:length(selLMH)
    plot(cInterp,cEmean{ii},'color',cols.stimLevelColors{selLMH(ii)},'linewidth',1.25)
end
axes(ff.ha);
ylabel('Probability Density Function','FontSize',ylabelFontSize);
xlabel('Glu Signal Amplitude(\DeltaF/F_0 %)','FontSize',ylabelFontSize);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
xlim([minX maxX]);
ylim([minY maxY]);
axes(ff.ha);
thisCols = {cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{3}};
x1 = 1.5; x2 = x1+0.12; y1 = (1:0.5:3); y1 = y1(1:4); y2 = y1;
legs = meta.stimNames([1,2,4]);
for ii = 1:length(legs)
    plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols{ii},'linewidth',2);
    text(x2+0.02,y1(ii),legs{ii},'Color',thisCols{ii},'FontSize',legendFontSize+2);
end
set(gca,'linewidth',1.25,'FontSize',axesFontSize);
save2pdf(makeName(sprintf('fig_distributions_glut_LoMedCombined_%d.pdf',ss),getpdffolder),ff.hf,600);

function [cEi, cSCCi, cSi] = calcCounts(bVarE, bVarSCC, bVarS, xbins, cInterp)
[cE,centers] = hist(bVarE,xbins);
[cSCC,~] = hist(bVarSCC,xbins);
[cS,~] = hist(bVarS,xbins);

cEi = interp1(centers, cE, cInterp);
cSCCi = interp1(centers, cSCC, cInterp);
cSi = interp1(centers, cS, cInterp);

dx = cInterp(2)-cInterp(1);
cEi = 1*cEi/sum(cEi)/dx;
cSCCi = 1*cSCCi/sum(cSCCi)/dx;
cSi = 1*cSi/sum(cSi)/dx;

