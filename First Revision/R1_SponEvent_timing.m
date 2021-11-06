function R1_SponEvent_timing
clc;close all;
match = {'FL' 'HL','VC'};
pdfFileName = sprintf('%s.pdf',mfilename);
pdfFileName = makeName(pdfFileName,getpdffolder)
if exist(pdfFileName,'file')
    delete(pdfFileName);
end
if exist(makeName('temp.pdf',getpdffolder),'file')
    delete(makeName('temp.pdf',getpdffolder));
end

% distribution calculation parameters
nbins = 20;
ninterp = 60;
minX = 0;
maxX = 5;
Lx = maxX - minX;
maxY = 1.25;
xbins = linspace(minX,maxX,nbins);
cInterp = linspace(minX,maxX,ninterp);

% FL and HL
% load data from base work space where all the data has to be pre loaded
mdataE = evalin('base','mdataE');
mdataS = evalin('base','mdataSDFth');
meta = evalin('base','meta');
% data related settings
animals = [1 2 3 4 5];
roi = [2 3];
dfth = 2;
CCTh = 1;
cols = getColors;
selSCC = 3;
amplitude_threshold = meta.amplitude_threshold;
% mfilename
reWrite = 0;
if reWrite
    for ss = 1:2
        fileName = sprintf('%s_%d.mat',mfilename,ss);
        for an = 1:length(animals)
            valsE{an} = getDistributions_simpleAmps_artur_2(mdataE{ss},animals(an),roi(ss),ss,amplitude_threshold);
            valsS{an} = getDistributionsSpon_simple_artur(mdataS{1},animals(an),roi(ss),ss,dfth,amplitude_threshold);
            for LMH = 3%length(ssLMH)
                valsSCC{an,LMH} = getDistributionsSpon_motifs_simple_artur (mdataE,mdataS{1},animals(an),roi(ss),match{ss},'Max',CCTh,LMH,amplitude_threshold);
            end
        end
        save(fileName,'valsE','valsS','valsSCC');
    end
end
for ss = 1:2
    fileName = sprintf('%s_%d.mat',mfilename,ss);
    load(fileName);
    ieiS = []; distS = []; durS = [];
    for an = 1:length(animals)
        % inter event interval calculation
        iei = diff(valsS{an}.allFrames(:,3)); % inter event interval
        fd = diff(valsS{an}.allFrames(:,2));  % file difference
        fd = ~fd; % to remove iei between files.
        iei = iei(fd);
        ieiS{an} = iei/150; % to sec
        distS(an,:) = calcDist(ieiS{an}, xbins, cInterp);
        % event rate calculation
        nS(an) = size(valsS{an}.allFrames,1);
        dur = 0;
        for n = 1:size(mdataS{1, 1}.d,2)
            if ~isempty(mdataS{1}.d{an,n})
                dur = dur + length(mdataS{1}.d{an,n}{1});
            end
        end
        durS(an) = dur/150; % to sec
    end
    if ss == 1
        ieiFL = ieiS;
        distFL = distS;
        erFL = nS./durS; % event rate for FL
    else
        ieiHL = ieiS;
        distHL = distS;
        erHL = nS./durS; % event rate for HL
    end
end
% VC
ss = 3;
mdataE = evalin('base','mdataE_VC');
mdataS = evalin('base','mdataS_VC');
meta = evalin('base','meta');
animalNumber = mdataE{ss}.animalNumber;
animals = [1 2 3 4];
amplitude_threshold = meta.amplitude_threshold;
roi = 1;
dfth = 2;
reWrite = 0;
fileName = sprintf('%s_%d.mat',mfilename,ss);
if reWrite
    for an = 1:length(animals)
        valsE{an} = getDistributions_simpleAmps_artur_2(mdataE{ss},animals(an),roi,ss,amplitude_threshold);
        valsS{an} = getDistributionsSpon_simple_artur(mdataS{1},animals(an),roi,ss,dfth,amplitude_threshold);
        for LMH = 3%length(ssLMH)
            valsSCC{an,LMH} = getDistributionsSpon_motifs_simple_artur (mdataE,mdataS{1},animals(an),roi,match{ss},'Max',CCTh,LMH,amplitude_threshold,0);
        end
    end
    save(fileName,'valsE','valsS','valsSCC');
end
load(fileName);
nS = []; durS = [];
for an = 1:length(animals)
    % inter event interval calculation
    iei = diff(valsS{an}.allFrames(:,3)); % inter event interval
    fd = diff(valsS{an}.allFrames(:,2));  % file difference
    fd = ~fd; % to remove iei between files.
    iei = iei(fd);
    ieiVC{an} = iei/150; % to sec
    distVC(an,:) = calcDist(ieiVC{an}, xbins, cInterp);
    % event rate calculation
    nS(an) = size(valsS{an}.allFrames,1);
    dur = 0;
    for n = 1:size(mdataS{1, 1}.d,2)
        if ~isempty(mdataS{1}.d{an,n})
            dur = dur + length(mdataS{1}.d{an,n}{1});
        end
    end
    durS(an) = dur/150; % to sec
end
erVC = nS./durS;  % event rate for VC
        
%% AC Glut
match = {'AC'};
ss = 1;
mdataE = evalin('base','mdataE_glut_BPF_050Hz_6Hz');
mdataS = evalin('base','mdataS_glut_BPF_050Hz_6Hz');
meta = evalin('base','meta');
animals = [1 2 3 4 5];
amplitude_threshold = meta.amplitude_threshold;
roi = 1;
dfth = 0;
reWrite = 0;
fileName = sprintf('%s_%d.mat',mfilename,ss+3);
if reWrite
    for an = 1:length(animals)
        valsE{an} = getDistributions_simpleAmps_artur_2_glut(mdataE{ss},animals(an),roi(ss),ss,amplitude_threshold);
        valsS{an} = getDistributionsSpon_simple_artur_glut(mdataS{1},animals(an),roi(ss),ss,dfth,amplitude_threshold);
        for LMH = 3%length(ssLMH)
            valsSCC{an,LMH} = getDistributionsSpon_motifs_simple_artur_glut_2 (mdataE,mdataS{1},animals(an),roi(ss),'AC','Max',CCTh,LMH,amplitude_threshold);
        end
    end
    save(fileName,'valsE','valsS','valsSCC');
end
load(fileName);
nS = []; durS = [];
for an = 1:length(animals)
    % inter event interval calculation
    iei = diff(valsS{an}.allFrames(:,3)); % inter event interval
    fd = diff(valsS{an}.allFrames(:,2));  % file difference
    fd = ~fd; % to remove iei between files.
    iei = iei(fd); iei(iei<=10) = [];
    ieiAC{an} = iei/100; % to sec
    distAC(an,:) = calcDist(ieiAC{an}, xbins, cInterp);
    % event rate calculation
    nS(an) = size(valsS{an}.allFrames,1);
    dur = 0;
    for n = 1:size(mdataS{1, 1}.d,2)
        if ~isempty(mdataS{1}.d{an,n})
            dur = dur + length(mdataS{1}.d{an,n}{1});
        end
    end
    durS(an) = dur/100; % to sec
end
erAC = nS./durS; % event rate for AC Glut

% stat on event rate and plotting bar graph
annovaVar = [erFL, erHL, erVC(1:3), erAC];
g1 = {'FL','FL','FL','FL','FL',...
      'HL','HL','HL','HL','HL',...
      'VC','VC','VC',...
      'AC','AC','AC','AC','AC'};
g2 = {'VSD1','VSD1','VSD1','VSD1','VSD1',...
      'VSD1','VSD1','VSD1','VSD1','VSD1',...
      'VSD2','VSD2','VSD2',...
      'Glut','Glut','Glut','Glut','Glut'};
% p = anovan(annovaVar,{g1,g2}); % in case you want to check for animal groups
thisCols = {'k','k','k','k'};
thisCols = {'b','r','g','k'};
[p,tbl,stats] = anova1(annovaVar,g1,'off');
[c,~,~,gnames] = multcompare(stats,'CType','hsd');
prT = c(:,6);
hrT = prT<0.05;
mVals = [mean(erFL), mean(erHL), mean(erVC(1:3)), mean(erAC)];
semVals = [std(erFL)/sqrt(5), std(erHL)/sqrt(5), std(erVC(1:3))/sqrt(3), std(erAC)/sqrt(5)];
combs = nchoosek(1:size(mVals,2),2);

close all
axesFontSize = 6;
ff = makeFigureWindow__one_axes_only(4,[1 1 1.5 1.75],[0.25 0.2 0.7 0.68]);
plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols,'ySpacingFactor',10);
xlim([0.4 0.6+size(mVals,2)]);
hy = ylabel('Spon Event Rate (event/s)');
pos = get(hy,'Position');pos(1) = pos(1) + 0.0;set(hy,'Position',pos);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
set(ff.ha,'XTick',1:size(mVals,2));
set(ff.ha,'XTickLabel',{'FL','HL', 'VC','AC'});
set(ff.ha,'Fontweight','bold');
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))

% plot inter-event-interval distributions for all modalities
ff = makeFigureWindow__one_axes_only(4,[1 1 1.5 1.75],[0.25 0.2 0.7 0.68]);
hold on;
linewidth = 1; faceAlpha = 0.5;
% FL
ss = 1;
meanS = mean(distFL);
semS = std(distFL)/sqrt(5); % 5 animals
shadedErrorBar(cInterp,meanS,semS,{'color',thisCols{ss},'linewidth',linewidth},faceAlpha);
% HL
ss = 2;
meanS = mean(distHL);
semS = std(distHL)/sqrt(5); % 5 animals
shadedErrorBar(cInterp,meanS,semS,{'color',thisCols{ss},'linewidth',linewidth},faceAlpha);
% VC
ss = 3;
meanS = mean(distVC);
semS = std(distVC)/sqrt(4); % 4 animals
shadedErrorBar(cInterp,meanS,semS,{'color',thisCols{ss},'linewidth',linewidth},faceAlpha);
% AC
ss = 4;
meanS = mean(distAC);
semS = std(distAC)/sqrt(5); % 5 animals
shadedErrorBar(cInterp,meanS,semS,{'color',thisCols{ss},'linewidth',linewidth},faceAlpha);

xlim([minX,maxX]); ylim([0,maxY]);
set(ff.ha,'XTick',[0:maxX],'YTick',[0:0.5:maxY])
xlabel('Inter-Event-Interval (s)');
hy = ylabel('Probability Density Function');
pos = get(hy,'Position');pos(1) = pos(1) + 0.0;set(hy,'Position',pos);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
set(ff.ha,'Fontweight','bold');

x1 = 3; x2 = x1+0.5; y1 = (0.5:0.17:1.25); y1 = y1(1:4); y2 = y1;
legendFontSize = 6;
%legs = meta.stimNames;
legs = {'FL','HL', 'VC','AC'};
for ii = 1:length(legs)
    plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols{ii},'linewidth',2);
    text(x2+0.1,y1(ii),legs{ii},'Color',thisCols{ii},'FontSize',legendFontSize+2);
end

save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))


if exist(makeName('temp.pdf',getpdffolder),'file')
    delete(makeName('temp.pdf',getpdffolder));
end



function [ci] = calcDist(bVar, xbins, cInterp)
[cE,centers] = hist(bVar,xbins);
ci = interp1(centers, cE, cInterp);
dx = cInterp(2)-cInterp(1);
ci = ci/sum(ci)/dx;
