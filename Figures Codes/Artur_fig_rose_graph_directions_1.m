function Artur_fig_rose_graph_directions_1(ss,roi)
close all
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
mdataE = evalin('base','mdataE');
mdataS = evalin('base','mdataSDFth');
meta = evalin('base','meta');
% data related settings
animals = [1 2 3 4 5];
selLMH_all = mdataE{ss}.selLMH.selLMH_all;
dfth = 2;
match = {'FL' 'HL'};
CCTh = 1;
cols = meta.colors;
selSCC = 3;
amplitude_threshold = meta.amplitude_threshold;
for an = 1:length(animals)
    valsE{an} = getDistributions_simpleAmps_artur_2(mdataE{ss},animals(an),roi,ss,amplitude_threshold);
    valsS{an} = getDistributionsSpon_simple_artur(mdataS{1},animals(an),roi,ss,dfth,amplitude_threshold);
    for LMH = 3%length(ssLMH)
        valsSCC{an,LMH} = getDistributionsSpon_motifs_simple_artur (mdataE,mdataS{1},animals(an),roi,match{ss},'Max',CCTh,LMH,amplitude_threshold);
    end
end
ff = makeFigureRowsCols(111,[3 3.2 5.9 1.75],'RowsCols',[1 4],'spaceRowsCols',[0.1 0.03],'rightUpShifts',[0.06 0.01],'widthHeightAdjustment',[-40 -80]);
if ss == 1
    sFactors = [0.85 0.85 0.85 0.85];
    sFactor = sFactors(roi);
end
if ss == 2
    sFactors = [0.85 0.8 0.85 0.85];
    sFactor = sFactors(roi);
end
% Spon motifs
rr = 1; cc = 1;
axes(ff.h_axes(rr,cc));
set(ff.h_axes(rr,cc),'Visible','off');hold on;
directionDistributions_spon (ff,rr,cc,ss,animals,valsSCC,selLMH_all,cols,meta,selSCC,sFactor)
% Evoked Lo
rr = 1; cc = 2;
axes(ff.h_axes(rr,cc));
set(ff.h_axes(rr,cc),'Visible','off');hold on;
directionDistributions_evoked (ff,rr,cc,ss,animals,valsE,selLMH_all,cols,meta,1,sFactor);
% Evoked Med
rr = 1; cc = 3;
axes(ff.h_axes(rr,cc));
set(ff.h_axes(rr,cc),'Visible','off');hold on;
directionDistributions_evoked (ff,rr,cc,ss,animals,valsE,selLMH_all,cols,meta,2,sFactor);
% Evoked Hi
rr = 1; cc = 4;
axes(ff.h_axes(rr,cc));
set(ff.h_axes(rr,cc),'Visible','off');hold on;
directionDistributions_evoked (ff,rr,cc,ss,animals,valsE,selLMH_all,cols,meta,3,sFactor);
% Coordinates
axes(ff.h_axes(1,1));
midX = 0.8;
midY = 1;
hL = 0.3; % half length of the cross
sp = 0.03; % space between letters and the ends of cross
orientation.top = 'A'; orientation.bottom = 'P'; orientation.left = 'M'; orientation.right = 'L';
FontSize = 9;
col = 'k';
CranialWinOrient(midX, midY, hL, sp, orientation, FontSize, col);
% appending pdf
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder));

ff = makeFigureWindow__one_axes_only(2,[9 3.2 1 1.75],[0.38 0.2 0.6 0.69]);
plotBarGraphDirection (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder));


function plotBarGraphDirection (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
% figure related settings
axesFontSize = 6; ylabelFontSize = 6;
for an = length(animals):-1:1
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
post_hoc = 'hsd';
annovaVar = [avgValsSCC' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3)];
annovaVar = abs(annovaVar);
disp('Diretion')
prT = ranovaFLHL(annovaVar);
thisCols = {cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3},cols.stimColors{4}};

% [~,~,stats] = anova1(annovaVar,{'SpontCC','Evok_1','Evok_2','Evok_3'},'off');
% [c,~,~,~] = multcompare(stats,'CType',post_hoc,'Display','off');
% prT = c(:,6);

hrT = prT<0.05;
annovaVar = [avgValsSCC' avgValsE(:,1) avgValsE(:,2) avgValsE(:,3)];
annovaVar = abs(annovaVar);
axes(ff.ha);
[mVals, semVals] = findMeanAndStandardError(annovaVar,1);
combs = nchoosek(1:size(annovaVar,2),2);
plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols,'ySpacingFactor',10);
xlim([0.4 0.6+size(annovaVar,2)]);
ylim([0 1.49]);
ylabel('Avg Vector Magnitude (A.U.)','FontSize',ylabelFontSize);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
set(ff.ha,'XTick',1:size(annovaVar,2));
set(ff.ha,'XTickLabel',meta.stimNames,'FontSize',axesFontSize);
xtickangle(ff.ha,45);
set(ff.ha,'Fontweight','bold','FontSize',axesFontSize);
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);


function directionDistributions_evoked (ff,rr,cc,ss,animals,valsE,selLMH_all,cols,meta,sL,sFactor)
nbin = 15*2;
r_max = 1;
ti = 0:0.0001:1;
r_max_lim = 1.25;
xcirc = r_max_lim*cos(2*pi*ti);
ycirc = r_max_lim*sin(2*pi*ti);
rE_all = zeros(length(animals),nbin*4);
for an = 1:length(animals)
    selLMH = selLMH_all{ss}(an,sL);
    angE = [];
    angE = angle(conj(valsE{an}.allang{selLMH}));
    if an == 4
        angE = angle(conj(-valsE{an}.allang{selLMH}));
    end
    [t,r] = rose(angE,nbin);
    r = r_max*r/max(r);
    rE_all(an,:) = r;
end
rE_all = sFactor*(r_max*rE_all/max(mean(rE_all)));
mE = mean(rE_all);
sdE = std(rE_all)/sqrt(length(animals));
% [ka ka0] = circ_kurtosis(t,mE)
% p_val = circ_symtest(mE)
Runique = mE(2:4:end); Runique(end+1) = Runique(1);
sdE = sdE(2:4:end); sdE(end+1) = sdE(1);
Rlow = Runique;
Rhigh = Runique+sdE;
Tunique = t(2:4:end); Tunique(end+1) = Tunique(1);
col = cols.stimLevelColors{sL};
%                     patchColor = col/1.2;%col+(1-col)*(1-patchSaturation);
for k = 1:numel(Runique)-1
    patch([cos(Tunique(k))*Rlow(k) cos(Tunique(k))*Rhigh(k) cos(Tunique(k+1))*Rhigh(k) cos(Tunique(k+1))*Rlow(k)], ...
        [sin(Tunique(k))*Rlow(k) sin(Tunique(k))*Rhigh(k) sin(Tunique(k+1))*Rhigh(k) sin(Tunique(k+1))*Rlow(k)],...
        1, 'facecolor', col, 'edgecolor', 'none','Facealpha',0.25);
end
p1h = polar(t,mean(rE_all));
set(p1h,'color',col)
hold on;
plot(xcirc,ycirc,'k')
set(ff.h_axes(rr,cc),'xlim',[-1 1]*r_max_lim)
set(ff.h_axes(rr,cc),'ylim',[-1 1]*r_max_lim)
axis equal
textString = meta.stimNames{sL+1};
htxt = text(0.5,0.95,textString,'unit','normalized','fontsize',meta.legendFontSize+3,'Rotation',0,'FontWeight','Bold',...
                    'HorizontalAlignment','center','VerticalAlignment','middle');

function directionDistributions_spon (ff,rr,cc,ss,animals,valsSCC,selLMH_all,cols,meta,selSCC,sFactor)   
nbin = 15*2;
r_max = 1;
ti = 0:0.0001:1;
r_max_lim = 1.25;
xcirc = r_max_lim*cos(2*pi*ti);
ycirc = r_max_lim*sin(2*pi*ti);
rSCC_all = zeros(length(animals),nbin*4);
aa=[];
for an = 1:length(animals)%[1 2 3 5]%
    angSCC = angle(conj(valsSCC{an,selSCC}.allang));
    pval(an) = circ_rtest(angSCC);
    if an == 4
        angSCC = angle(conj(-valsSCC{an,selSCC}.allang));
    end
    [t,r] = rose(angSCC,nbin);
    r = r_max*r/max(r);
    rSCC_all(an,:) = r;
    aa = [aa,angSCC];
end
disp(pval)
circ_rtest(aa)
rSCC_all = sFactor*(r_max*rSCC_all/max(mean(rSCC_all)));
mSCC = mean(rSCC_all);
sdSCC = std(rSCC_all)/sqrt(length(animals));

% [ka ka0] = circ_kurtosis(t,mSCC)
p_val = circ_symtest(mSCC);

Runique = mSCC(2:4:end); Runique(end+1) = Runique(1);
sdSCC = sdSCC(2:4:end); sdSCC(end+1) = sdSCC(1);
Rlow = Runique;
Rhigh = Runique+sdSCC;
Tunique = t(2:4:end); Tunique(end+1) = Tunique(1);
col = cols.sponColorsN{1};
patchColor = col;%col+(1-col)*(1-patchSaturation);
for k = 1:numel(Runique)-1
    patch([cos(Tunique(k))*Rlow(k) cos(Tunique(k))*Rhigh(k) cos(Tunique(k+1))*Rhigh(k) cos(Tunique(k+1))*Rlow(k)], ...
        [sin(Tunique(k))*Rlow(k) sin(Tunique(k))*Rhigh(k) sin(Tunique(k+1))*Rhigh(k) sin(Tunique(k+1))*Rlow(k)],...
        1, 'facecolor', patchColor, 'edgecolor', 'none','Facealpha',0.25);
end
p1h = polar(t,mean(rSCC_all));
set(p1h,'color',col)
n=0;
plot(xcirc,ycirc,'k')

set(ff.h_axes(rr,cc),'xlim',[-1 1]*r_max_lim)
set(ff.h_axes(rr,cc),'ylim',[-1 1]*r_max_lim)
axis equal
textString = meta.stimNames{1};
yl = ylim;
htxt = text(0.5,0.95,textString,'unit','normalized','fontsize',meta.legendFontSize+3,'Rotation',0,'FontWeight','Bold',...
                    'HorizontalAlignment','center','VerticalAlignment','middle');

textString = sprintf('%d%%',round(100*1/sFactor));
htxt = text(0.1,0.8,textString,'unit','normalized','fontsize',meta.legendFontSize,'Rotation',45,'FontWeight','Bold',...
                    'HorizontalAlignment','center','VerticalAlignment','middle');
              