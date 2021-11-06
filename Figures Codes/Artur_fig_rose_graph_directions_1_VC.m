function Artur_fig_rose_graph_directions_1_VC(ss,roi)

if ~exist('ss','var')
    ss = 3; roi = 5;
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
mdataE = evalin('base','mdataE_VC');
mdataS = evalin('base','mdataS_VC');
meta = evalin('base','meta');

% data related settings
animals = [1 2 3 4];
selLMH_all = mdataE{ss}.selLMH.selLMH_all;
dfth = 2;
match = {'FL' 'HL' 'VC'};
CCTh = 1;
cols = getColors;
selSCC = 3;
amplitude_threshold = meta.amplitude_threshold;
selLMH = 1;

for an = 1:length(animals)
    valsE{an} = getDistributions_simpleAmps_artur_2(mdataE{ss},animals(an),roi,ss,amplitude_threshold);
    valsS{an} = getDistributionsSpon_simple_artur(mdataS{1},animals(an),roi,ss,dfth,amplitude_threshold);
    for LMH = 3%length(ssLMH)
        valsSCC{an,LMH} = getDistributionsSpon_motifs_simple_artur (mdataE,mdataS{1},animals(an),roi,match{ss},'Max',CCTh,LMH,amplitude_threshold,0);
    end
end

ff = makeFigureRowsCols(111,[6 2 3 1.75],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.03],'rightUpShifts',[0.03 0.01],'widthHeightAdjustment',[-40 -80]);
sFactor = 0.95;
% Spon motifs
rr = 1; cc = 1;
axes(ff.h_axes(rr,cc));
set(ff.h_axes(rr,cc),'Visible','off');hold on;
directionDistributions_spon (ff,rr,cc,ss,animals,valsSCC,selLMH_all,cols,meta,selSCC,sFactor)
% Evoked
rr = 1; cc = 2;
axes(ff.h_axes(rr,cc));
set(ff.h_axes(rr,cc),'Visible','off');hold on;
directionDistributions_evoked (ff,rr,cc,ss,animals,valsE,selLMH_all,cols,meta,1,sFactor);
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
% save as pdf
save2pdf(pdfFileName,ff.hf,600);

function directionDistributions_evoked (ff,rr,cc,ss,animals,valsE,selLMH_all,cols,meta,sL,sFactor)
nbin = 15*2;
r_max = 1;
ti = 0:0.0001:1;
r_max_lim = 1.25;
xcirc = r_max_lim*cos(2*pi*ti);
ycirc = r_max_lim*sin(2*pi*ti);
rE_all = zeros(length(animals),nbin*4);
selLMH = sL;%selLMH_all{ss}(an,sL);
for an = 1:length(animals)
    angE = angle(conj(valsE{an}.allang{selLMH}));
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
col = cols.stimLevelColors{3};
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
textString = meta.stimNamesVC{2};
htxt = text(0.5,1,textString,'unit','normalized','fontsize',meta.legendFontSize+3,'Rotation',0,'FontWeight','Bold',...
                    'HorizontalAlignment','center','VerticalAlignment','middle');


function directionDistributions_spon (ff,rr,cc,ss,animals,valsSCC,selLMH_all,cols,meta,selSCC,sFactor)   
nbin = 15*2;
r_max = 1;
ti = 0:0.0001:1;
r_max_lim = 1.25;
xcirc = r_max_lim*cos(2*pi*ti);
ycirc = r_max_lim*sin(2*pi*ti);
rSCC_all = zeros(length(animals),nbin*4);
for an = 1:length(animals)%[1 2 3 5]%
    angSCC = angle(conj(valsSCC{an,selSCC}.allang));
    [t,r] = rose(angSCC,nbin);
    r = r_max*r/max(r);
    rSCC_all(an,:) = r;
end
rSCC_all = sFactor*(r_max*rSCC_all/max(mean(rSCC_all)));
mSCC = mean(rSCC_all);
sdSCC = std(rSCC_all)/sqrt(length(animals));

%                     [ka ka0] = circ_kurtosis(t,mSCC)
p_val = circ_symtest(mSCC)

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
htxt = text(0.5,1,textString,'unit','normalized','fontsize',meta.legendFontSize+3,'Rotation',0,'FontWeight','Bold',...
                    'HorizontalAlignment','center','VerticalAlignment','middle');

textString = sprintf('%d%%',round(100*1/sFactor));
htxt = text(0.1,0.85,textString,'unit','normalized','fontsize',meta.legendFontSize,'Rotation',45,'FontWeight','Bold',...
                    'HorizontalAlignment','center','VerticalAlignment','middle');
