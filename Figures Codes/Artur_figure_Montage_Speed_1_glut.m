function Artur_figure_Montage_Speed_1_glut (ss)

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


% run the loadData.m file before this
mdataE = evalin('base','mdataE_glut_HPF');
mdataS = evalin('base','mdataS_glut_HPF');
meta = evalin('base','meta');
cols = meta.colors;
if ss == 1
    an = 2; roi = 2; match = 'AC';
end
stepVF = 7;

Cmap = meta.Cmap;
amplitude_threshold = meta.amplitude_threshold;
CCTh = 1;
selLMH_all = mdataE{ss}.selLMH.selLMH_all;
ln = selLMH_all{ss}(an,:);
% bregma = getBregma(mdataE{ss}.animalNumber{an});


temp = getDistributionsSpon_motifs_simple_artur_glut_2_HPF (mdataE,mdataS{1},an,roi(ss),'AC','Max',CCTh,3,amplitude_threshold);
valsS = getdfbyfspon27_2_glut(mdataE{ss}.animalNumber(an),temp);

valsEE = getDistributions_simpleAmps_artur_2_glut_HPF(mdataE{ss},an,roi(ss),ss,amplitude_threshold);
for ii = 1:3
    [lists names]= getEvokedListsDiffStimAmps_glut(mdataE{ss}.animalNumber(an));
    valsE{ii} = getdfbyfEvoked27_1_glut(mdataE{ss}.animalNumber(an),names{4}(ln(ii)),valsEE.avg_onsets(selLMH_all{ss}(an,1)));
end
legs = meta.stimNames;
maskI = mdataE{ss}.masks{an}{1};
ff = makeFigureRowsCols(110,[1 1 7 3],'RowsCols',[4 9],'spaceRowsCols',[0.011 0.0009],'rightUpShifts',[0.03 0.05],'widthHeightAdjustment',[-10 -25]);
thisCol = {cols.sponColorsN{1},cols.stimLevelColors{1},cols.stimLevelColors{2},cols.stimLevelColors{3}};
rr = 1;
Imin = 0; 
if ss == 1
    Imax = 0.9; scale_factor = [20 100 20 20];
end
for cc = 1:9
    axes(ff.h_axes(rr,cc));hold on;
    thisFrame = applyMask(valsS.df(:,:,cc),maskI,'shiftToZero',1);
    imagesc(thisFrame,[Imin Imax]);
    colormap(Cmap);
    axis equal; axis off;
    if cc == 1
        text(-25,65,legs{rr},'rotation',90,'FontSize',meta.axesFontSize,...
            'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','cap','color',thisCol{rr});
        text(5,100,num2str(scale_factor(rr)),'rotation',90,'FontSize',meta.axesFontSize,...
                'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','cap','color','w');
    end
    [u,v,X,Y] = trim_uv(valsS.uv(:,:,cc),stepVF,maskI);
    quiver(X,Y,u*scale_factor(rr),v*scale_factor(rr),'AutoScale','off','color','w');
    set(ff.h_axes(rr,cc),'Ydir','reverse');
    [xMask,yMask] = findMaskBorder(maskI);
    plot(xMask, yMask, '-w','linewidth',0.25);
end
for rr = 2:4
    for cc = 1:9
        axes(ff.h_axes(rr,cc));hold on;
        thisFrame = applyMask(valsE{rr-1}.df(:,:,cc),maskI,'shiftToZero',1);
        imagesc(thisFrame,[Imin Imax]);
        colormap(Cmap);
        axis equal; axis off;
        if cc == 1
            text(-25,65,legs{rr},'rotation',90,'FontSize',meta.axesFontSize,...
                'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','cap','color',thisCol{rr});
            text(5,100,num2str(scale_factor(rr)),'rotation',90,'FontSize',meta.axesFontSize,...
                'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','cap','color','w');
        end
        [u,v,X,Y] = trim_uv(valsE{rr-1}.uv(:,:,cc),stepVF,maskI);
        quiver(X,Y,u*scale_factor(rr),v*scale_factor(rr),'AutoScale','off','color','w');
        set(ff.h_axes(rr,cc),'Ydir','reverse');
        [xMask,yMask] = findMaskBorder(maskI);
        plot(xMask, yMask, '-w','linewidth',0.25);
    end
end

for rr = 1:4
    ha = axes('Position',ff.axesPos{rr,9}+[0.065 0 0 0]);
    axis off;
    h_cbar = colorbar('eastoutside');
    set(h_cbar,'ytick',[])
    yl = ylim; xl = xlim; dxl = xl(2)-xl(1);
    text(xl(2)+3*dxl,yl(1)+0.05*dxl,sprintf('%.1f',0.1*floor(10*Imin)),'color','k','FontSize',meta.axesFontSize-1,'HorizontalAlignment','center','rotation',0);
    text(xl(2)+3*dxl,yl(2)-0.05*dxl,num2str(round(Imax,2)),'color','k','FontSize',meta.axesFontSize-1,'HorizontalAlignment','center','rotation',0);
    text(xl(2)+3*dxl,yl(2)-0.5*dxl,'\DeltaF/F_0 %','color','k','FontSize',meta.axesFontSize-1,'HorizontalAlignment','center','rotation',270)
end
sFF = 1;
axes(ff.h_axes(1,1));hold on;
plot(80+[0 30], 10*[1 1], 'w','linewidth',0.5)
% text(7,110,'2 mm','color','w','FontSize',axesFontSize-sFF)
midX = 73;
midY = 67;
hL = 17; % half length of the cross
sp = 3; % space between letters and the ends of cross
plot(midX+[-hL hL], midY*[1 1], 'w','linewidth',0.5)
plot(midX*[1 1], midY+[-hL hL], 'w','linewidth',0.5)
text(midX,midY-hL-sp,'A','color','w','FontSize',meta.axesFontSize-sFF,'HorizontalAlignment','center','VerticalAlignment','baseline');
text(midX,midY+hL+sp,'P','color','w','FontSize',meta.axesFontSize-sFF,'HorizontalAlignment','center','VerticalAlignment','cap');
text(midX+hL+sp,midY,'L','color','w','FontSize',meta.axesFontSize-sFF,'HorizontalAlignment','left','VerticalAlignment','middle');
text(midX-hL-sp,midY,'M','color','w','FontSize',meta.axesFontSize-sFF,'HorizontalAlignment','right','VerticalAlignment','middle');


for cc = 1:9
    axes(ff.h_axes(4,cc));
    imgT = round(10*((cc-2)*1000/150))/10;
    text(10,140,[num2str(imgT),' ms'],'FontSize',meta.axesFontSize,'color','k');
end

axes(ff.h_axes(4,3));hold on;
thisROI = mdataE{ss}.masks{an}{roi};
bo = bwboundaries(thisROI);
xs = bo{1}(:,1);
ys = bo{1}(:,2);
plot(ys,xs,'m','linewidth',0.5);

n=0;

save2pdf(pdfFileName,ff.hf,600);

function [u,v,X,Y] = trim_uv(uv,step,mask)
u = real(uv);        v = imag(uv);
lenx = size(u,2);        leny = size(u,1);
[X Y] = meshgrid(1:lenx,1:leny);
idxs = 1:step:lenx;        idys = 1:step:leny;
X = X(idxs,idys);        Y = Y(idxs,idys);
u = u(idxs,idys);        v = v(idxs,idys);
bo = bwboundaries(mask);
xs = bo{1}(:,1);        ys = bo{1}(:,2);
in = inpolygon(X,Y,ys,xs);
Y = Y(in); X = X(in);
u = u(in);        v = v(in);