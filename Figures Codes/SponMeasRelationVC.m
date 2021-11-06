function SponMeasRelationVC(roi)
close all
clc

if ~exist('roi')
    roi = 1;
end
match = {'FL' 'HL' 'VC'};
ss = 3;
mfn = mfilename;
pdfFileName = sprintf('%s_%s_roi_%d.pdf',mfn(1:end-2),match{ss},roi);
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

cols = meta.colors;
% data related settings
animals = [1 2 3 4];
selLMH_all = mdataE{ss}.selLMH.selLMH_all;
dfth = 2;
CCTh = 1;
cols = getColors;
selSCC = 3;
amplitude_threshold = meta.amplitude_threshold;
% mfilename

for an = 1:length(animals)
    valsE{an} = getDistributions_simpleAmps_artur_2(mdataE{ss},animals(an),roi,ss,amplitude_threshold);
    valsS{an} = getDistributionsSpon_simple_artur(mdataS{1},animals(an),roi,ss,dfth,amplitude_threshold);
    for LMH = 3%length(ssLMH)
        valsSCC{an,LMH} = getDistributionsSpon_motifs_simple_artur (mdataE,mdataS{1},animals(an),roi,match{ss},'Max',CCTh,LMH,amplitude_threshold,0);
    end
end

% Amplitude and Speed relation
ff = makeFigureWindow__one_axes_only(2,[1 4 1.75 1.5],[0.21 0.25 0.7 0.69]);
plotAmpSpeed (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder));

% Amplitude and directionality relation
ff = makeFigureWindow__one_axes_only(2,[1 4 1.75 1.5],[0.21 0.25 0.7 0.69]);
plotAmpDirection (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder));

% Amplitude and Fractal Dimension
ff = makeFigureWindow__one_axes_only(2,[1 4 1.75 1.5],[0.21 0.25 0.7 0.69]);
plotAmpFD (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder));


function plotAmpSpeed (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; 
xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;

ampS = []; spdS = []; ampE = []; spdE = [];
for an = length(animals):-1:1
    ampS = [ampS,valsSCC{an,selSCC}.alldf];
    spdS = [spdS,valsSCC{an,selSCC}.absalluv];
    ampE = [ampE,valsE{an}.alldf{1}];
    spdE = [spdE,valsE{an}.absalluv{1}];
end
thL = mean(spdS)-3*std(spdS); thH = mean(spdS)+3*std(spdS); 
idxS = (spdS >= thL) & (spdS <= thH);
ampS = ampS(idxS);
spdS = spdS(idxS);
thL = mean(spdE)-3*std(spdE); thH = mean(spdE)+3*std(spdE); 
idxS = (spdE >= thL) & (spdE <= thH);
ampE = ampE(idxS);
spdE = spdE(idxS);
% linear model
% spon
mdlSpeedS = fitlm(ampS,spdS);
RsquaredE = mdlSpeedS.Rsquared.Adjusted;
pS = mdlSpeedS.anova.pValue(1);
bS = mdlSpeedS.Coefficients.Estimate(1);
mS = mdlSpeedS.Coefficients.Estimate(2);
xM = max([ampS,ampE]);
xM = 0.01*ceil(100*xM);
xm = 0.15;
xS = [xm,xM];
yS = mS*xS+bS;
% evoked
mdlSpeedE = fitlm(ampE,spdE);
mdlSpeedE = fitlm(ampE,spdE,'linear','RobustOpts','on');
RsquaredE = mdlSpeedE.Rsquared.Adjusted;
pE = mdlSpeedE.anova.pValue(1);
bE = mdlSpeedE.Coefficients.Estimate(1);
mE = mdlSpeedE.Coefficients.Estimate(2);
xE = [xm,xM];
yE = mE*xE+bE;
yM = 100;
% display slopes and pvaleus
pS,pE
fprintf('spon speed: m = %0.3f; pval = %0.10f\n',mS,pS)
fprintf('evok speed: m = %0.3f; pval = %0.10f\n\n',mE,pE)
% plotting
hold on;
% spon
scatter(ampS,spdS,'b.')
plot(xS,yS,'color','b','linewidth',2)
% evoked
scatter(ampE,spdE,'r.')
plot(xE,yE,'color','r','linewidth',2)
% significancy
[~,sig] = isSignificant(pS);
xt = xM*0.6; a = 0.008; yt = yM*a; 
str = sprintf('m=%0.2f, %s',mS,sig);
text(xt,yt,str,'color','b','FontSize',ylabelFontSize+2,'fontweight','bold')
[~,sig] = isSignificant(pE);
yt = yM*a*0.6; 
str = sprintf('m=%0.2f, %s',mE,sig);
text(xt,yt,str,'color','r','FontSize',ylabelFontSize+2,'fontweight','bold')

xlim([0 xM]);
ylim([0 yM]);
xlabel('VSD Signal Amplitude (\DeltaF/F_0 %)','FontSize',ylabelFontSize);
ylabel('Speed (mm/sec)','FontSize',ylabelFontSize);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
xs = 0.1*round(10*xM/4);
set(ff.ha,'XTick',0:xs:xM);
yL = round(logspace(log10(1),log10(yM),3));
set(ff.ha,'YTick',yL);
set(ff.ha,'Fontweight','bold','FontSize',axesFontSize);
set(gca,'yscale','log')
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);

function plotAmpDirection (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; 
xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;

ampS = []; dS = []; ampE = []; dE = [];
for an = length(animals):-1:1
    ampS = [ampS,valsSCC{an,selSCC}.alldf];
    dS = [dS,abs(valsSCC{an,selSCC}.allfd)];
    ampE = [ampE,valsE{an}.alldf{1}];
    dE = [dE,valsE{an}.allfd{1}];
end
idxS = ampS<=0.7943;
ampS = ampS(idxS);
dS = dS(idxS);
idxS = ampE<=0.8853;
ampE = ampE(idxS);
dE = dE(idxS);
% linear model
% spon
mdldS = fitlm(ampS,dS);
mdldS = fitlm(ampS,dS,'linear','RobustOpts','on');
Rsquared = mdldS.Rsquared.Adjusted;
pS = mdldS.anova.pValue(1);
bS = mdldS.Coefficients.Estimate(1);
mS = mdldS.Coefficients.Estimate(2);
xM = max([ampS,ampE]);
xM = 0.01*ceil(100*xM);
xS = [0,xM];
yS = mS*xS+bS;
% evoked
mdldE = fitlm(ampE,dE);
mdldE = fitlm(ampE,dE,'linear','RobustOpts','on');
Rsquared = mdldE.Rsquared.Adjusted;
pE = mdldE.anova.pValue(1);
bE = mdldE.Coefficients.Estimate(1);
mE = mdldE.Coefficients.Estimate(2);
xE = [0,xM];
yE = mE*xE+bE;
yM = 1; 
ym = -1;
% display slopes and pvaleus
pS,pE
fprintf('spon speed: m = %0.3f; pval = %0.10f\n',mS,pS)
fprintf('evok speed: m = %0.3f; pval = %0.10f\n\n',mE,pE)
% plotting
hold on;
% spon
scatter(ampS,dS,'b.')
plot(xS,yS,'color','b','linewidth',2)
% evoked
scatter(ampE,dE,'r.')
plot(xE,yE,'color','r','linewidth',2)
% significancy
[~,sig] = isSignificant(pS);
xt = xM*0.6; a = 0.2; yt = ym+(yM-ym)*a; 
str = sprintf('m=%0.2f, %s',mS,sig);
text(xt,yt,str,'color','b','FontSize',ylabelFontSize+2,'fontweight','bold')
[~,sig] = isSignificant(pE);
yt = ym+(yM-ym)*a*0.4; 
str = sprintf('m=%0.2f, %s',mE,sig);
text(xt,yt,str,'color','r','FontSize',ylabelFontSize+2,'fontweight','bold')

xlim([0 xM]);
ylim([ym yM]);
xlabel('VSD Signal Amplitude (\DeltaF/F_0 %)','FontSize',ylabelFontSize);
ylabel('Direction Stability','FontSize',ylabelFontSize);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
xs = 0.1*round(10*xM/4);
set(ff.ha,'XTick',0:xs:xM);
ys = 0.5;
set(ff.ha,'YTick',-1:ys:yM);
set(ff.ha,'Fontweight','bold','FontSize',axesFontSize);

save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);

function plotAmpFD (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; 
xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;

match = {'FL' 'HL' 'VC'};
StrajName = 'Straj_';
pathPrefix = 'Y:\homes\navvab.afrashteh\CloudStation\Navvab-Sam\NeuroPhotonic Paper\Test\';
saveStraj_FileName = sprintf('%sFigures_new_code\\%s%s_allAn_04corr.mat',pathPrefix ,StrajName,match{ss});
load(saveStraj_FileName,'Straj')
EtrajName = 'Etraj_';
saveEtraj_FileName = sprintf('%sFigures_new_code\\%s%s_allAn_04corr.mat',pathPrefix ,EtrajName,match{ss});
load(saveEtraj_FileName,'Etraj')

ampS = []; fdS = []; ampE = []; fdE = [];
for an = length(animals):-1:1
    mask  = mdataE{ss}.masks{an}{1};
    % spon
    Lscc = ismember(valsSCC{an,selSCC}.allFrames(:,[2,5]),Straj{an}.TrajFrames(:,[2,5]),'row');
    Ltrj = ismember(Straj{an}.TrajFrames(:,[2,5]),valsSCC{an,selSCC}.allFrames(:,[2,5]),'row');
    sX = Straj{an}.Xinterp(Ltrj,:);
    sY = Straj{an}.Yinterp(Ltrj,:);
    sigma = 2;
    FD = calcFD(sX,sY,mask,sigma);
    fdS = [fdS,FD];
    ampS = [ampS,valsSCC{an,selSCC}.alldf(Lscc)];
    % evok
    ii = 1; % evoked level
    Le = ismember(valsE{an}.Trials(ii,:),Etraj{an}.Trials{ii});
    Let = ismember(Etraj{an}.Trials{ii},valsE{an}.Trials(ii,:));
    eX = Etraj{an}.Xinterp{ii}(Let,:);
    eY = Etraj{an}.Yinterp{ii}(Let,:);
    FD = calcFD(eX,eY,mask,sigma);
    fdE = [fdE,FD];
    Le = Le(1:length(valsE{an}.alldf{ii}));
    ampE = [ampE,valsE{an}.alldf{ii}(Le)];
end
% linear model
% spon
mdlFDS = fitlm(ampS,fdS);
mdlFDS = fitlm(ampS,fdS,'linear','RobustOpts','on');
Rsquared = mdlFDS.Rsquared.Adjusted;
pS = mdlFDS.anova.pValue(1);
bS = mdlFDS.Coefficients.Estimate(1);
mS = mdlFDS.Coefficients.Estimate(2);
xM = max([ampE,ampS]);
xM = 0.01*ceil(xM*100);
xS = [0,xM];
yS = mS*xS+bS;
% evoked
mdlFDE = fitlm(ampE,fdE);
mdlFDE = fitlm(ampE,fdE,'linear','RobustOpts','on');
RsquaredE = mdlFDE.Rsquared.Adjusted;
pE = mdlFDE.anova.pValue(1);
bE = mdlFDE.Coefficients.Estimate(1);
mE = mdlFDE.Coefficients.Estimate(2);
xE = [0,xM];
yE = mE*xE+bE;
yM = 1.6;
ym = 0.8;
% display slopes and pvaleus
pS,pE
fprintf('spon speed: m = %0.3f; pval = %0.10f\n',mS,pS)
fprintf('evok speed: m = %0.3f; pval = %0.10f\n\n',mE,pE)
% plotting
hold on;
% spon
scatter(ampS,fdS,'b.')
plot(xS,yS,'color','b','linewidth',2)
% evoked
scatter(ampE,fdE,'r.')
plot(xE,yE,'color','r','linewidth',2)
% significancy
[~,sig] = isSignificant(pS);
xt = xM*0.6; a = 0.2; yt = ym+(yM-ym)*a; 
str = sprintf('m=%0.2f, %s',mS,sig);
text(xt,yt,str,'color','b','FontSize',ylabelFontSize+2,'fontweight','bold')
[~,sig] = isSignificant(pE);
yt = ym+(yM-ym)*a*0.4;  
str = sprintf('m=%0.2f, %s',mE,sig);
text(xt,yt,str,'color','r','FontSize',ylabelFontSize+2,'fontweight','bold')

xlim([0 xM]);
ylim([ym yM]);
xlabel('VSD Signal Amplitude (\DeltaF/F_0 %)','FontSize',ylabelFontSize);
ylabel('Fractal Dimension','FontSize',ylabelFontSize);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
xs = 0.1*round(10*xM/4);
set(ff.ha,'XTick',0:xs:xM);
ys = 0.3;
set(ff.ha,'YTick',1:ys:yM);
set(ff.ha,'Fontweight','bold','FontSize',axesFontSize);
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);

function FD = calcFD(sX,sY,mask,sigma)
% calculate fractal dimension for each trajectory line in (sX,sY)
[dim1,dim2] = size(mask);
method = 'bicubic';
param.Sigma = sigma;
param.GaussWindow = 4*param.Sigma+1;
Center  = fix([param.GaussWindow/2,param.GaussWindow/2])+1;
[R,C]   = ndgrid(1:param.GaussWindow, 1:param.GaussWindow);
Gauss2d = gauss2dC(R,C,param.Sigma,Center);
for jj = 1:size(sX,1)
    I = zeros(dim1, dim2);
    xi = sX(jj,:);
    yi = sY(jj,:);
    xd = diff(xi);
    yd = diff(yi);
    TrajLen = sqrt(xd.^2 + yd.^2);
    TrajLen = sum(TrajLen);
    np = round(TrajLen*3);
    [cx,cy,~] = improfile(I,xi,yi,np,method);
    cx = round(cx); cy = round(cy);
    I((cx-1)*dim1+cy) = 1;
    I = conv2(I,Gauss2d,'same');
    I = I.*mask;
    level = graythresh(I);
    imbw = im2bw(I,level);
    FD(jj) = BoxCountfracDim(imbw);
end