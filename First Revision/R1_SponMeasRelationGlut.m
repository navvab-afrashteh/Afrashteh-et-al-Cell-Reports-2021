function R1_SponMeasRelationGlut(roi,selE)
close all
clc

if ~exist('roi')
    roi = 1; selE = [1,2,3];
end
match = {'Glut'};
ss = 1;
mfn = mfilename;
pdfFileName = sprintf('%s_%s_roi_%d.pdf',mfn(1:end-4),match{ss},roi);
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
dfth = 0;
CCTh = 1;
cols = getColors;
selSCC = 3;
amplitude_threshold = meta.amplitude_threshold;
% mfilename
fileName = sprintf('%s_ss%d_roi%d.mat',mfilename,ss+3,roi);
rewrite = 0;
if rewrite
    for an = 1:length(animals)
        valsE{an} = getDistributions_simpleAmps_artur_2_glut(mdataE{ss},animals(an),roi,ss,amplitude_threshold);
        valsS{an} = getDistributionsSpon_simple_artur_glut(mdataS{1},animals(an),roi,ss,dfth,amplitude_threshold);
        for LMH = 3%length(ssLMH)
            valsSCC{an,LMH} = getDistributionsSpon_motifs_simple_artur_glut_2 (mdataE,mdataS{1},animals(an),roi,'AC','Max',CCTh,LMH,amplitude_threshold);
        end
    end
    save(fileName,'valsE','valsS','valsSCC');
else
    load(fileName,'valsE','valsS','valsSCC');
end

% Amplitude and Speed relation
ff = makeFigureWindow__one_axes_only(2,[1 4 1.75 1.5],[0.21 0.25 0.7 0.69]);
plotAmpSpeed (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta,selE);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder));

% Amplitude and directionality relation
ff = makeFigureWindow__one_axes_only(2,[1 4 1.75 1.5],[0.21 0.25 0.7 0.69]);
plotAmpDirection (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta,selE);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder));

% Amplitude and Fractal Dimension
ff = makeFigureWindow__one_axes_only(2,[1 4 1.75 1.5],[0.21 0.25 0.7 0.69]);
plotAmpFD (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta,selE);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder));


function plotAmpSpeed (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta,selE)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; 
xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;

ampS = []; spdS = []; ampE = []; spdE = [];
for an = length(animals):-1:1
    ampS = [ampS,valsSCC{an,selSCC}.alldf];
    spdS = [spdS,valsSCC{an,selSCC}.absalluv];
    for iii = selE
        ii = selLMH_all{ss}(an,iii);
        ampE = [ampE,valsE{an}.alldf{ii}];
        spdE = [spdE,valsE{an}.absalluv{ii}];
    end
end
% spdS = log10(spdS);
% spdE = log10(spdE);
ei = linspace(0,25,100);
[cdfS,ci] = CDF(spdS,ei);
cdfE = CDF(spdE,ei);
ci = [0,ci]; cdfS = [0,cdfS]; cdfE = [0,cdfE];
% linear model
%% bootstraping to get a robust relation
Nb = 100; % number of bootstrapping
f = 0.5;  % fraction of original data set used in each iteration
rng(0);   % for reproducibility
% spon
n = length(spdS); k = round(n*f);
Ke = permuteMat(n,k,Nb);
pS = []; bS = []; mS = []; RsquaredS = [];
for ni = Nb:-1:1
    mdlSpeedS = fitlm(ampS(Ke(ni,:)),spdS(Ke(ni,:)));
    RsquaredS(ni) = mdlSpeedS.Rsquared.Adjusted;
    pS(ni) = mdlSpeedS.anova.pValue(1);
    bS(ni) = mdlSpeedS.Coefficients.Estimate(1);
    mS(ni) = mdlSpeedS.Coefficients.Estimate(2);
end
bS = median(bS); mS = median(mS); RsquaredS = median(RsquaredS);
xM = max([ampS,ampE]);
xM = 0.01*ceil(100*xM);
xS = [0,xM];
yS = mS*xS+bS;
% evoked
n = length(spdE); k = round(n*f);
Ke = permuteMat(n,k,Nb);
pE = []; bE = []; mE = []; RsquaredE = [];
for ni = Nb:-1:1
    mdlSpeedE = fitlm(ampE(Ke(ni,:)),spdE(Ke(ni,:)));
    RsquaredE(ni) = mdlSpeedE.Rsquared.Adjusted;
    pE(ni) = mdlSpeedE.anova.pValue(1);
    bE(ni) = mdlSpeedE.Coefficients.Estimate(1);
    mE(ni) = mdlSpeedE.Coefficients.Estimate(2);
end
bE = median(bE); mE = median(mE); RsquaredE = median(RsquaredE);
xE = [0,xM];
yE = mE*xE+bE;
% display slopes and pvaleus
pS = median(pS); pE = median(pE);
fprintf('spon speed: m = %0.3f; pval = %0.33f; R2 = %0.2f\n',mS,pS,RsquaredS)
fprintf('evok speed: m = %0.3f; pval = %0.15f; R2 = %0.2f\n\n',mE,pE,RsquaredE)
% plotting
hold on;
scatter(ampS,spdS,'b.')
scatter(ampE,spdE,'r.')
plot(xS,yS,'color','b','linewidth',2)
plot(xE,yE,'color','r','linewidth',2)
yM = 25; ym = -1; xM = 3.5;
% significancy
[~,sig] = isSignificant(pS);
xt = xM*0.55;
a = 0.2; if ss==2; a = 0.033; end
yt = (yM-ym)*a+ym; 
str = sprintf('m=%0.2f, %s',mS,sig);
text(xt,yt,str,'color','b','FontSize',ylabelFontSize+2,'fontweight','bold')
[~,sig] = isSignificant(pE);
yt = (yM-ym)*a*0.4+ym; 
str = sprintf('m=%0.2f, %s',mE,sig);
text(xt,yt,str,'color','r','FontSize',ylabelFontSize+2,'fontweight','bold')

xlim([0 xM]);
ylim([ym yM]);
xlabel('Glu Signal Amplitude (\DeltaF/F_0 %)','FontSize',ylabelFontSize);
hy = ylabel('Speed (mm/sec)','FontSize',ylabelFontSize);
pos = get(hy,'Position');pos(1) = pos(1) - 0.0;set(hy,'Position',pos);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
xs = 0.1*round(10*xM/4);
xL = 0:xs:xM; xL = xL(1:end-1);
set(ff.ha,'XTick',xL);
yL = 0:5:25;
set(ff.ha,'YTick',yL);
set(ff.ha,'Fontweight','bold','FontSize',axesFontSize);
% CDFs
hcdf = axes('Position',[0.72,0.7,0.2,0.2]);
hold on
plot(ci,cdfE,'color','r','linewidth',1)
plot(ci,cdfS,'color','b','linewidth',1)
xM = round(max(ci));
xlabel('Speed','FontSize',ylabelFontSize);
hy = ylabel('CDF','FontSize',ylabelFontSize);
pos = get(hy,'Position');pos(1) = pos(1) + 10;set(hy,'Position',pos);
set(hcdf,'TickDir','out','FontSize',axesFontSize);
xs = 10;
set(hcdf,'XTick',0:xs:xM);
yL = 0:1;
set(hcdf,'YTick',yL);
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);

function plotAmpDirection (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta,selE)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; 
xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;

ampS = []; dS = []; ampE = []; dE = [];
for an = length(animals):-1:1
    ampS = [ampS,valsSCC{an,selSCC}.alldf];
    dS = [dS,abs(valsSCC{an,selSCC}.allfd)];  
    for iii = selE
        ii = selLMH_all{ss}(an,iii);
        ampE = [ampE,valsE{an}.alldf{ii}];
        dE = [dE,valsE{an}.allfd{ii}];
    end
end
ei = linspace(-1,1,100);
[cdfS,ci] = CDF(dS,ei);
cdfE = CDF(dE,ei);
ci = [-1,ci,1]; cdfS = [0,cdfS,1]; cdfE = [0,cdfE,1];
%% bootstraping to get a robust relation
Nb = 100; % number of bootstrapping
f = 0.75;  % fraction of original data set used in each iteration
rng(0);   % for reproducibility
% spon
n = length(dS); k = round(n*f);
Ke = permuteMat(n,k,Nb);
pS = []; bS = []; mS = []; RsquaredS = [];
for ni = Nb:-1:1
    mdlSpeedS = fitlm(ampS(Ke(ni,:)),dS(Ke(ni,:)));
    RsquaredS(ni) = mdlSpeedS.Rsquared.Adjusted;
    pS(ni) = mdlSpeedS.anova.pValue(1);
    bS(ni) = mdlSpeedS.Coefficients.Estimate(1);
    mS(ni) = mdlSpeedS.Coefficients.Estimate(2);
end
bS = median(bS); mS = median(mS); RsquaredS = median(RsquaredS);
xM = max([ampS,ampE]);
xM = 0.01*ceil(100*xM);
xS = [0,xM];
yS = mS*xS+bS;
% evoked
n = length(dE); k = round(n*f);
Ke = permuteMat(n,k,Nb);
pE = []; bE = []; mE = []; RsquaredE = [];
for ni = Nb:-1:1
    mdlSpeedE = fitlm(ampE(Ke(ni,:)),dE(Ke(ni,:)));
    RsquaredE(ni) = mdlSpeedE.Rsquared.Adjusted;
    pE(ni) = mdlSpeedE.anova.pValue(1);
    bE(ni) = mdlSpeedE.Coefficients.Estimate(1);
    mE(ni) = mdlSpeedE.Coefficients.Estimate(2);
end
bE = median(bE); mE = median(mE); RsquaredE = median(RsquaredE);
xE = [0,xM];
yE = mE*xE+bE;
% display slopes and pvaleus
pS = median(pS); pE = median(pE);
fprintf('spon dir stability: m = %0.3f; pval = %0.10f; R2 = %0.4f\n',mS,pS,RsquaredS)
fprintf('evok dir stability: m = %0.3f; pval = %0.15f; R2 = %0.4f\n\n',mE,pE,RsquaredE)

% plotting
hold on;
scatter(ampS,dS,'b.')
scatter(ampE,dE,'r.')
plot(xS,yS,'color','b','linewidth',2)
plot(xE,yE,'color','r','linewidth',2)
yM = 1; ym = -1; xM = 2.6;
% significancy
[~,sig] = isSignificant(pS);
xt = xM*0.04;
a = 0.2; if ss==2; a = 0.033; end
yt = (yM-ym)*a+ym; 
str = sprintf('m=%0.2f, %s',mS,sig);
text(xt,yt,str,'color','b','FontSize',ylabelFontSize+2,'fontweight','bold')
[~,sig] = isSignificant(pE);
yt = (yM-ym)*a*0.4+ym; 
str = sprintf('m=%0.2f, %s',mE,sig);
text(xt,yt,str,'color','r','FontSize',ylabelFontSize+2,'fontweight','bold')

xlim([0 xM]);
ylim([ym yM]);
xlabel('Glu Signal Amplitude (\DeltaF/F_0 %)','FontSize',ylabelFontSize);
hy = ylabel('Direction Stability','FontSize',ylabelFontSize);
pos = get(hy,'Position');pos(1) = pos(1) - 0.0;set(hy,'Position',pos);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
xs = 0.1*round(10*xM/4);
xL = 0:xs:xM; %xL = xL(1:end-1);
set(ff.ha,'XTick',xL);
yL = -0.5:0.5:1;
set(ff.ha,'YTick',yL);
set(ff.ha,'Fontweight','bold','FontSize',axesFontSize);

% CDFs
hcdf = axes('Position',[0.68,0.42,0.2,0.2]);
hold on
plot(ci,cdfE,'color','r','linewidth',1)
plot(ci,cdfS,'color','b','linewidth',1)
xlabel('Direction Stability','FontSize',ylabelFontSize);
hy = ylabel('CDF','FontSize',ylabelFontSize);
pos = get(hy,'Position');pos(1) = pos(1) + 0.7;set(hy,'Position',pos);
set(hcdf,'TickDir','out','FontSize',axesFontSize);
xs = 1;
set(hcdf,'XTick',-1:xs:1);
yL = 0:1;
set(hcdf,'YTick',yL);
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);

function plotAmpFD (ff,animals,mdataE,valsE,valsS,valsSCC,selSCC,ss,selLMH_all,cols,meta,selE)
% figure related settings
annotationFontSize = 8; sF = 1; axesFontSize = 6; legendFontSize = 5; 
xlabelFontSize = 6; ylabelFontSize = 6;
titleFontSize = 8;

match = {'AC'};
pathPrefix = 'Y:\homes\navvab.afrashteh\CloudStation\Navvab-Sam\NeuroPhotonic Paper\Test\';
StrajName = 'Straj_BPF_05_6Hz_';
saveStraj_FileName = sprintf('%sFigures_new_code_AC\\%s%s_glut_allAn_04corr.mat',pathPrefix ,StrajName,match{ss});
load(saveStraj_FileName,'Straj')
EtrajName = 'Etraj_BPF_05_6Hz_';
saveEtraj_FileName = sprintf('%sFigures_new_code_AC\\%s%s_glut_allAn_04corr.mat',pathPrefix ,EtrajName,match{ss});
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
    for iii = selE
        ii = selLMH_all{ss}(an,iii); % evoked level
        Le = ismember(valsE{an}.Trials(ii,:),Etraj{an}.Trials{ii});
        Let = ismember(Etraj{an}.Trials{ii},valsE{an}.Trials(ii,:));
        eX = Etraj{an}.Xinterp{ii}(Let,:);
        eY = Etraj{an}.Yinterp{ii}(Let,:);
        FD = calcFD(eX,eY,mask,sigma);
        fdE = [fdE,FD];
        Le = Le(1:length(valsE{an}.alldf{ii}));
        ampE = [ampE,valsE{an}.alldf{ii}(Le)];
    end
end
thL = mean(ampS)-3*std(ampS); thH = mean(ampS)+3*std(ampS); 
idxA = (ampS >= thL) & (ampS <= thH);
ampS = ampS(idxA);
fdS = fdS(idxA);
thL = mean(ampE)-3*std(ampE); thH = mean(ampE)+3*std(ampE); 
idxA = (ampE >= thL) & (ampE <= thH);
ampE = ampE(idxA);
fdE = fdE(idxA);

% CDFs
ei = linspace(0.8,1.6,100);
[cdfS,ci] = CDF(fdS,ei);
cdfE = CDF(fdE,ei);
ci = [0.8,ci,1.6]; cdfS = [0,cdfS,1]; cdfE = [0,cdfE,1];
% linear model
%% bootstraping to get a robust relation
Nb = 100; % number of bootstrapping
f = 0.75;  % fraction of original data set used in each iteration
rng(0);   % for reproducibility
% spon
n = length(fdS); k = round(n*f);
Ke = permuteMat(n,k,Nb);
pS = []; bS = []; mS = []; RsquaredS = [];
for ni = Nb:-1:1
    mdlSpeedS = fitlm(ampS(Ke(ni,:)),fdS(Ke(ni,:)));
    RsquaredS(ni) = mdlSpeedS.Rsquared.Adjusted;
    pS(ni) = mdlSpeedS.anova.pValue(1);
    bS(ni) = mdlSpeedS.Coefficients.Estimate(1);
    mS(ni) = mdlSpeedS.Coefficients.Estimate(2);
end
bS = median(bS); mS = median(mS); RsquaredS = median(RsquaredS);
xM = max([ampS,ampE]);
xM = 0.01*ceil(100*xM);
xS = [0,xM];
yS = mS*xS+bS;
% evoked
n = length(fdE); k = round(n*f);
Ke = permuteMat(n,k,Nb);
pE = []; bE = []; mE = []; RsquaredE = [];
for ni = Nb:-1:1
    mdlSpeedE = fitlm(ampE(Ke(ni,:)),fdE(Ke(ni,:)));
    RsquaredE(ni) = mdlSpeedE.Rsquared.Adjusted;
    pE(ni) = mdlSpeedE.anova.pValue(1);
    bE(ni) = mdlSpeedE.Coefficients.Estimate(1);
    mE(ni) = mdlSpeedE.Coefficients.Estimate(2);
end
bE = median(bE); mE = median(mE); RsquaredE = median(RsquaredE);
xE = [0,xM];
yE = mE*xE+bE;
% display slopes and pvaleus
pS = median(pS); pE = median(pE);
fprintf('spon fractal dim: m = %0.3f; pval = %0.10f; R2 = %0.2f\n',mS,pS,RsquaredS)
fprintf('evok fractal dim: m = %0.3f; pval = %0.15f; R2 = %0.2f\n\n',mE,pE,RsquaredE)

% plotting
hold on;
scatter(ampS,fdS,'b.')
scatter(ampE,fdE,'r.')
plot(xS,yS,'color','b','linewidth',2)
plot(xE,yE,'color','r','linewidth',2)
yM = 1.6; ym = 0.8; xM = 3.2;
% significancy
[~,sig] = isSignificant(pS);
xt = xM*0.5;
a = 0.2; if ss==2; a = 0.033; end
yt = (yM-ym)*a+ym; 
str = sprintf('m=%0.2f, %s',mS,sig);
text(xt,yt,str,'color','b','FontSize',ylabelFontSize+2,'fontweight','bold')
[~,sig] = isSignificant(pE);
yt = (yM-ym)*a*0.4+ym; 
str = sprintf('m=%0.2f, %s',mE,sig);
text(xt,yt,str,'color','r','FontSize',ylabelFontSize+2,'fontweight','bold')

xlim([0 xM]);
ylim([ym yM]);
xlabel('Glu Signal Amplitude (\DeltaF/F_0 %)','FontSize',ylabelFontSize);
hy = ylabel('Fractal Dimension','FontSize',ylabelFontSize);
pos = get(hy,'Position');pos(1) = pos(1) - 0.0;set(hy,'Position',pos);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
xs = 0.1*round(10*xM/4);
xL = 0:xs:xM; xL = xL(1:end-1);
set(ff.ha,'XTick',xL);
yL = 1:0.3:1.6;
set(ff.ha,'YTick',yL);
set(ff.ha,'Fontweight','bold','FontSize',axesFontSize);

% CDFs
hcdf = axes('Position',[0.67,0.7,0.2,0.2]);
hold on
plot(ci,cdfE,'color','r','linewidth',1)
plot(ci,cdfS,'color','b','linewidth',1)
xM = (max(ci));
xlabel('Fractal Dimension','FontSize',ylabelFontSize);
hy = ylabel('CDF','FontSize',ylabelFontSize);
pos = get(hy,'Position');pos(1) = pos(1) + 0.3;set(hy,'Position',pos);
set(hcdf,'TickDir','out','FontSize',axesFontSize);
xs = 0.5;
set(hcdf,'XTick',1:xs:xM);
yL = 0:1;
set(hcdf,'YTick',yL);
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

function Ke = permuteMat(n,k,N)
% create permutation matrix
% k = 20; N = 500;
Ke = zeros(N,k);
i = 0;
while i<N
    Ke(i+1,:) = sort(randperm(n,k));
    i = i + 1;
    if i==N
        Ke = unique(Ke,'rows','stable');
        i = size(Ke,1);
        Ke = [Ke;zeros(N-i,k)];
    end
end

function [cdf,ci] = CDF(x,ei)
cdf = histogram(x,ei,'Visible','off'); cdf = cdf.Values;
cdf = cumsum(cdf);
cdf = cdf/max(cdf);
ci = (ei(1:end-1)+ei(2:end))/2;