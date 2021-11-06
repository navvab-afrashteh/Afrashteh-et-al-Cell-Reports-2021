function R1_SponPCCEvent_TimePercentage
clc;close all;

pdfFileName = sprintf('%s.pdf',mfilename);
pdfFileName = makeName(pdfFileName,getpdffolder)
if exist(pdfFileName,'file')
    delete(pdfFileName);
end
if exist(makeName('temp.pdf',getpdffolder),'file')
    delete(makeName('temp.pdf',getpdffolder));
end

stim = {'FL' 'HL','VC','AC'};
folName = 'Y:\homes\navvab.afrashteh\CloudStation\Navvab-Sam\NeuroPhotonic Paper\Test\First Revision';
for ss = 1:length(stim)
    match = stim{ss};
    fileName = sprintf('%s\\sponDurPercentage_%s',folName,match);
    load(fileName,'sponDurPercentage')
    eval(['tp',match, '= sponDurPercentage;']);
end

% stat on event rate and plotting bar graph
annovaVar = [tpFL, tpHL, tpVC, tpAC];
g1 = {'FL','FL','FL','FL','FL',...
      'HL','HL','HL','HL','HL',...
      'VC','VC','VC','VC',...
      'AC','AC','AC','AC','AC'};
g2 = {'VSD1','VSD1','VSD1','VSD1','VSD1',...
      'VSD1','VSD1','VSD1','VSD1','VSD1',...
      'VSD2','VSD2','VSD2','VSD2',...
      'Glut','Glut','Glut','Glut','Glut'};
% p = anovan(annovaVar,{g1,g2}); % in case you want to check for animal groups
thisCols = {'k','k','k','k'};
thisCols = {'b','r','g','k'};
[p,tbl,stats] = anova1(annovaVar,g1,'off');
[c,~,~,gnames] = multcompare(stats,'CType','hsd');
prT = c(:,6);
hrT = prT<0.05;
mVals = [mean(tpFL), mean(tpHL), mean(tpVC(1:3)), mean(tpAC)];
semVals = [std(tpFL)/sqrt(5), std(tpHL)/sqrt(5), std(tpVC)/sqrt(4), std(tpAC)/sqrt(5)];
combs = nchoosek(1:size(mVals,2),2);

etta2 = tbl{2,2}/tbl{4,2};
fprintf('F-value = %0.2f; etta2 = %0.2f; p-value = %f \n',tbl{2, 5},etta2,p)
fprintf('post-hoc comparisons: \n')
for ii = 1:length(prT)
    fprintf('%s vs %s : p = %f \n',gnames{combs(ii,2)}, gnames{combs(ii,1)}, prT(ii))
end

close all
axesFontSize = 6;
ff = makeFigureWindow__one_axes_only(4,[1 1 1.5 1.75],[0.25 0.2 0.7 0.68]);
plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols,'ySpacingFactor',10);
xlim([0.4 0.6+size(mVals,2)]);
hy = ylabel('PCC Events Time Percentage (%)');
pos = get(hy,'Position');pos(1) = pos(1) + 0.0;set(hy,'Position',pos);
set(ff.ha,'TickDir','out','FontSize',axesFontSize);
set(ff.ha,'XTick',1:size(mVals,2));
set(ff.ha,'XTickLabel',{'FL','HL', 'VC','AC'});
set(ff.ha,'Fontweight','bold');
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);
append_pdfs(pdfFileName,makeName('temp.pdf',getpdffolder))
