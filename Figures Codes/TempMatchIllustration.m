function TempMatchIllustration()
linewidth = 1.5; fontsize = 10;
Cmap = evalin('base','Cmap');
fde = evalin('base','mdataE');
thn = 1;
CrossOrMax = 'Max';
LMH = 3;
multLevels = 1;

an = 1;
listR = listOfRecordingsSpon(fde{1}.animalNumber(an));
Mask = fde{an}.masks{1}{1};
ii = 6;
figure(1); clf;
set(gcf,'units','inches','pos',[1,1,6.9,4],'color','w');

ax1 = axes('units','inches','pos',[0.1,2.9+0.55,3,0.5]);
ax2 = axes('units','inches','pos',[0.1,2.9,3,0.5]);
ax3 = axes('units','inches','pos',[3.5,2.9,3,1]);
ax4 = axes('units','inches','pos',[3.2,2.9-0.95,3.6,0.5]);
ax5 = axes('units','inches','pos',[3.2,2.9-1.5,3.6,0.5]);
ax6 = axes('units','inches','pos',[0.5,1.9-1.3,6,0.5]);
ax7 = axes('units','inches','pos',[0.5,1.9-1.35,6,0.5]);
ax8 = axes('units','inches','pos',[1.3,1.6,1,1]);
hold on;
for ss = 1:2
    if ss==1
        match = 'FL'; col = 'r';
    else
        match = 'HL'; col = 'b';
    end
    [~, names] = getEvokedListsDiffStimAmps(fde{ss}.animalNumber(an),multLevels);
    ssLMH = cell2mat(fde{ss}.anGroup{an});
    LMHname = names{ss}{ssLMH(LMH)};
    motifs = getMotifsFromCCResultsSam_simple(fde{ss}.animalNumber(an),listR,ii,thn,match,CrossOrMax,LMHname,multLevels);
    axes(ax3); hold on;
    plot(motifs.allCCs,'color',col,'linewidth',linewidth)
    % get template
    anPath = makeName(fde{ss}.animalNumber(an),getMainDataFolder);
    bregmaPath = makeName('bregmaXY.txt',anPath);
    bregma = load(bregmaPath);
    pEvoked = makeName('pEvoked',anPath);
    folPath = makeName(LMHname,pEvoked);
    tempPath = makeName('stimTemplate.mat',folPath);
    stimTemplate = load(tempPath);
    stimTemplate = stimTemplate.stimTemplate;
    frameNumbers = stimTemplate.frameNumbers;
    % add one frame before and one after just for display purposes
    frameNumbers = [frameNumbers(1)-1, frameNumbers, frameNumbers(end)+1];
    ImgSeqPath = makeName('ImgSeq.mat',folPath);
    ImgSeq = load(ImgSeqPath);
    ImgSeq = ImgSeq.ImgSeq;
    tempSeq = ImgSeq(:,:,frameNumbers);
    tempSeq = applyMask (tempSeq,Mask); % apply mask
    tempSeq = tempSeq(:);
    tempSeq = reshape(tempSeq,128,128*5); % put frames besides each other
    axes(ax1)
    if ss==2
        axes(ax2)
    end
    imagesc(tempSeq,[-0.2,0.7]); hold on
    plot(bregma(1)+(0:128:128*4),bregma(2),'o','MarkerSize',2,'MarkerEdgeColor','w','MarkerFaceColor','w')
    colormap(Cmap)
    axis image
    set(gca,'xtick',[],'ytick',[])
    % load spon image sequence and display example matches
    ImgSeqPathS = makeName([listR{ii},'.mat'],anPath);
    matObj = matfile(ImgSeqPathS);
    exFrame = 4;
    f1 = motifs.frames(exFrame)-6;
    if ss==2
        exFrame = 1;
        f1 = motifs.frames(exFrame)-13;
    end
    f2 = f1+16;
    exDF = matObj.DF_F0(:,:,f1:3:f2);
    exDF = applyMask (exDF,Mask); % apply mask
    exDF = exDF(:);
    exDF = reshape(exDF,128,128*6); % put frames besides each other
    axes(ax4)
    if ss==2
        axes(ax5)
    end
    imagesc(exDF,[0,0.3]); hold on
    plot(bregma(1)+(0:128:128*5),bregma(2),'o','MarkerSize',2,'MarkerEdgeColor','w','MarkerFaceColor','w')
    colormap(Cmap)
    axis image
    set(gca,'xtick',[],'ytick',[])
    % save CCs for dispay
    m = 250;
    cc = motifs.allCCs(400:4899);
    cc = repmat(cc,m,1);
    CCs((ss-1)*m+1:ss*m,:) = cc;
end
a = 10;
CCs(m-a:m+a+1,:) = -1;
axes(ax3)
xm = 400; xM = xm+250; xd = 40;
th = 0.48;
xlim([xm,xM])
plot(xm:xM,th*ones(1,xM-xm+1),'k--')
set(gca,'xtick',(xm:xd:xM));
set(gca,'xticklabel',(0:xd:xM-xm)*1000/100);
ym = -0.6; yM = -ym; yd = 0.2;
ylim([ym,yM]);
set(gca,'ytick',(ym+0.0:yd:yM));
set(gca,'yticklabel',(ym+0.0:yd:yM));
xlabel('Time (ms)','fontsize',fontsize)
ylabel('Correlation index','fontsize',fontsize)
set(ax3,'ycolor','k','xcolor','k','tickdir','out')
txtStr = sprintf('Threshold = %0.2f',th);
text(xm*1.01,th*1.2,txtStr,'fontsize',fontsize,'color','k')
% display correlation values for all sequence
axes(ax7)
set(ax7,'xtick',linspace(0,1,7),'xticklabel',linspace(0,30,7))
set(ax7,'ycolor','w','xcolor','k','tickdir','out')
xlabel('Time (s)','fontsize',fontsize)
axes(ax6)
imagesc(CCs,[ym,yM]); hold on
colormap(Cmap)
set(gca,'xtick',[],'ytick',[])
%% example frame for spon
axes(ax8)
ff = 20;
exFrame = matObj.DF_F0(:,:,ff);
exFrame = applyMask (exFrame,Mask); % apply mask
imagesc(exFrame,[0.0,0.6]); hold on
plot(bregma(1),bregma(2),'o','MarkerSize',4,'MarkerEdgeColor','w','MarkerFaceColor','w')
colormap(Cmap)
axis image
set(gca,'xtick',[],'ytick',[],'ydir','rev')
%%
% save to pdf
pdfFileName = 'Template Matching Illustration Matlab.pdf';
pdfFileName = makeName(pdfFileName,getpdffolder);
save2pdf(pdfFileName,gcf,600)
n=0;