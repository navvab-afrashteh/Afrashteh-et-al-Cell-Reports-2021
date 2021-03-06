function R1_EvokedPlot_GFP ()
close all
clc

[avgDF_Glu,GFPall,GFP2,GFP1,t_GFP,t_Glut] = calcAvgs;

pdfFileName = sprintf('%s1_vs_GFP2.pdf',mfilename)
pdfFileName = makeName(pdfFileName,getpdffolder);
ff = makeFigureWindow__one_axes_only(2,[1 4 2 1.75],[0.25 0.2 0.6 0.69]);
plotHiLo(ff,GFP2,GFP1,t_GFP);
movefile(makeName('temp.pdf',getpdffolder),pdfFileName)

pdfFileName = sprintf('%s_vs_GLut.pdf',mfilename)
pdfFileName = makeName(pdfFileName,getpdffolder);
ff = makeFigureWindow__one_axes_only(2,[1 4 2 1.75],[0.25 0.2 0.6 0.69]);
plotGFPvsGlu(ff,avgDF_Glu,GFPall,t_GFP,t_Glut);
movefile(makeName('temp.pdf',getpdffolder),pdfFileName)

function [avgDF_Glu,GFPall,GFP2,GFP1,t_GFP,t_Glut] = calcAvgs
% calculate averages for GFP blocks 1 and 2 and Glut Hi stim.
reWrite = 0;
fileName = sprintf('%s.mat',mfilename);
if reWrite
    Fpass1 = 0.5; Fpass2 = 6;
    Fpass1 = round(Fpass1*100);
    Fpass2 = round(Fpass2);
    % GFP animals
    mainDataFolder = getMainDataFolderYFP;
    animalNumber =  {'190533'; '190545'};
    for an = 1:length(animalNumber)
        dataFolder = makeName(animalNumber{an},mainDataFolder);
        fileNameROI = makeName('ROI_BC.mat',dataFolder);
        maskBC = load(fileNameROI,'ROI'); maskBC = maskBC.ROI;
        mask = maskBC;
        peDataFolder = makeName('pEvoked',dataFolder);
        ignoreFolders = [];
        ofolders = dir(peDataFolder);
        folders = cleanFolderList(ofolders,ignoreFolders);
        for ii = 1%:length(folders)
            meFolder = makeName(folders(ii).name,peDataFolder);
            badTrials = [];
            osFolders = dir(meFolder);
            sFolders = cleanFolderList(osFolders,badTrials);
            for jj = length(sFolders):-1:1
                trialFolder = makeName(sFolders(jj).name,meFolder);
                fileNameDF = sprintf('%s\\ImgSeq_BPF_0%dHz_%dHz.mat',trialFolder, Fpass1, Fpass2);
                mFile = matfile(fileNameDF);
                ImgSeq = mFile.ImgSeq;
                maskdf = getMaskedValues(ImgSeq,mask);
                allDFb(jj,:,an) = mean(maskdf,2);
            end
        end
    end
    % Glut animals
    mainDataFolder = getMainDataFolderGlut;
    animalNumber =  {'170193'; '170196'; '163287'; '163146'; '163388' };
    for an = length(animalNumber):-1:1
        dataFolder = makeName(animalNumber{an},mainDataFolder);
        fileNameROI = makeName('ROI_AC.mat',dataFolder);
        maskAC = load(fileNameROI,'ROI'); maskAC = maskAC.ROI;
        mask = maskAC;
        peDataFolder = makeName('pEvoked',dataFolder);
        ignoreFolders = [];
        ofolders = dir(peDataFolder);
        folders = cleanFolderList(ofolders,ignoreFolders);
        ii = contains({folders.name},'HI');
        meFolder = makeName(folders(ii).name,peDataFolder);
        fileNameDF = sprintf('%s\\ImgSeq_AllFrames_BPF_0%dHz_%dHz.mat',meFolder, Fpass1, Fpass2);
        load(fileNameDF,'ImgSeq');
        maskdf = getMaskedValues(ImgSeq,mask);
        avgDF_Glu(an,:) = mean(maskdf,2);
    end
    % saving files
    save(fileName,'allDFb','avgDF_Glu')
else
    load(fileName,'allDFb','avgDF_Glu')
end
% GFP animals
mainDataFolder = getMainDataFolderYFP;
animalNumber =  {'190533'; '190545'};
for an = length(animalNumber):-1:1
    dataFolder = makeName(animalNumber{an},mainDataFolder);
    peDataFolder = makeName('pEvoked',dataFolder);
    ignoreFolders = [];
    ofolders = dir(peDataFolder);
    folders = cleanFolderList(ofolders,ignoreFolders);
    for ii = 1%:length(folders)
        meFolder = makeName(folders(ii).name,peDataFolder);
        allDFm = mean(allDFb(:,75:90,an),2); allDFm = repmat(allDFm,1,300);
        allDF = allDFb(:,:,an) - allDFm;
        trials2 = 1:10; trials1 = 11:20;
        badTrials = load(makeName('badTrials.txt',meFolder));
        g2trials = trials2(~ismember(trials2,badTrials));
        g1trials = trials1(~ismember(trials1,badTrials));
        
        GFP2(an,:) = mean(allDF(g2trials,:),1);
        GFP1(an,:) = mean(allDF(g1trials,:),1);
        GFPall(an,:) = mean(allDF([g1trials,g2trials],:),1);
    end
end
% Glut
mDF_Glu = mean(avgDF_Glu(:,70:135),2);
avgDF_Glu = avgDF_Glu - mDF_Glu;

t_GFP = (1:size(GFP1,2))/150; % time in seconds. 150 is frame rate
t_GFP = t_GFP - 0.6; % 0.6s is stim onset in GFP experiments
t_Glut = (1:size(avgDF_Glu,2))/150; % time in seconds. 150 is frame rate
t_Glut = t_Glut - 0.9; % 0.9s is stim onset in Glut experiments


function plotHiLo(ff,GFP2,GFP1,t_GFP)
%% GFP2 vs GFP1

shadedErrorBar(t_GFP,mean(GFP2,1),std(GFP2,[],1)/sqrt(size(GFP2,1)),{'color','r','linewidth',2},0.5);
shadedErrorBar(t_GFP,mean(GFP1,1),std(GFP1,[],1)/sqrt(size(GFP1,1)),{'color','b','linewidth',2},0.5); 

xlim([-0.1,1.4-0.05]); set(ff.ha,'XTick',[0,0.4,0.8,1.2])
ylim([-0.5 1.7]); set(ff.ha,'YTick',[-0.5,0,0.5,1])
set(ff.ha,'TickDir','out','FontSize',6);
set(ff.ha,'linewidth',1.25,'Fontweight','bold','FontSize',6);
xlabel('Time (s)');
hy = ylabel('Avg GFP Signal \DeltaF/F_0 %');
pos = get(hy,'Position'); pos(1) = pos(1) - 0.02; set(hy,'Position',pos);

x1 = 0; dx = 0.6/3; y1 = 1.45; dy = 0.15; fontsz = 7;
for nb = 1:2
    plot([x1-0.02 x1+3*dx],[y1,y1],'color','k','linewidth',1);
    plot([x1 x1],[y1,y1+dy],'color','b','linewidth',2);
    text(x1,y1+dy,'#1','Color','b','FontSize',fontsz,'VerticalAlignment','bottom','HorizontalAlignment','center');
    plot([x1 x1]+dx,[y1,y1+dy],'color','r','linewidth',2);
    text(x1+dx,y1+dy,'#2','Color','r','FontSize',fontsz,'VerticalAlignment','bottom','HorizontalAlignment','center');
    text(x1+dx/2,y1,'10s','Color','k','FontSize',fontsz,'VerticalAlignment','top','HorizontalAlignment','center');
    text(x1+2*dx,y1,'20s','Color','k','FontSize',fontsz,'VerticalAlignment','top','HorizontalAlignment','center');
    x1 = x1+3*dx;
end
text(x1+dx/2,y1+dy/2,' ....','Color','k','FontWeight','bold','FontSize',fontsz+3,'VerticalAlignment','middle','HorizontalAlignment','center');
text(x1+dx/2,y1+dy,'10 trials each','Color','k','FontSize',fontsz,'VerticalAlignment','baseline','HorizontalAlignment','center');

save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);


function plotGFPvsGlu(ff,avgDF_Glu,GFPall,t_GFP,t_Glut)
% plots
avgGlut = mean(avgDF_Glu);
semGlut = std(avgDF_Glu)/sqrt(size(avgDF_Glu,1));
shadedErrorBar(t_Glut,avgGlut,semGlut,{'color','k','linewidth',2},0.75);
avgGFP = mean(GFPall);
semGFP = std(GFPall)/sqrt(size(GFPall,1));
shadedErrorBar(t_GFP,avgGFP,semGFP,{'color','m','linewidth',2},0.75);
% ticks, limits, and labels
xlim([-0.1,1.4-0.05]); set(ff.ha,'XTick',[0,0.4,0.8,1.2])
ylim([-0.5 1.]); set(ff.ha,'YTick',[-0.5,0,0.5,1])
set(ff.ha,'TickDir','out','FontSize',6);
set(ff.ha,'linewidth',1.25,'Fontweight','bold','FontSize',6);
xlabel('Time (s)');
hy = ylabel('Avg Optical Signal \DeltaF/F_0 %');
pos = get(hy,'Position');pos(1) = pos(1) - 0.02; set(hy,'Position',pos);
% legend
thisCols = {'k','m'};
x1 = 1; x2 = x1+0.15; y1 = (0.75:0.2:7); y2 = y1;
legendFontSize = 6-1;
legs = {'Glut','GFP'};
for ii = [1,2]
    plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols{ii},'linewidth',2);
    text(x2+0.05,y1(ii),legs{ii},'Color',thisCols{ii},'FontSize',legendFontSize+2);
end
% saving as pdf
save2pdf(makeName('temp.pdf',getpdffolder),ff.hf,600);
