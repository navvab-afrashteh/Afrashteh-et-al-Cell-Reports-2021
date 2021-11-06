function manualEvokedTrialsVerify_YFP ()
close all
clc
mainDataFolder = getMainDataFolderYFP; 
animalNumber =  {'190533'; '190545'};
stim = {'BC' 'HL'};
quitthis = 0;
Glut = 0; YFP = 1;
Fpass1 = 0.5; Fpass2 = 6;
Fpass1 = round(Fpass1*100);
Fpass2 = round(Fpass2);
for an = 1:length(animalNumber)
    if quitthis
        break;
    end
    dataFolder = makeName(animalNumber{an},mainDataFolder);
    mask.bigMask = getMask(animalNumber(an),'ec',1.1,'glut',Glut,'yfp',YFP);
    mask.mask = getMask(animalNumber(an),'glut',Glut,'yfp',YFP);
    peDataFolder = makeName('pEvoked',dataFolder);
    ignoreFolders = [];
    ofolders = dir(peDataFolder);
    folders = cleanFolderList(ofolders,ignoreFolders);
    for ii = 1:length(folders)
        [an ii]
        meFolder = makeName(folders(ii).name,peDataFolder);
        badTrials = [];
        osFolders = dir(meFolder);
        sFolders = cleanFolderList(osFolders,badTrials);
        for jj = length(sFolders):-1:1
            trialFolder = makeName(sFolders(jj).name,meFolder);
            fileName = sprintf('%s\\ImgSeq_BPF_0%dHz_%dHz.mat',trialFolder, Fpass1, Fpass2);
            mFile = matfile(fileName);
            ImgSeq = mFile.ImgSeq;
            maskdf = getMaskedValues(ImgSeq,mask.mask);
            allDF(jj,:) = mean(maskdf,2);
            jj
        end
        jj=1;
        badtrials = false(1,length(sFolders));
        while jj <= length(sFolders)
            maskdf = allDF(jj,:);
            frame = 91; % stim onset
            figure(2);clf;
            plot(maskdf);hold on;
            plot([frame frame],[min(maskdf) max(maskdf)],'r');
            title(['trial #',num2str(jj)],'interpreter','none');
            
             while 1
                ch = getkey;
                if ch == 113 % q
                    quitthis = 1;
                    break;
                end
                if ch == 29 % -->
                    jj = jj + 1;
                    break;
                end
                if ch == 28 % <--
                    jj = jj - 1;
                    break;
                end
                
                if ch == 115 % s
                    badtrials(jj) = true;
                end
            end
        end
        
        avgdf = allDF(~badtrials,:);
        avgdf = nanmean(avgdf);
        figure(2);clf;
        plot(avgdf);hold on;
        plot([frame frame],[min(avgdf) max(avgdf)],'r');
        title(['avg of ',num2str(sum(~badtrials))],'interpreter','none');
        clc
        disp(find(badtrials))
    end
end
n=0;

gCols = 3;
gRows = 5;
totalPlots = gRows * gCols;
grcs = reshape(1:totalPlots,gCols,gRows)';
figure(1);clf;
for an = 1:5
    for ii = 1:3
%         titleText = sprintf('
        subplot(gRows,gCols,grcs(an,ii));
        thisdf = allACdf{an,ii};
        plot(thisdf);hold on;plot([136 136],[-2.5 1.5],'r');plot(mean(thisdf,2),'k','linewidth',2);
        ylim([-1 2]);
    end
end
