function motifs = getMotifsFromCCResultsSam_simple_glut(animalNumber,listR,ii,thn,match,CrossOrMax,LMHname,fde,varargin)
multLevels = 1;
if nargin > 8
   multLevels = varargin{1}; 
end

fr10 = 0;

mainDataFolder = getMainDataFolderGlut;
dataFolder = makeName(animalNumber{1},mainDataFolder);
psDataFolder = makeName('pSpon',dataFolder);
peDataFolder = makeName('pEvoked',dataFolder);
if fr10
    fileName = makeName('evokCorrSpon_10fr.mat',psDataFolder);
else
    if ~isempty(fde.suffix)
        fileName = makeName(sprintf('evokCorrSpon_%s.mat',fde.suffix),psDataFolder);
    else
        fileName = makeName('evokCorrSpon.mat',psDataFolder);
    end
end
if ~exist(fileName,'file') == 2
    motifs = [];
    return;
end
load(fileName);

% fileName = makeName('evokCorrevok.mat',peDataFolder);
% load(fileName);
% typeCC = 'mean';
% cmdText = sprintf('selectedEvokTh = evokCorrevok.%sCC(end);',typeCC);
% % selectedEvokTh = evokCorrevok.maxCC(end);
% eval(cmdText);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testString = sprintf('evokCorrevok.%sCCFLname',typeCC);
% for jj = 1:length(evokCorrSpon)
%     if strcmp(eval(testString),evokCorrSpon{jj}.eFolder)
%         break;
%     end
% end
% allCCsRowFL = jj;
% 
% testString = sprintf('evokCorrevok.%sCCHLname',typeCC);
% for jj = 1:length(evokCorrSpon)
%     if strcmp(eval(testString),evokCorrSpon{jj}.eFolder)
%         break;
%     end
% end
% allCCsRowHL = jj;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for jj = 1:length(evokCorrSpon)
    if strcmp(LMHname,evokCorrSpon{jj}.eFolder)
        break;
    end
end
allCCsRowMatch = jj;

[lists, names]= getEvokedListsDiffStimAmps_glut(animalNumber(1),multLevels);
if strcmp(match,'FL')
    %     allCCsRowMatch = allCCsRowFL;
    allCCsRowFL = allCCsRowMatch;
    for jj = 1:length(names{1})
        if strcmp(names{1}{jj},evokCorrSpon{allCCsRowFL}.eFolder)
            break;
        end
    end
    stimAmpNumber = jj;
elseif strcmp(match,'HL')
    %     allCCsRowMatch = allCCsRowHL;
    allCCsRowHL = allCCsRowMatch;
    for jj = 1:length(names{2})
        if strcmp(names{2}{jj},evokCorrSpon{allCCsRowHL}.eFolder)
            break;
        end
    end
    stimAmpNumber = jj;
elseif strcmp(match,'AC')
    %     allCCsRowMatch = allCCsRowHL;
    allCCsRowHL = allCCsRowMatch;
    for jj = 1:length(names{4})
        if strcmp(names{4}{jj},evokCorrSpon{allCCsRowHL}.eFolder)
            break;
        end
    end
    stimAmpNumber = jj;
elseif strcmp(match,'VC')
    %     allCCsRowMatch = allCCsRowHL;
    allCCsRowVC = allCCsRowMatch;
    for jj = 1:length(names{3})
        if strcmp(names{3}{jj},evokCorrSpon{allCCsRowVC}.eFolder)
            break;
        end
    end
    stimAmpNumber = jj;
end

thmin = 0.4;
% if selectedEvokTh < thmin
selectedEvokTh = thmin;
% end
temp = strfind(listR,fde.suffix);
ind = find(not(cellfun('isempty', temp)));
if isempty(ind)
    folderName = makeName(listR{3},psDataFolder);
    if fr10
        fileName = makeName('allCCs_10fr.mat',folderName);
    else
        fileName = makeName('allCCs.mat',folderName);
    end
else
    folderName = makeName(listR{ind},psDataFolder);
    fileName = makeName(sprintf('allCCs_%s.mat',fde.suffix),folderName);
end
if exist(fileName,'file') == 2
    load(fileName);
    for jj = allCCsRowMatch 
        oneCC = allCCs(jj,:);
        frames = cleanCCFramesList_glut(oneCC,selectedEvokTh*thn,CrossOrMax);
%         frames = cleanCCFramesList(oneCC,thn,CrossOrMax);
    end
    motifs.frames = frames;
    motifs.stimAmpNum = stimAmpNumber;
    motifs.allCCs = oneCC;
else
    motifs.frames = [];
end

if isempty(motifs.frames)
    motifs = [];
end
