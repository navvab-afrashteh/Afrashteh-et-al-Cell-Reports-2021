function motifs = getMotifsFromCCResultsSam_simple(animalNumber,listR,ii,thn,match,CrossOrMax,LMHname,varargin)
multLevels = 1;
if nargin > 7
   multLevels = varargin{1}; 
end

mainDataFolder = getMainDataFolder;
dataFolder = makeName(animalNumber{1},mainDataFolder);
psDataFolder = makeName('pSpon',dataFolder);
peDataFolder = makeName('pEvoked',dataFolder);
fileName = makeName('evokCorrSpon.mat',psDataFolder);
if ~exist(fileName,'file') == 2
    motifs = [];
    return;
end
load(fileName);

fileName = makeName('evokCorrevok.mat',peDataFolder);
load(fileName);
typeCC = 'mean';
cmdText = sprintf('selectedEvokTh = evokCorrevok.%sCC(end);',typeCC);
% selectedEvokTh = evokCorrevok.maxCC(end);
eval(cmdText);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
testString = sprintf('evokCorrevok.%sCCFLname',typeCC);
for jj = 1:length(evokCorrSpon)
    if strcmp(eval(testString),evokCorrSpon{jj}.eFolder)
        break;
    end
end
allCCsRowFL = jj;

testString = sprintf('evokCorrevok.%sCCHLname',typeCC);
for jj = 1:length(evokCorrSpon)
    if strcmp(eval(testString),evokCorrSpon{jj}.eFolder)
        break;
    end
end
allCCsRowHL = jj;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for jj = 1:length(evokCorrSpon)
    if strcmp(LMHname,evokCorrSpon{jj}.eFolder)
        break;
    end
end
allCCsRowMatch = jj;

[lists, names]= getEvokedListsDiffStimAmps(animalNumber(1),multLevels);
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
if selectedEvokTh < thmin
    selectedEvokTh = thmin; 
end

folderName = makeName(listR{ii},psDataFolder);
fileName = makeName('allCCs.mat',folderName);
if exist(fileName,'file') == 2
    load(fileName);
    for jj = allCCsRowMatch 
        oneCC = allCCs(jj,:);
        frames = cleanCCFramesList(oneCC,selectedEvokTh*thn,CrossOrMax);
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
