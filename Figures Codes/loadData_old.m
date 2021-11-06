clear all
clc
addpath('Y:\homes\navvab.afrashteh\CloudStation\Navvab-Sam\NeuroPhotonic Paper\Test')
add_to_path;

mainDataFolder = getMainDataFolder;

%% FL and HL

fileName = makeName('eData_artur.mat',mainDataFolder);
eData = load(fileName);
mdataE = eData.mdata;

fileName = makeName('sData.mat',mainDataFolder);
sData = load(fileName);
mdataS = sData.mdata;

fileName = makeName('sDataDFth.mat',mainDataFolder);
sDataDFth = load(fileName);
mdataSDFth = sDataDFth.mdata;

% for an = 1:length(mdataE{1}.animalNumber)
%     [old_lists{an} old_names{an}]= getEvokedListsDiffStimAmps(mdataE{1}.animalNumber(an));
% end

selLMH.selLMH_all_majid {1} = [1 4 6;1 3 5;2 3 5;1 2 4;1 2 4];
selLMH.selLMH_all_majid {2} = [1 2 5;2 3 5;1 2 3;1 3 4;1 2 3];
selLMH.selLMH_all {1} = [3 4 6;2 3 4;2 3 5;2 3 4;2 3 4];
selLMH.selLMH_all {2} = [1 3 5;2 3 5;1 2 3;2 3 4;1 2 3];

anGroup{1} = {[1],[2],[4]}; anGroup{2} = {[1],[2],[3]};
anGroup{3} = {[1],[2],[4]}; anGroup{4} = {[1],[2],[3]};
anGroup{5} = {[1],[2],[3]};
mdataE{1}.anGroup = anGroup;
mdataE{1}.selLMH = selLMH;

% for HL
anGroup{1} = {[1],[3],[5]}; anGroup{2} = {[1],[2],[4]};
anGroup{3} = {[1],[2],[3]}; anGroup{4} = {[1],[2],[3]};
anGroup{5} = {[1],[2],[3]};
mdataE{2}.anGroup = anGroup;
mdataE{2}.selLMH = selLMH;

%% GLUT simple

fileName = makeName('eData_glut_artur.mat',mainDataFolder);
eData = load(fileName);eData.mdata{1}.suffix = '';
mdataE_glut = eData.mdata;
mdataE_glut{1}.selLMH.selLMH_all {1} = [1 2 3;1 2 3;1 2 3;1 2 3;1 2 3];

fileName = makeName('sData_glut.mat',mainDataFolder);
sData = load(fileName);
mdataS_glut = sData.mdata;


%% GLUT HPF

fileName = makeName('eData_glut_artur_HPF.mat',mainDataFolder);
eData = load(fileName); eData.mdata{1}.suffix = 'HPF';
mdataE_glut_HPF = eData.mdata;
mdataE_glut_HPF{1}.selLMH.selLMH_all {1} = [1 2 3;1 2 3;1 2 3;1 2 3;1 2 3];
% mdataE_glut_HPF{1} = rmfield(mdataE_glut_HPF{1},'davg');
% mdataE_glut_HPF{1} = rmfield(mdataE_glut_HPF{1},'d1avg');
% saveAvg_eData_glut_HPF_HS(mdataE_glut_HPF);

fileName = makeName('sData_glut_HPF.mat',mainDataFolder);
sData = load(fileName);
mdataS_glut_HPF = sData.mdata;


%% GLUT HPF HS
fileName = makeName('eData_glut_artur_HPF_HS.mat',mainDataFolder);
eData = load(fileName);
eData.mdata{1}.suffix = 'HPF_HS';
mdataE_glut_HPF_HS = eData.mdata;
mdataE_glut_HPF_HS{1}.selLMH.selLMH_all {1} = [1 2 3;1 2 3;1 2 3;1 2 3;1 2 3];

fileName = makeName('sData_glut_HPF_HS.mat',mainDataFolder);
sData = load(fileName);
mdataS_glut_HPF_HS = sData.mdata;


%% GLUT BPF 050Hz_6Hz

fileName = makeName('eData_glut_BPF_050Hz_6Hz.mat',mainDataFolder);
eData = load(fileName); eData.mdata{1}.suffix = 'BPF_050Hz_6Hz';
mdataE_glut_BPF_050Hz_6Hz = eData.mdata;
mdataE_glut_BPF_050Hz_6Hz{1}.selLMH.selLMH_all {1} = [1 2 3;1 2 3;1 2 3;1 2 3;1 2 3];

fileName = makeName('sData_glut_BPF_050Hz_6Hz.mat',mainDataFolder);
sData = load(fileName);
mdataS_glut_BPF_050Hz_6Hz = sData.mdata;


%% GLUT BPF 0100Hz_6Hz

fileName = makeName('eData_glut_BPF_0100Hz_6Hz.mat',mainDataFolder);
eData = load(fileName); eData.mdata{1}.suffix = 'BPF_0100Hz_6Hz';
mdataE_glut_BPF_0100Hz_6Hz = eData.mdata;
mdataE_glut_BPF_0100Hz_6Hz{1}.selLMH.selLMH_all {1} = [1 2 3;1 2 3;1 2 3;1 2 3;1 2 3];

fileName = makeName('sData_glut_BPF_0100Hz_6Hz.mat',mainDataFolder);
sData = load(fileName);
mdataS_glut_BPF_0100Hz_6Hz = sData.mdata;

%% VC
fileName = makeName('eData_artur_VC.mat',mainDataFolder);
eData = load(fileName);
mdataE_VC = eData.mdata;
mdataE_VC{3}.selLMH.selLMH_all {1} = [1 2 3;1 2 3;1 2 3;1 2 3;1 2 3];

fileName = makeName('sDataDFThVC.mat',mainDataFolder);
sData = load(fileName);
mdataS_VC = sData.mdata;

%% AC

% AC anaesthesia
mainDataFolder = getMainDataFolderAC;

fileName = makeName('sData_AC.mat',mainDataFolder);
sData = load(fileName);
mdataS_AC = sData.mdata;

fileName = makeName('eData_AC.mat',mainDataFolder);
eData = load(fileName);
mdataE_AC = eData.mdata;
mdataE_AC{1}.selLMH.selLMH_all = mdataE_AC{1}.anGroup;

%% Meta Data Related to making figures
meta.stimNames = {'Spon','Lo-Ev','Med-Ev','Hi-Ev'};
meta.stimNamesVC = {'Spon','Ev'};
meta.stimNames1 = {'Sp','L-E','M-E','H-E'};
meta.stimNamesComplete = {'Spontaneous','Lo-Evoked','Med-Evoked','Hi-Evoked'};
meta.colors = getColors;
% figure related settings
meta.annotationFontSize = 8; meta.axesFontSize = 6; meta.legendFontSize = 5; meta.xlabelFontSize = 6; meta.ylabelFontSize = 6;
meta.titleFontSize = 8;
meta.amplitude_threshold = 0;
try
    Cmap = load('T:\CloudStation\common\selfmade3.txt');
catch
    Cmap = load('selfmade3.txt');
end
Cmap = Cmap(:,2:4)/255;
meta.Cmap = Cmap;

