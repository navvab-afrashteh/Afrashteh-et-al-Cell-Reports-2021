clear;
clc;
mainDataFolder = getMainDataFolderYFP;

% All animals 
animalNumber =  {'190533'; '190545'};
Rotate = 0;
theta = 0;
filtOptions.Fstop1 = 0.25; filtOptions.Fpass1 = 0.5; 
filtOptions.Fpass2 = 6; filtOptions.Fstop2 = 6.5;
filtOptions.dur = 5;
%% convert Masks to .mat files
% convertMask_TiffMat(animalNumber,'yfp');

%% Preprocessing spontaneous activity and run optical flow on the results
% runBandPassFilter_Spon(animalNumber,Rotate,theta,'yfp',filtOptions);
% R1_Spon_SpatioTemp_YFP;

%% Preprocessing evoked activity and run optical flow on the results
% runBandPassFilter_Evoked(animalNumber,Rotate,theta,'yfp',filtOptions);
% manualEvokedTrialsVerify_YFP;
% findTrialAverageAll_BPF (animalNumber,'yfp',filtOptions)

%% save template and ROIs
saveEvokedTemplate_GFP_BPF(animalNumber,'gfp',filtOptions);
saveROIs_GFP_BPF;






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% findTrialAverageAll_750Frames_BPF (animalNumber,'glut',filtOptions)
% glut_down_sample_evoked_BPF(animalNumber,'glut',filtOptions)
% saveEvokedTemplate_diffAmps_glut_BPF(animalNumber,'glut',filtOptions)

% filtOptions.Fstop1 = 0.75; filtOptions.Fpass1 = 1; 
% filtOptions.Fpass2 = 6; filtOptions.Fstop2 = 6.5;
% filtOptions.dur = 5;
% runBandPassFilter_Evoked(animalNumber,Rotate,theta,'glut',filtOptions);
% findTrialAverageAll_BPF (animalNumber,'glut',filtOptions)
% findTrialAverageAll_750Frames_BPF (animalNumber,'glut',filtOptions)
% glut_down_sample_evoked_BPF(animalNumber,'glut',filtOptions)
% 
% runOF_evok_CLG_BPF(animalNumber,'glut');

%% calculate cross-correlation between evoked and spontaneous activity
% filtOptions.Fstop1 = 0.25; filtOptions.Fpass1 = 0.5; 
% filtOptions.Fpass2 = 6; filtOptions.Fstop2 = 6.5;
% filtOptions.dur = 5;
% evokCorrCoeffSpon_new_glut_BPF(animalNumber,'glut',filtOptions);

% filtOptions.Fstop1 = 0.75; filtOptions.Fpass1 = 1; 
% filtOptions.Fpass2 = 6; filtOptions.Fstop2 = 6.5;
% filtOptions.dur = 5;
% evokCorrCoeffSpon_new_glut_BPF(animalNumber,'glut',filtOptions);

%% save ROIs
% saveROIs_glut_BPF % not completed - did not run. Used Previously saved ROIs.

%% save all spont
% filtOptions.Fstop1 = 0.25; filtOptions.Fpass1 = 0.5; 
% filtOptions.Fpass2 = 6; filtOptions.Fstop2 = 6.5;
% filtOptions.dur = 5;
% saveSponAll_glut_BPF(animalNumber,'glut',filtOptions)
% save_sData_glut_BPF(animalNumber,'glut',filtOptions)
% filtOptions.Fstop1 = 0.75; filtOptions.Fpass1 = 1; 
% filtOptions.Fpass2 = 6; filtOptions.Fstop2 = 6.5;
% filtOptions.dur = 5;
% saveSponAll_glut_BPF(animalNumber,'glut',filtOptions)
% save_sData_glut_BPF(animalNumber,'glut',filtOptions)

%% save all evoked
% filtOptions.Fstop1 = 0.25; filtOptions.Fpass1 = 0.5; 
% filtOptions.Fpass2 = 6; filtOptions.Fstop2 = 6.5;
% filtOptions.dur = 5;
% save_eData_glut_BPF(animalNumber,'glut',filtOptions);

% filtOptions.Fstop1 = 0.75; filtOptions.Fpass1 = 1; 
% filtOptions.Fpass2 = 6; filtOptions.Fstop2 = 6.5;
% filtOptions.dur = 5;
% save_eData_glut_BPF(animalNumber,'glut',filtOptions);

