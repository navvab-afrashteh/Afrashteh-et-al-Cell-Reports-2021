% First load the data
loadData_old

%% Figure 1
% plot time series, time to peak hist, and DF_F0 hist for VSDI FL and HL.
% Also montages
ss = 1; % FL
Artur_fig_bar_graph_amplitude_riseTime_1(ss);
Artur_figure_Montage_Amplitude_1 (ss)
ss = 2; % HL
Artur_fig_bar_graph_amplitude_riseTime_1(ss);
Artur_figure_Montage_Amplitude_1 (ss)

% plot DF_F0 distribution (PDF) for VSDI signal FL and HL.
ss = 1; % FL
Artur_figure_1_Distributions_2(ss);
ss = 2; % HL
Artur_figure_1_Distributions_2(ss);

%% Figure 2
% plot time series, time to peak hist, and DF_F0 hist for Glut AC.
Artur_fig_bar_graph_amplitude_riseTime_1_glut;
% Here lo and Med are combined
Artur_fig_bar_graph_amplitude_riseTime_1_glut_LoMedCombined;
% plot DF_F0 distribution (PDF) for Glut AC.
Artur_figure_1_Distributions_2_glut;
Artur_figure_1_Distributions_2_glut_LoMedCombined;
% montage
Artur_figure_Montage_Amplitude_1_glut
%% Figure 3 and S3
% plot time series and Speed hist for VSDI FL and HL.
% plot Speed distribution (PDF) for VSDI FL and HL.
for ss = 1:2 % ss=1 for FL stim and ss=2 for HL stim
    roi = 2; % FL roi
    Artur_fig_bar_graph_speed_directions_2(ss,roi);
    Artur_figure_Speed_Distributions_2 (ss,roi);
    roi = 3; % HL roi
    Artur_fig_bar_graph_speed_directions_2(ss,roi);
    Artur_figure_Speed_Distributions_2 (ss,roi);
    roi = 4; % PtA roi
    Artur_fig_bar_graph_speed_directions_2(ss,roi);
    Artur_figure_Speed_Distributions_2 (ss,roi);
    % speed montages
    Artur_figure_Montage_Speed_1(ss)
end
%% Figure 4
% plot time series and Speed hist for Glut AC.
Artur_fig_bar_graph_speed_directions_2_glut_average;
Artur_fig_bar_graph_speed_directions_2_glut_avg_LoMedCombined(1,2);
% plot Speed distribution (PDF) for Glut AC.
Artur_figure_Speed_Distributions_2_glut_average;
Artur_figure_Speed_Distributions_2_glut_average_LoMedCombined;
% Speed montage
Artur_figure_Montage_Speed_1_glut
%% Figure 5 - trajectory and permutation
for ss = 1:2 % ss=1 for FL stim and ss=2 for HL stim
    for roi = 2:4  % FL, HL, PtA ROIs
        Artur_fig_rose_graph_directions_1(ss,roi)
    end
end
% generate trajectories and example trajectories and fractal dimension bar graph
Artur_Trajectory_VSDI;
% permutation and distributions for trajectories
Artur_Trajectory_Permutation_VSDI;

%% Figure 6 - Rose plots, trajectory and permutation 
ss = 1;
for roi = 2:4
    Artur_fig_rose_graph_directions_1_glut(ss,roi);
end
% generate trajectories and example trajectories and fractal dimension bar graph
Artur_Trajectory_glut;
% permutation and distributions for trajectories
Artur_Trajectory_Permutation_glut;

%% New Figure - Spontaneous Measurments Relations
% speed, flow directionality, and trajectory fractal dimension relation to
% signal amplitude
selE = [1,2,3]; % selected evoked levels
for ss = 1:2 % for FL and HL 
    for roi = 1:4
        SponMeasRelation(ss,roi,selE);
    end
end
for roi = [1,5]
    SponMeasRelationVC(roi); % for VC
end
for roi = 1:4
    SponMeasRelationGlut(roi,selE); % for AC glut
end

%% S1 - trajectory illustration
traj_illustration;
% template matching illustration
TempMatchIllustration
%% Supplementary for VC
% plot Speed distribution (PDF) for VSDI VC.
Artur_figure_Speed_Distributions_2_VC
% vsd time series, FD bar graph, time to peak, VSD Amp bar graph
Artur_fig_bar_graph_amplitude_riseTime_1_VC
% Speed time series, Speed bar graph, Direction A.U.
Artur_fig_bar_graph_speed_directions_2_VC
% Rose plot for direction
Artur_fig_rose_graph_directions_1_VC
% generate trajectories and example trajectories
Artur_Trajectory_VSDI_VC