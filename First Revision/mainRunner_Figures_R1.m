% First load the data
run('Y:\homes\navvab.afrashteh\CloudStation\Navvab-Sam\NeuroPhotonic Paper\Test\Figures_new_code_AC\loadData_old.m')

%% Figure 5 R1
% generate trajectories and example trajectories and fractal dimension bar graph
R1_Artur_Trajectory_VSDI;
% permutation and distributions for trajectories
R1_Artur_Trajectory_Permutation_VSDI;

%% Figure 6 R1
R1_Artur_Trajectory_glut
R1_Artur_Trajectory_Permutation_glut

%% New figures R1
% Distribution of time duration os spon and evok and calc spon motifs duration ratio to all spon recordings
ss = 1; R1_EVokSponDuration_VSD(ss);
ss = 2; R1_EVokSponDuration_VSD(ss);
ss = 3; R1_EVokSponDuration_VSD(ss);
ss = 1; R1_EVokSponDuration_Glut(ss);
R1_SponPCCEvent_TimePercentage;

% PSD for spon recordings. The signal is avgeraged in the entire imaging window.
ss=1; R1_SponPSD_VSD(ss); % for FL and HL animals
ss=3; R1_SponPSD_VSD(ss); % for VC animals
ss=1; R1_SponPSD_Glu(ss); % for AC Glut animals

% spontaneous events rate bar graph for all modalities. Also, PDF of
% inter-event-interval times for all modalities.
R1_SponEvent_timing;

%% Figure 1 R1
% plot time series, time to peak hist, and DF_F0 hist for VSDI FL and HL.
ss = 1; % FL
R1_Artur_fig_bar_graph_amplitude_riseTime_1(ss);
ss = 2; % HL
R1_Artur_fig_bar_graph_amplitude_riseTime_1(ss);

%% Figure 2 R1
% plot time series, time to peak hist, and DF_F0 hist for Glut AC.
R1_Artur_fig_bar_graph_amplitude_riseTime_1_glut_LoMedCombined

%% Figure S4 R1
R1_Artur_fig_bar_graph_amplitude_riseTime_1_VC

%% Figure S6 R1
R1_Artur_Trajectory_VSDI_VC
R1_Artur_Trajectory_Permutation_VSDI_VC

%% Figure 7 R1 - Spontaneous Measurments Relations
% speed, flow directionality, and trajectory fractal dimension relation to
% signal amplitude
selE = [1,2,3]; % selected evoked levels
for ss = 1:2 % for FL and HL 
    for roi = 1:4
        R1_SponMeasRelation(ss,roi,selE);
    end
end
for roi = [1,5]
    R1_SponMeasRelationVC(roi); % for VC
end
for roi = 1:4
    R1_SponMeasRelationGlut(roi,selE); % for AC glut
end

%% R1: YFP figures
% spatio-temporal characteristics of spon activity
R1_Spon_SpatioTemp_YFP



%% Figure 6 - Rose plots, trajectory and permutation 
ss = 1;
for roi = 2:4
    Artur_fig_rose_graph_directions_1_glut(ss,roi);
end
% generate trajectories and example trajectories and fractal dimension bar graph
Artur_Trajectory_glut;
% permutation and distributions for trajectories
Artur_Trajectory_Permutation_glut;
