% Demo_model_prediction.m
%   
% This script demonstrates generation of model predictions for AWA neural
% response and behavior probability given input stimulation conditions.
% 20250531, DRA

%% Clear workspace and Load model
clc; clear all; close all;
if exist('20sOptimization.mat','file'); 
    load('20sOptimization.mat');
else 
    disp('20sOptimization.mat not found. Please run Optimization.m in the folder first.'); 
end

%% Stimulation Parameters 
stim.dt = 0.1;          % neural time step (seconds)
stim.dt_ds = 2;         % downsampled time step (for behavior) (seconds) 

stim.num_trials = 24;   % Number of repeated trials
stim.trial_dur = 60;    % Trial duration (seconds)
stim.t_init = 5;        % Stimulus onset (seconds)
stim.t_on = 50;         % Stimulus duration (seconds)

stim.conc = 0.6;        % relative concentration (log scale) 
                        % 1 = moderate stimulus; 1.25 = strong; 0.75 = weak

%% Predict Neural and Behavioral Responses

response = simulate_neural_activity_fn(params,stim);
behavior = simulate_behavioral_activity_fn(params,response);
beh_idx = summarize_behavior(behavior,stim);

%% Plot results

figure(1); clf;
tiledlayout(2,1);

stim_highlights = []; 
for i=1:stim.num_trials
    stim_highlights = [stim_highlights; (i-1)*stim.trial_dur + stim.t_init + [0 stim.t_on]];
end

% Plot neural response
nexttile;
plot(response.t_index / 60,response.calcium,'LineWidth',2);
ylabel({'AWA Neural','Response (dF/F0)'});
hilite(stim_highlights / 60,[],[0.9 0.9 1]);
xlabel('Experiment Time (min)');
xlim('tight');

title(sprintf('Stimulus on %ds every %ds, %d repeats, concentration: %0.2f', ...
    stim.t_on,stim.trial_dur,stim.num_trials,stim.conc));

% Plot behavior response
nexttile;
PlotSpacing = [1,0.5,1]; PlotCodes = 2:4; stair = 1;
behmat = [zeros(1,size(behavior.prob_ds,2));behavior.prob_ds];
stateplot3(behmat,NaN,(1:length(behmat))/60*stim.dt_ds,PlotSpacing,stair,PlotCodes);
set(findobj('FaceColor',[0 0 0]),'FaceColor',[0.5 0 0]);
ylabel('Behavior Probability');
xlabel('Experiment Time (min)');
ylim([0.9,3.5]);
hilite(stim_highlights / 60,[],[0.9 0.9 1]);

