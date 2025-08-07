% Optimization.m
%   
% This script generated optimized neural model parameters given an input 
% experimental dataset. It relies on the follwoing functions:
% - cost_fn.m: the cost function for optimization
% - simulate_neural_activity.m: the prediction of neural activity given the
%     current model parameters

%% Some Setup Steps
clc; clear all; close all;

%% Set params
stim.dt = 0.1;        % Sampling interval (s)
stim.t_init = 5;      % Stimulus onset (s)
stim.t_on = 20;       % stimulus duration (s)
stim.trial_dur = 60;  % time in (s) per trial
stim.num_trials = 24; % Number of trials to simulate
stim.conc = 1.0;     
stim.dt_ds = 2;       % 2 s interval for downsampling

% penalty parameters
params.q = 20;  % Q: inaccuracy of fast latent rep
params.r = 0.1; % R: rapid fluctuations
params.s = 2;   % S: energy expenditure

params.m_dim = 2;
params.n_dim = 20;
params.overlap = 5;

% parameters for z_c psychometric curve
params.c0 = 0.35;  
params.l1 = 0.1;
params.a1 = 8;
params.l2 = 0.35;
params.a2 = 2;

% parameters for calcium trace
params.baseline = 0.25; 
params.low_asymp = 0;
params.up_asymp = 1;
params.curve = 7.5;
params.lowpass = -0.145;
params.gain = 1;
params.intensity = 4.0;

params.SEED = 7;
params.noise_level = 0.01;

% Optimization parameters
params.lb = [-1, 0, 0, -0.1, 0, 0.25];
params.ub = [0, 1, 1.5, 0, 20, 0.5];
params.parameter_init = [-0.1, 0.30, 0.75, -0.075, 13.5, 0.3];
%params.lb = [-1, 0, 0, -0.1, 0, 0.25];
%params.ub = [0, 1, 1.5, 0, 20, 0.5];
%params.parameter_init = [-0.1, 0.30, 0.75, -0.075, 13.5, 0.3];
params.options = optimoptions(@fmincon,'OutputFcn',@outfun_vF,'Display', ...
    'iter','MaxIterations', 100);

% behavior model parameters
params.var.forward = 0.1;
params.var.pause = 0.6;
params.var.reverse = 0.3;
params.stim_on_mean = [1; 0];
params.stim_off_mean = [0; 0];
params.mat.transition = [0.7 0.15 0.15; ...
                         0.15 0.7 0.15; ...
                         0.15 0.15 0.7];
params.mat.valence = [0.1 0.8; ...
                      0.1 0.1; ...
                      0.8 0.1];
params.weights = [0.5, 0.5];

%% Load Training Data
load('ExpData_conc.mat');
% This dataset contains 25 trials with 20s (short) or 50s (long) pulses
% every 60s, at 4x10-6, 4x10-7, 4x10^-8, 4x10^-9 DA dilutions 
% (index 1,2,3,4 = 20s, 10^-6,-7,-8,-9);(index 5,6,7,8 = 50s, 10^-6,-7,-8,-9)
% Stimulus begions at 5s for 20s, i.e. from 5 - 25s or frames 50-250. 

%% Extract the 20s neural data
neural_act_20s = exp_neural.calcium(2,:);  % 20s pulse, 60s int, 10^-7 dil
t_min = exp_neural.t_s / 60; % time in min

% interpolate neural trace between trials (exponential)
dF = neural_act_20s; 
fr_per_trial = (stim.trial_dur)/stim.dt;
dF(1:fr_per_trial:end)=dF(2:fr_per_trial:end); % fix first frame error
exp_trials = floor(size(dF,2)/stim.trial_dur*stim.dt);
stim.num_trials = exp_trials;

exp_frames = 260:300; % frames per trial after odor-off (i.e. 26-30s)
dF_exp = dF;
for n = 0:(exp_trials-1)
    dFfit = dF(n*fr_per_trial + exp_frames)'; % calculate exp fit curve from last 4s
    tfit = (0:(length(dFfit)-1))';
    f=fit(tfit,dFfit,'exp1'); % fit exponential curve dropping to zero
    fill_fr = (max(exp_frames)+1):fr_per_trial;
    dF_exp(n*fr_per_trial + fill_fr) = f(fill_fr - min(exp_frames));
end
neural_act_20s_exp = dF_exp;

figure(1); clf; 

% plot data with interpolated points
nexttile;
    plot(t_min, neural_act_20s_exp);
    xlabel('Time (min)'); ylabel('dF/F0, exp.'); xlim([0, 24])

% plot just trials 1, 11, 21
nexttile;
    tps = [1:600,(1:600)+6000,(1:600)+12000];
    plot((1:length(tps))*stim.dt, neural_act_20s_exp(tps));
    xlabel('Time (s)'); ylabel('dF/F0, exp.'); title('Trials 1,10,20');

%% solve optimization
tic;
data_train_20s = neural_act_20s_exp(:,1:stim.num_trials*stim.trial_dur/stim.dt);
figure(2);
history_20s = runfmincon_vF(data_train_20s, params, stim);
toc;

%% Save Relevant Data from Optimization Run
params.model = history_20s.parameters(end,:);
save('20sOptimization.mat');
disp('Optimized parameters:');
disp(params.model);

%% Calculate and Plot fit
test_stim = stim;
%test_stim.t_on = 50;
response = simulate_neural_activity_fn(params,test_stim);
err = norm(response.calcium - data_train_20s);

figure(3); clf;
plot(response.t_index,[response.calcium;response.gamma;response.nu;data_train_20s]);
xlabel('Time (s)'); ylabel('dF/F0 or Latent Variable');

behavior = simulate_behavioral_activity_fn(params, response);
figure(4); clf;
plot(response.t_index,behavior.forward);
xlabel('Time (s)'); ylabel('Forward behavior probability');

%% Plot optimization evolution
figure(5); clf; hold on; 
iter_list = [1 2 4 size(history_20s.parameters,1)];
for jj = 1:length(iter_list)
    params.model = history_20s.parameters(iter_list(jj),:);
    response = simulate_neural_activity_fn(params,stim);
    plot(response.t_index,response.calcium);
end
%legend(sprintf('Iteration %d',iter_list))

%% Functions
function [history] = runfmincon_vF(exp_data,params,stim)
   
    history.parameters = [];
    history.fval = [];
    
    % Set up options 
    options = optimoptions(@fmincon,'OutputFcn',@outfun_vF, 'Display', 'iter',...
                                 'MaxIterations', 100); %with information about each iteration.
    
    % Define cost function to minimize
    fun = @(parameters) cost_fn(parameters,exp_data,params,stim);
    
    % Constraints and Initial conditions
    lb = params.lb;
    ub = params.ub;
    parameters0 = params.parameter_init;
    
    % Iterate parameters
    parameters = fmincon(fun,parameters0,[],[],[],[],lb,ub,[],options);
    
    % Display output during interation
    function stop = outfun_vF(parameters, optimValues, state)
        stop = false;
        switch state
            case 'init'
                hold on;
            case 'iter'
                %Concatenate current point and objective function value with history.
                history.fval = [history.fval; optimValues.fval];
                history.parameters = [history.parameters; parameters];
                plot(optimValues.iteration, optimValues.fval,'o--','Color','k');
            case 'done'
                %plot(0:optimValues.iteration, history.fval,'k');
                hold off
            otherwise
        end
    end

end