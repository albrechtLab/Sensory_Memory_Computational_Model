function response = simulate_neural_activity_fn (params,stim,conc,pulse_dur)
%
% function response = simulate_neural_activity_fn (params,stim)
%
%   OUTPUT/INPUTS:
%
%	response.  	calcium, x, nu, gamma, t_index  
%
%	stim.		dt, t_init, t_on, trial_dur, num_trials, conc   {stimulation parameters}
%
%	params.		model [drift_nu, drift_nu_gamma, tau_nu, drift_gamma, tau_gamma, beta]  {optimization}
%			    q, r, s,					{penalty values}
%			    m_dim, n_dim, overlap,		{dimensions}
%			    c0, l1, a1, l2, a2 			{z_c}
%			    baseline, low_asymp, up_asymp, curve, lowpass, gain, intensity		{calcium transform}
%			    SEED, noise_level           {computational constants}

% update parameters
if exist('conc','var') stim.conc = conc; end
if exist('pulse_dur','var') stim.t_on = pulse_dur; end

rng(params.SEED); 

%% Set parameters along the time axis
t_end = stim.trial_dur * stim.num_trials;
t_index = 0:stim.dt:t_end - stim.dt;
response.t_index = t_index;

%% Set dimensions of the model
m_dim = params.m_dim;
n_dim = params.n_dim;
overlap = params.overlap;
b_matrix = weighting_matrix(m_dim, n_dim, overlap);

%% Set parameters of the model

drift_nu = params.model(1); 	
drift_nu_gamma = params.model(2); 	
tau_nu = params.model(3); 		
drift_gamma = params.model(4); 	
tau_gamma = params.model(5); 	
beta = params.model(6); 		


%% Set parameters for target tracking
zc_fun = @(c)  (params.l1 + (1-params.l1).*exp(-params.a1.*(c - params.c0).^2)).*(c<=params.c0) + ...
    (params.l2 + (1-params.l2).*exp(-params.a2.*(c - params.c0).^2)).*(c>params.c0);

z_c = zc_fun(stim.conc);
z_target = zeros(params.m_dim, stim.trial_dur/stim.dt);
z_target(1, stim.t_init/stim.dt + 1: floor((stim.t_init+stim.t_on)/stim.dt)) = z_c;
z = repmat(z_target, 1, stim.num_trials);

%% Set parameters for optimization
null_mn = zeros(m_dim, n_dim);
null_nm = zeros(n_dim, m_dim);

A_modified = [(drift_gamma/tau_gamma)*eye(m_dim), 	(1/tau_gamma)*beta*eye(m_dim), 	null_mn;...
    (drift_nu_gamma/tau_nu)*eye(m_dim), 	(drift_nu/tau_nu)*eye(m_dim), 	(1/tau_nu)*b_matrix;...
    null_nm,               		null_nm,          	zeros(n_dim)];
aux_eig_A = -1e-5 * ones(1, m_dim);
aux_A_matrix = diag(aux_eig_A);
A_final = blkdiag(A_modified, aux_A_matrix);
B_final = [zeros(2*m_dim, n_dim); eye(n_dim); zeros(m_dim, n_dim)];

%% Penalty matrices

Q = params.q * eye(m_dim);
R = params.r * eye(n_dim);
S = params.s * eye(n_dim);

Q_modified = [  	   Q	null_mn   -Q	; ...
                    null_nm    S 	null_nm	; ...
                      -Q 	null_mn	   Q	];
Q_final = blkdiag(1e-5*eye(m_dim), Q_modified);
R_final = R;

%% Solve Riccati Equation

K0 = zeros(size(A_final));
tspan = t_index;
[T,K] = ode45(@(t,K) mRiccati_F(t,K,A_final,B_final,R_final,Q_final), tspan, K0);

%% Set parameters for scaling response from a.u. to a positive range
stim_A_num = ceil((overlap + n_dim)/ 2);
pure_stim_A = stim_A_num - overlap;

%% Forward simulate the model

% Initialize response variables
response.x = zeros(n_dim, length(t_index));
response.nu = zeros(m_dim, length(t_index));
response.gamma = zeros(m_dim, length(t_index));
x_calcium_raw = ones(n_dim, length(t_index));
x_scaled = zeros(n_dim, length(t_index));
x_scaled(:,1) = params.low_asymp + (params.up_asymp - params.low_asymp)*sigmoid(x_scaled(:,1), params.baseline, params.curve);

%% Euler's approx of diff eq

for i = 1:length(t_index)-1
    Kt = reshape(K(i+1,:),size(A_final));
    W = inv(R_final)*B_final'*Kt; % feedback matrix
    W_mat{i} = -W;
    
    W_1 = -W(:,1:m_dim);
    W_2 = -W(:,m_dim+1:2*m_dim);
    W_3 = -W(:,2*m_dim+1:2*m_dim+n_dim);
    W_4 = -W(:,n_dim+2*m_dim+1:end);
    
    response.gamma(:,i+1) = response.gamma(:,i)+stim.dt*((drift_gamma/tau_gamma)*response.gamma(:,i) + (1/tau_gamma)*beta*response.nu(:,i)+ params.noise_level*randn);
    response.nu(:,i+1) = response.nu(:,i)+stim.dt*((drift_nu/tau_nu)*response.nu(:,i)+ (drift_nu_gamma/tau_nu)*response.gamma(:,i)+(1/tau_nu)*b_matrix*response.x(:,i)+ params.noise_level*randn);
    response.x(:,i+1) = response.x(:,i)+stim.dt*(W_1*response.gamma(:,i)+W_2*response.nu(:,i)+W_3*response.x(:,i)+ W_4*z(:,i)+ params.noise_level*randn);
    
    %   Scale to firing rate
    x_scaled(:,i+1) = params.low_asymp + (params.up_asymp - params.low_asymp)*sigmoid(response.x(:,i+1), params.baseline, params.curve);
    
end
response.zeta = 0.5 * response.nu + 0.5 * response.gamma;


%% Use a low pass filter to obtain Calcium dynamics
for i = 1:length(t_index)-1
    x_calcium_raw(:, i+1) = x_calcium_raw(:, i)+stim.dt*(params.lowpass*x_calcium_raw(:, i)+ params.gain*x_scaled(:, i));
end

%% Baseline Correction 
x_calcium_trace = mean(x_calcium_raw(1:stim_A_num, :));
F0 = prctile(x_calcium_trace(floor(0.5*size(x_calcium_trace,2))+1:end),5);
x_calcium_trace(:,1:stim.t_init/stim.dt) = F0; % block pre-stimulus
response.calcium = params.intensity*(x_calcium_trace/F0 - 1);
response.t_index = t_index;

%% Find peaks per trial
response.mat = reshape(response.calcium,[],stim.num_trials)';
response.peaks = peaksDetermine(response.mat, (stim.t_init/stim.dt) + (1:stim.t_on/stim.dt), 1:size(response.mat,2),'max');
response.peaks_norm = response.peaks ./ max(response.peaks);
response.adaptationLevel = nanmean(response.peaks_norm(end-3:end));

end


%% Calc the Riccati eqn
function dXdt = mRiccati_F(t, X, A, B, R, Q)
    X = reshape(X, size(A)); % Convert from "n^2"-by-1 to "n"-by-"n"
    dXdt =( A'*X + X*A - X*B*inv(R)*B'*X + Q); % Determine derivative
    dXdt = dXdt(:); % Convert from "n"-by-"n" to "n^2"-by-1
end

%% Calc the sigmoid fxn
function y = sigmoid(x,c,a)
% sigmoid evaluates a simple sigmoid function along x:
%
%         ______1________
%     y =        -a(x-c)
%          1 + e^
%
% Parse Inputs:
narginchk(1,3)
if nargin<3
    a = 1;
else
    assert(isscalar(a)==1,'a must be a scalar.')
end
if nargin<2
    c = 0;
else
    assert(isscalar(c)==1,'c must be a scalar.')
end
% Perform mathematics:
y = 1./(1 + exp(-a.*(x-c)));

end


%% Weighting matrix
function [ b_Matrix ] = weighting_matrix( m_dim, n_dim, overlap, spread, intensity )
 
    if nargin<5, intensity = 10; end
    if nargin<4, spread = 10; end

    stim_A_num = ceil((overlap + n_dim)/ 2);

    stim_A_mean =  ceil(stim_A_num/2);
    std_A = spread;
    stim_A_weights = pdf('Normal', 1:n_dim, stim_A_mean, std_A);

    stim_B_mean =  n_dim - stim_A_mean + 1; 
    std_B = spread;
    stim_B_weights = pdf('Normal', 1:n_dim, stim_B_mean, std_B);

    b_Matrix = intensity * [stim_A_weights; stim_B_weights];

end

%% Determine peaks
function [Mxx_Values,peakTime] =  peaksDetermine(mat,range,timeVar,type)

    if ~exist('type','var') || isempty(type) type = []; end
    
    if type == 'min'
        for i = 1:size(mat,1)
            [Mxx_Values(i), pidx] = nanmin(mat(i,range));
            Mxx_Values = Mxx_Values';
            peakTime(i) = timeVar(range(pidx));
            peakTime = peakTime';
        end
    elseif type == 'max'
        for i = 1:size(mat,1)
            [Mxx_Values(i), pidx] = nanmax(mat(i,range));
            Mxx_Values = Mxx_Values';
            peakTime(i) = timeVar(range(pidx));
            peakTime = peakTime';
        end
    else
        disp("Error: Please enter a valid {'max','min'} value for the type variable to determine peaks")
    end
end
