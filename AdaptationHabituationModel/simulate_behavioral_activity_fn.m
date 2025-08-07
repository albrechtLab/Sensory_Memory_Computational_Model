function behavior = simulate_behavioral_activity_fn(params, response)
%
% function behavior = simulate_behavioral_activity_fn(params, response)
%
%	behavior.  	forward, pause, reverse, sens_on_prob, sens_off_prob 	
%
%	response.  	calcium, x, nu, gamma 					{from simulate_neural_activity}
%
%	params.		model [drift_nu, drift_nu_gamma, tau_nu, drift_gamma, tau_gamma, beta]  {optimization}
%			    q, r, s,						{penalty values}
%			    m_dim, n_dim, overlap,			{dimensions}
%			    c0, l1, a1, l2, a2 				{z_c}
%			    baseline, low_asymp, up_asymp, curve, lowpass, gain, intensity		{calcium transform}
%			    SEED, noise_level						
%			    var. forward, pause, reverse
%			    mat. transition, valence
%			    weights, stim_on_mean, stim_off_mean
%

  %% Set parameters 
  m_dim = params.m_dim;
  n_dim = params.n_dim;
      
  if (false)

    params.var.forward = 0.15; 
    params.var.pause = 0.2;
    params.var.reverse = 0.3; 
    
    params.stim_on_mean = [1; 0];
    params.stim_off_mean = [0; 0];
    
    params.mat.transition = [0.75, 0.15, 0.1; 
                              0.1, 0.7, 0.2; 
                             0.15, 0.15, 0.75]; % or 0.25, 0.15, 0.6];
                  
    params.mat.valence = [0.3,  0.6;    % or [0.3, 0.75;
                          0.35, 0.2;    % or 0.35, 0.125;
                          0.35, 0.2];   % or 0.35, 0.125];
               
    params.weights = [0.5, 0.5]; % or [0.75, 0.25];

  end
              
    
  %% Additional visualizations
    
  odor_on_forward = mvnrnd(params.stim_on_mean, params.var.forward*eye(m_dim), 100);
  odor_on_pause = mvnrnd(params.stim_on_mean, params.var.pause*eye(m_dim), 100);
  odor_on_reverse = mvnrnd(params.stim_on_mean, params.var.reverse*eye(m_dim), 100);
    
  odor_off_forward = mvnrnd(params.stim_off_mean, params.var.forward*eye(m_dim), 100);
  odor_off_pause = mvnrnd(params.stim_off_mean, params.var.pause*eye(m_dim), 100);
  odor_off_reverse = mvnrnd(params.stim_off_mean, params.var.reverse*eye(m_dim), 100);
    
    
  %% Initialize vectors
    
  t_index = response.t_index;
  cum_sensory_input = zeros(m_dim, length(t_index));
    
  % Initialize for 3 behavioral states (forward, pause, reverse)
  prob_sensory_off = zeros(3, length(t_index));
  prob_sensory_on = zeros(3, length(t_index));
    
  behavior.prob = zeros(3, length(t_index));
  behavior.prob(:,1) = 1/3;
    
  for i = 2:length(t_index) 
    cum_sensory_input(:,i) = params.weights(1)*response.nu(:,i) + params.weights(2)*response.gamma(:,i) + ...
       	     mvnrnd([0,0],params.noise_level*eye(m_dim))';
    
    prob_sensory_off(1, i) = mvnpdf(cum_sensory_input(:,i), params.stim_off_mean, params.var.forward*eye(m_dim));
    prob_sensory_off(2, i) = mvnpdf(cum_sensory_input(:,i), params.stim_off_mean, params.var.pause*eye(m_dim));
    prob_sensory_off(3, i) = mvnpdf(cum_sensory_input(:,i), params.stim_off_mean, params.var.reverse*eye(m_dim));
        
    prob_sensory_on(1, i) = mvnpdf(cum_sensory_input(:,i), params.stim_on_mean, params.var.forward*eye(m_dim));
    prob_sensory_on(2, i) = mvnpdf(cum_sensory_input(:,i), params.stim_on_mean, params.var.pause*eye(m_dim));
    prob_sensory_on(3, i) = mvnpdf(cum_sensory_input(:,i), params.stim_on_mean, params.var.reverse*eye(m_dim));
        
    behavior.prob(1,i) = ((params.mat.valence(1,1)*prob_sensory_off(1,i) + params.mat.valence(1,2)*prob_sensory_on(1,i)) / ...
                        sum(params.mat.valence(1,:)))*params.mat.transition(1,:)*behavior.prob(:,i-1);
    
    behavior.prob(2,i) = ((params.mat.valence(2,1)*prob_sensory_off(2,i) + params.mat.valence(2,2)*prob_sensory_on(2,i)) / ...
                        sum(params.mat.valence(2,:)))*params.mat.transition(2,:)*behavior.prob(:,i-1);
    
    behavior.prob(3,i) = ((params.mat.valence(3,1)*prob_sensory_off(3,i) + params.mat.valence(3,2)*prob_sensory_on(3,i)) / ...
                        sum(params.mat.valence(3,:)))*params.mat.transition(3,:)*behavior.prob(:,i-1);
  
    
    probsum= sum(behavior.prob(:,i));    
    behavior.prob(1,i) = behavior.prob(1,i) / probsum;
    behavior.prob(2,i) = behavior.prob(2,i) / probsum;     
    behavior.prob(3,i) = behavior.prob(3,i) / probsum;

  end
    
  behavior.forward = behavior.prob(1, :);
  behavior.pause = behavior.prob(2, :);
  behavior.reverse = behavior.prob(3, :);

  % Downsample to match experiment data for behavior
  behavior.prob_ds = squeeze(median(reshape(behavior.prob,3,20,[]),2)); % downsample by 20x
  
end