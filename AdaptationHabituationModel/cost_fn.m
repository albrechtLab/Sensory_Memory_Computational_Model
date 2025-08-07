function err = cost_fn(parameters, exp_data, params, stim)

    params.model = parameters;
    response = simulate_neural_activity_fn(params,stim);

    % find peaks
    trial_pts = stim.trial_dur / stim.dt;
    exp_peaks = max(reshape(exp_data, trial_pts,[]));
    exp_peaks_norm = exp_peaks ./max(exp_peaks);

    % calculate error
    err = 0.75 * norm(exp_peaks_norm - response.peaks_norm) ./ sqrt(length(response.peaks)) + ...
          0.25 * norm(exp_data - response.calcium) ./ sqrt(length(response.calcium));

end