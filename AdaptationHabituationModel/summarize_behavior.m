function beh_index = summarize_behavior(behavior,stim)

    % input behavior can be either model output structure, OR experiment
    % data as a linear vector of forward probability
    if (isstruct(behavior)) 
        fwdMat = reshape(behavior.prob_ds(1,:),[],stim.num_trials)';
    else
        fwdMat = reshape(behavior,[],stim.num_trials)';
    end

    % Calculate Behavior index and Habituation index for downsampled data
    [onBins, offBins, baseBins] = getOnOffBins(stim);
    [BehIdxMat,BehIdx,RespIdx,HabIdx] = behaviorIndices(fwdMat, baseBins, onBins, offBins);
    beh_index.mat = BehIdxMat;
    beh_index.linear = BehIdx;
    beh_index.per_trial = RespIdx;
    beh_index.habituation = HabIdx;
    beh_index.bins.on = onBins;
    beh_index.bins.off = offBins;
    beh_index.bins.base = baseBins;

end

%% Determine time bins for response classification 
function [onBins, offBins, baseBins] = getOnOffBins(stim)
    t_mat = 1:stim.dt_ds:stim.trial_dur;

    onBins = find((t_mat > stim.t_init) & (t_mat <= stim.t_init + stim.t_on)); %add 1 since t_mat is the ending time in s
    offBins = find((t_mat > stim.t_init + stim.t_on + 4) & (t_mat < stim.trial_dur - 6)); 
    baseBins = [1,2,length(t_mat)+[-1,0]];

    % augment as needed to ensure averaging at least 4 off bins
    offvals = length(offBins);
    if offvals < 4 % augment extra time bins as needed
        offBins = [offBins, baseBins(1:(4-offvals))];
    end
end

%% Calculate behavior indices 
function [BehIdxMat,BehIdx,RespIdx,HabIdx] = behaviorIndices(fwdMat, baseBins, onBins, offBins)

    % Behavior Index (instantaneous over time) = fwd fraction subtracting baseline per trial, defined by baseBins 
    bases = nanmean(fwdMat(:,baseBins),2);  % baselines per trial
    BehIdxMat = fwdMat - repmat(bases,1,size(fwdMat,2));
    BehIdx = reshape(BehIdxMat',[],1)';
    
    % Response Index (per trial) = Behavior Index in stimulus (onBins) - outside
    % stimulus (offBins)
    RespIdx = nanmean(BehIdxMat(:,onBins),2) - nanmean(BehIdxMat(:,offBins),2);
    
    %HabIdx = nanmean(RespIdx(end-3:end)) / RespIdx(1); % compare last 4 to first trial
    HabIdx = nanmean(RespIdx(end-3:end)) / nanmax(RespIdx(1:5)); % compare last 4 to first trial
end