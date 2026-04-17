function [samples, completed_runs, checkpoint_file] = loadCheckpoint(output_dir, n_samples)
    % loadCheckpoint  Load checkpoint data if available, to resume interrupted runs
    %
    % Returns:
    %   samples          - LHS sample matrix (n_samples x n_params)
    %   completed_runs   - Logical array indicating which runs completed
    %   checkpoint_file  - Path to checkpoint file (for saving later)
    
    checkpoint_file = fullfile(output_dir, 'MC_checkpoint.mat');
    
    if isfile(checkpoint_file)
        load(checkpoint_file, 'samples', 'completed_runs');
        fprintf('MC: Resuming from checkpoint. %d/%d runs already completed.\n', sum(completed_runs), n_samples);
    else
        samples = [];
        completed_runs = false(n_samples, 1);
        fprintf('MC: No checkpoint found. Starting fresh run.\n');
    end
end

function saveCheckpoint(samples, completed_runs, checkpoint_file)
    % saveCheckpoint  Save current progress to checkpoint file
    
    try
        save(checkpoint_file, 'samples', 'completed_runs', '-v7.3');
    catch ME
        warning('MC: Could not save checkpoint: %s', ME.message);
    end
end
