% Test Logger in parfor context to isolate the "Dot indexing" error
%
% This script tests whether the Logger class works correctly when called
% from within a parfor loop, which is the context where MC simulations fail.

fprintf('Testing Logger in parfor context...\n');

% Test 1: Logger creation
fprintf('\n=== Test 1: Creating Logger ===\n');
try
    logConfig = struct('logDir', fullfile(pwd, 'logs'), 'level', 'INFO', 'enableDebug', false);
    logger = Logger(logConfig.logDir, logConfig.level, logConfig.enableDebug);
    fprintf('✓ Logger created successfully\n');
catch ME
    fprintf('✗ Logger creation FAILED: %s\n', ME.message);
    fprintf('  Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('    %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
    return;
end

% Test 2: Simple logging
fprintf('\n=== Test 2: Simple Logging ===\n');
try
    logger.info('test_event', struct('value', 42, 'name', 'test'));
    fprintf('✓ Simple logging works\n');
catch ME
    fprintf('✗ Simple logging FAILED: %s\n', ME.message);
    fprintf('  Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('    %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
    return;
end

% Test 3: Logging in parfor (the critical test)
fprintf('\n=== Test 3: Logging in parfor ===\n');
try
    results = cell(4, 1);
    parfor i = 1:4
        try
            % Create a unique logger for each worker to avoid file conflicts
            workerLogDir = fullfile(pwd, 'logs', sprintf('worker_%d', i));
            if ~exist(workerLogDir, 'dir')
                mkdir(workerLogDir);
            end
            workerLogger = Logger(workerLogDir, 'INFO', false);
            
            % Try logging some data
            workerLogger.info('parfor_test', struct('iteration', i, 'value', i*10));
            
            results{i} = 'SUCCESS';
        catch ME
            results{i} = sprintf('FAILED: %s at %s:%d', ME.message, ME.stack(1).name, ME.stack(1).line);
        end
    end
    
    % Check results
    failures = 0;
    for i = 1:4
        if contains(results{i}, 'FAILED')
            fprintf('  Worker %d: %s\n', i, results{i});
            failures = failures + 1;
        else
            fprintf('  Worker %d: SUCCESS\n', i);
        end
    end
    
    if failures == 0
        fprintf('✓ All parfor workers logged successfully\n');
    else
        fprintf('✗ %d/%d parfor workers failed\n', failures, 4);
    end
catch ME
    fprintf('✗ parfor test FAILED: %s\n', ME.message);
    fprintf('  Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('    %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
end

fprintf('\n=== Test Complete ===\n');
fprintf('If you see "Dot indexing not supported" errors, the issue is with Logger in parfor.\n');
fprintf('Otherwise, the issue is somewhere else in the simulation code.\n');
