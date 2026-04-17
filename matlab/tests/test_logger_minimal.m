% Minimal test to isolate Logger issue
fprintf('=== MINIMAL LOGGER TEST ===\n\n');

% Test 1: Create Logger with minimal args
fprintf('Test 1: Creating Logger with defaults...\n');
try
    logger1 = Logger();
    fprintf('  Logger created: %s\n', class(logger1));
    fprintf('  Is handle: %d\n', isa(logger1, 'handle'));
    fprintf('  Has logFile: %d\n', isprop(logger1, 'logFile'));
    fprintf('  logFile = %s\n', logger1.logFile);
catch ME
    fprintf('  FAILED: %s\n', ME.message);
    for i = 1:length(ME.stack)
        fprintf('    %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
end

% Test 2: Call info method
fprintf('\nTest 2: Calling logger.info()...\n');
try
    logger1.info('test_event', struct('value', 42));
    fprintf('  SUCCESS\n');
catch ME
    fprintf('  FAILED: %s\n', ME.message);
    fprintf('  logger1 type: %s\n', class(logger1));
    fprintf('  logger1 properties:\n');
    disp(logger1);
    for i = 1:length(ME.stack)
        fprintf('    %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
end

% Test 3: Create Logger as simulation does
fprintf('\nTest 3: Creating Logger as simulation does...\n');
try
    logConfig = struct();
    defaultLogDir = fullfile(pwd, 'logs');
    logConfig.logDir = defaultLogDir;
    logConfig.level = 'INFO';
    logConfig.enableDebug = false;
    
    fprintf('  logConfig fields: ');
    disp(fieldnames(logConfig));
    fprintf('  logConfig.logDir = %s\n', logConfig.logDir);
    fprintf('  logConfig.level = %s\n', logConfig.level);
    
    logger2 = Logger(logConfig.logDir, logConfig.level, logConfig.enableDebug);
    fprintf('  Logger created: %s\n', class(logger2));
    fprintf('  Calling logger2.info()...\n');
    logger2.info('simulation_initializing', struct('inputFile', 'test.csv', 'outputPrefix', 'test'));
    fprintf('  SUCCESS\n');
catch ME
    fprintf('  FAILED: %s\n', ME.message);
    fprintf('  At: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
end

fprintf('\n=== TEST COMPLETE ===\n');
