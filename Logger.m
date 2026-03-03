classdef Logger < handle
    % Lightweight structured text logger with sampling and levels
    % Outputs in standard format: YYYY-MM-DD HH:mm:ss.SSS [LEVEL] EventName: key1=value1 key2=value2
    properties
        logFile
        minLevel % string 'DEBUG'|'INFO'|'WARN'|'ERROR'|'FATAL'
        enableDebug logical = false;
        countsMap % containers.Map for sampled messages
        enabled logical = true;           % set to false if logfile cannot be written
        fallbacked logical = false;      % true if we fell back to tempdir
        hasWarnedOpenFailure logical = false; % avoid flooding warnings
    end
    methods
        function obj = Logger(logDir, level, enableDebug)
            if nargin < 1 || isempty(logDir)
                logDir = fullfile(pwd, 'logs');
            end
            % Try to create log directory, but fall back safely on failure
            try
                if ~exist(logDir, 'dir')
                    mkdir(logDir);
                end
            catch ME
                warning('Logger:DirCreateFailed', 'Could not create log directory %s: %s. Will attempt to use tempdir()', logDir, ME.message);
            end

            timestamp = datestr(now, 'yyyymmdd_HHMMSS');
            candidate = fullfile(logDir, sprintf('run5FU_%s.log', timestamp));
            % Attempt to open candidate file to validate writability
            fid = -1;
            try
                fid = fopen(candidate, 'a');
            catch
                fid = -1;
            end
            if fid ~= -1
                fclose(fid);
                obj.logFile = candidate;
                obj.fallbacked = false;
                obj.enabled = true;
            else
                % fallback to system tempdir
                fallback = fullfile(tempdir, sprintf('run5FU_%s.log', timestamp));
                fid2 = -1;
                try
                    fid2 = fopen(fallback, 'a');
                catch
                    fid2 = -1;
                end
                if fid2 ~= -1
                    fclose(fid2);
                    obj.logFile = fallback;
                    obj.fallbacked = true;
                    obj.enabled = true;
                    warning('Logger:FallbackToTemp','Could not open log file in %s, falling back to %s', logDir, fallback);
                else
                    % cannot write anywhere; disable logger to avoid repeated noisy warnings
                    obj.logFile = '';
                    obj.enabled = false;
                    obj.fallbacked = false;
                    warning('Logger:Disabled','Logger disabled: unable to open log files in %s or tempdir()', logDir);
                end
            end

            if nargin < 2 || isempty(level)
                level = 'INFO';
            end
            obj.minLevel = upper(level);
            if nargin >= 3
                obj.enableDebug = logical(enableDebug);
            else
                obj.enableDebug = false;
            end
            obj.countsMap = containers.Map('KeyType','char','ValueType','double');

            % Try to write an initial header entry (may be disabled if no write access)
            try
                if obj.enabled
                    % Directly write log line to avoid calling obj.info during construction
                    log_line = sprintf('%s [INFO] logger_init: logFile=%s minLevel=%s', ...
                        datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF'), obj.logFile, obj.minLevel);
                    fid = fopen(obj.logFile, 'a');
                    if fid ~= -1
                        fprintf(fid, '%s\n', log_line);
                        fclose(fid);
                    end
                end
            catch
                % swallow initialization logging errors
            end
        end
        function p = levelPriority(~, level)
            switch upper(string(level))
                case 'DEBUG'
                    p = 0;
                case 'INFO'
                    p = 1;
                case 'WARN'
                    p = 2;
                case 'ERROR'
                    p = 3;
                case 'FATAL'
                    p = 4;
                otherwise
                    p = 1;
            end
        end
        function enabled = isDebugEnabled(obj)
            enabled = obj.enableDebug && (obj.levelPriority('DEBUG') >= obj.levelPriority(obj.minLevel));
        end
        function enabled = isInfoEnabled(obj)
            enabled = obj.levelPriority('INFO') >= obj.levelPriority(obj.minLevel);
        end
        function log(obj, level, key, data)
            % No-op if logging disabled
            if ~obj.enabled
                return;
            end
            try
                if obj.levelPriority(level) < obj.levelPriority(obj.minLevel)
                    return; % below threshold
                end
                
                % Format timestamp: YYYY-MM-DD HH:mm:ss.SSS
                now_time = datetime('now');
                timestamp = sprintf('%04d-%02d-%02d %02d:%02d:%02d.%03d', ...
                    year(now_time), month(now_time), day(now_time), ...
                    hour(now_time), minute(now_time), floor(second(now_time)), ...
                    round(1000 * (second(now_time) - floor(second(now_time)))));
                
                % Format level in brackets
                level_str = sprintf('[%s]', upper(level));
                
                % Format payload as flattened key=value pairs
                payload_str = obj.flattenPayload(data);
                
                % Construct log line
                if isempty(payload_str)
                    log_line = sprintf('%s %s %s:', timestamp, level_str, key);
                else
                    log_line = sprintf('%s %s %s: %s', timestamp, level_str, key, payload_str);
                end
                
                % Append to log file
                fid = fopen(obj.logFile, 'a');
                if fid ~= -1
                    fprintf(fid, '%s\n', log_line);
                    fclose(fid);
                else
                    % Avoid repeatedly warning on inability to open the file
                    if ~obj.hasWarnedOpenFailure
                        warning('Logger:OpenFailed','Unable to open log file %s; disabling logger to avoid flood.', obj.logFile);
                        obj.hasWarnedOpenFailure = true;
                    end
                    % disable logging to prevent repeated failures
                    obj.enabled = false;
                end
            catch ME
                % On unexpected error, emit a single warning and disable
                if ~obj.hasWarnedOpenFailure
                    warning('Logger:WriteError', 'Logger write failed: %s. Disabling logger.', ME.message);
                    obj.hasWarnedOpenFailure = true;
                end
                obj.enabled = false;
            end
        end
        
        function payload_str = flattenPayload(~, data)
            % Convert struct/table data into flattened key=value string format
            % Rounds floating-point numbers to 4 decimal places
            payload_str = '';
            
            % Robust error handling - if anything fails, fall back gracefully
            try
                if isempty(data)
                    return;
                end
                
                if isstruct(data)
                    field_names = fieldnames(data);
                    pairs = cell(1, length(field_names));
                    for i = 1:length(field_names)
                        try
                            field = field_names{i};
                            value = data.(field);
                            pairs{i} = Logger.formatKeyValue(field, value);
                        catch
                            % If formatting fails for one field, use simple representation
                            pairs{i} = sprintf('%s=<error>', field_names{i});
                        end
                    end
                    payload_str = strjoin(pairs(~cellfun(@isempty, pairs)), ' ');
                elseif istable(data)
                    try
                        var_names = data.Properties.VariableNames;
                        pairs = cell(1, length(var_names));
                        for i = 1:length(var_names)
                            try
                                var = var_names{i};
                                value = data{1, i};
                                pairs{i} = Logger.formatKeyValue(var, value);
                            catch
                                pairs{i} = sprintf('%s=<error>', var_names{i});
                            end
                        end
                        payload_str = strjoin(pairs(~cellfun(@isempty, pairs)), ' ');
                    catch
                        payload_str = '<table_format_error>';
                    end
                else
                    % Unknown type - convert to string
                    try
                        payload_str = sprintf('%s', class(data));
                    catch
                        payload_str = '<unknown>';
                    end
                end
            catch ME
                % Ultimate fallback - just note there was data
                payload_str = sprintf('<format_error:%s>', ME.message);
            end
        end
    end
    methods (Static)
        function pair = formatKeyValue(key, value)
            % Format a single key=value pair with comprehensive error handling
            try
                if ischar(value) || isstring(value)
                    pair = sprintf('%s=%s', key, value);
                elseif isnumeric(value) && isscalar(value)
                    try
                        if floor(value) == value
                            % Integer value
                            pair = sprintf('%s=%d', key, value);
                        else
                            % Floating-point: round to 4 decimals
                            rounded = round(value, 4);
                            pair = sprintf('%s=%.4g', key, rounded);
                        end
                    catch
                        % Fallback for numeric formatting issues
                        pair = sprintf('%s=%s', key, mat2str(value));
                    end
                elseif islogical(value)
                    pair = sprintf('%s=%s', key, string(value));
                elseif isnumeric(value)
                    % Array - try to format safely
                    try
                        pair = sprintf('%s=%s', key, mat2str(value));
                    catch
                        pair = sprintf('%s=[%dx%d]', key, size(value, 1), size(value, 2));
                    end
                else
                    % Complex types - try multiple strategies
                    try
                        pair = sprintf('%s=%s', key, jsonencode(value));
                    catch
                        try
                            pair = sprintf('%s=%s', key, class(value));
                        catch
                            pair = sprintf('%s=<unknown>', key);
                        end
                    end
                end
            catch ME
                % Ultimate fallback
                pair = sprintf('%s=<error:%s>', key, ME.message);
            end
        end
    end
    
    % Instance methods for logging convenience
    methods
        function debug(obj, key, data)
            if nargin < 3
                data = struct();
            end
            obj.log('DEBUG', key, data);
        end
        function info(obj, key, data)
            if nargin < 3
                data = struct();
            end
            obj.log('INFO', key, data);
        end
        function warn(obj, key, data)
            if nargin < 3
                data = struct();
            end
            obj.log('WARN', key, data);
        end
        function error(obj, key, data)
            if nargin < 3
                data = struct();
            end
            obj.log('ERROR', key, data);
        end
        function fatal(obj, key, data)
            if nargin < 3
                data = struct();
            end
            obj.log('FATAL', key, data);
        end
        function warnSampled(obj, key, data, maxCount)
            if nargin < 4
                maxCount = 5;
            end
            if ~ischar(key)
                key = jsonencode(key);
            end
            if obj.countsMap.isKey(key)
                c = obj.countsMap(key) + 1;
            else
                c = 1;
            end
            obj.countsMap(key) = c;
            if c <= maxCount
                obj.warn(key, struct('occurrence', c, 'detail', data));
            end
        end
    end
end