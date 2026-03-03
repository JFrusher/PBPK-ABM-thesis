function ntfyNotify(title, message, tags, priority)
    % ntfyNotify  Send a notification via ntfy.sh service
    %
    % Usage:
    %   ntfyNotify('My Title', 'My message')
    %   ntfyNotify('Alert', 'Error occurred!', 'error', 'high')
    %
    % Args:
    %   title    - Notification title (will be used as ntfy topic)
    %   message  - Notification message body
    %   tags     - (optional) Tags for the notification (e.g., 'red_circle,warning')
    %   priority - (optional) Priority: 'low'|'default'|'high'|'urgent'
    %
    % Notes:
    %   - Requires internet connection
    %   - Topics are public on ntfy.sh (use a unique UUID for privacy)
    %   - To receive on phone: visit https://ntfy.sh/<topic>
    
    if nargin < 2
        return; % silently ignore if not enough args
    end
    if nargin < 3 || isempty(tags)
        tags = '';
    end
    if nargin < 4 || isempty(priority)
        priority = 'default';
    end
    
    try
        % Sanitize topic name (replace spaces, special chars)
        topic = regexprep(title, '[^a-zA-Z0-9_-]', '_');
        topic = lower(topic);
        
        % Build URL
        url = sprintf('https://ntfy.sh/%s', topic);
        
        % Prepare headers
        options = weboptions('Timeout', 5, 'MediaType', 'application/json');
        
        % Build JSON payload
        payload = struct('title', title, 'message', message);
        
        if ~isempty(tags)
            payload.tags = tags;
        end
        
        % Map priority to ntfy priority levels
        priority_map = containers.Map(...
            {'low', 'default', 'high', 'urgent'}, ...
            {1, 2, 4, 5});
        
        if priority_map.isKey(lower(priority))
            payload.priority = priority_map(lower(priority));
        else
            payload.priority = 2; % default
        end
        
        % Send via webwrite
        webwrite(url, payload, options);
        fprintf('[ntfy] Sent: %s | %s\n', title, message);
        
    catch ME
        % Silently fail for network issues - don't break the main script
        fprintf('[ntfy] Warning: Could not send notification (%s). Continuing...\n', ME.message);
    end
end
