function set_publication_figure_style(fig_handle, varargin)
% Apply consistent typography and formatting to figures.

    p = inputParser;
    addParameter(p, 'font_size', 11, @isnumeric);
    addParameter(p, 'font_name', 'Arial', @ischar);
    parse(p, varargin{:});

    font_size = p.Results.font_size;
    font_name = p.Results.font_name;

    ax_all = findall(fig_handle, 'Type', 'axes');
    for i = 1:length(ax_all)
        ax = ax_all(i);
        set(ax, 'FontName', font_name, 'FontSize', font_size);
        set(ax, 'Box', 'on', 'TickDir', 'in');
        set(ax, 'GridAlpha', 0.3);
    end

    text_all = findall(fig_handle, 'Type', 'text');
    for i = 1:length(text_all)
        set(text_all(i), 'FontName', font_name, 'FontSize', font_size);
    end

    set(fig_handle, 'Color', 'white');
end
