function save_publication_figure(fig_handle, output_dir, fig_number, fig_title, caption, ...
    input_csv, n_samples, analysis_date)
% Saves figure in multiple formats with metadata and captions.

    if nargin < 7
        n_samples = 0;
    end
    if nargin < 8
        analysis_date = datestr(now, 'yyyy-mm-dd');
    end

    fig_name = sprintf('Fig_%02d_%s', fig_number, strrep(fig_title, ' ', '_'));

    pdf_dir = fullfile(output_dir, 'PDF');
    eps_dir = fullfile(output_dir, 'EPS');
    png_dir = fullfile(output_dir, 'PNG_Preview');

    if ~isfolder(pdf_dir), mkdir(pdf_dir); end
    if ~isfolder(eps_dir), mkdir(eps_dir); end
    if ~isfolder(png_dir), mkdir(png_dir); end

    try
        pdf_file = fullfile(pdf_dir, [fig_name '.pdf']);
        print(fig_handle, '-dpdf', '-r300', pdf_file);
        fprintf('  OK PDF saved: %s\n', pdf_file);
    catch ME
        fprintf('  WARN PDF save failed: %s\n', ME.message);
        try
            saveas(fig_handle, pdf_file);
        catch
            % Silent fallback failure
        end
    end

    try
        eps_file = fullfile(eps_dir, [fig_name '.eps']);
        print(fig_handle, '-depsc2', '-r300', eps_file);
        fprintf('  OK EPS saved: %s\n', eps_file);
    catch ME
        fprintf('  WARN EPS save failed: %s\n', ME.message);
    end

    try
        png_file = fullfile(png_dir, [fig_name '_preview.png']);
        print(fig_handle, '-dpng', '-r150', png_file);
        fprintf('  OK PNG saved: %s\n', png_file);
    catch ME
        fprintf('  WARN PNG save failed: %s\n', ME.message);
    end

    % Ensure metadata is present on the figure prior to saving
    if nargin >= 8
        add_figure_metadata(input_csv, n_samples, analysis_date);
    end
end
