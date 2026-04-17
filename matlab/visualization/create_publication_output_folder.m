function output_folder = create_publication_output_folder(base_output_dir)
% Creates organized folder structure for publication figures.

    output_folder = fullfile(base_output_dir, 'Figures_Publication');

    if ~isfolder(output_folder)
        mkdir(output_folder);
    end

    subdirs = {'PDF', 'EPS', 'PNG_Preview', 'README'};
    for i = 1:length(subdirs)
        subdir_path = fullfile(output_folder, subdirs{i});
        if ~isfolder(subdir_path)
            mkdir(subdir_path);
        end
    end

    readme_file = fullfile(output_folder, 'README', 'FIGURE_ORGANIZATION.txt');
    fid = fopen(readme_file, 'w');
    fprintf(fid, 'PUBLICATION FIGURE ORGANIZATION\n');
    fprintf(fid, '================================\n\n');
    fprintf(fid, 'PDF/:    Use for journal submission (vector graphics, 300 dpi)\n');
    fprintf(fid, 'EPS/:    Use for editing in Office applications\n');
    fprintf(fid, 'PNG_Preview/: Use for quick preview and presentations\n');
    fprintf(fid, '\nAll files include:\n');
    fprintf(fid, '  - Input dosing scenario name\n');
    fprintf(fid, '  - Number of Monte Carlo simulations\n');
    fprintf(fid, '  - Analysis date\n');
    fprintf(fid, '  - Standardized Arial font (11 pt)\n');
    fprintf(fid, '  - Publication-ready captions\n');
    fclose(fid);

    fprintf('Publication output folder created: %s\n', output_folder);
end
