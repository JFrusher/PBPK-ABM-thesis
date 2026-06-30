function startup()
% Add all 5FU_PBPK subdirectories to path.
% Call once before running run5FU_PBPK_Simulation.
    base = fileparts(mfilename('fullpath'));
    subdirs = {'params','validation','dosing','pk_core','pd_model','analysis','output','utils'};
    addpath(base);
    for i = 1:numel(subdirs)
        addpath(fullfile(base, subdirs{i}));
    end
end
