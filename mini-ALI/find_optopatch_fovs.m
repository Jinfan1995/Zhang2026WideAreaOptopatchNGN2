function fov_table = find_optopatch_fovs(root_path)

%% Initialize output columns
plate_col = strings(0,1);
well_col = strings(0,1);
fov_col = uint16([]);
fov_name_col = strings(0,1);
fov_path_col = strings(0,1);


%% =========================================================
% Find plate folders
% ==========================================================

plate_dirs = dir(root_path);
plate_dirs = plate_dirs([plate_dirs.isdir]);

% Remove '.', '..', and Analysis
plate_dirs = plate_dirs(~ismember( ...
    {plate_dirs.name}, ...
    {'.','..','Analysis'}));


%% =========================================================
% Loop through plates
% ==========================================================

for p = 1:length(plate_dirs)

    plate_name = string(plate_dirs(p).name);

    plate_path = fullfile( ...
        plate_dirs(p).folder, ...
        plate_dirs(p).name);

    fprintf('\nChecking plate: %s\n', plate_name);


    %% =====================================================
    % Find well folders
    % ======================================================

    well_dirs = dir(plate_path);
    well_dirs = well_dirs([well_dirs.isdir]);

    well_dirs = well_dirs(~ismember( ...
        {well_dirs.name}, ...
        {'.','..'}));

    fprintf('  Found %d wells\n', length(well_dirs));


    %% =====================================================
    % Loop through wells
    % ======================================================

    for w = 1:length(well_dirs)

        well_folder_name = string(well_dirs(w).name);

        well_path = fullfile( ...
            well_dirs(w).folder, ...
            well_dirs(w).name);

        % Convert "Well_A1" -> "A1"
        well_name = erase(well_folder_name, "Well_");


        %% =================================================
        % Find ALL subfolders in this well
        % ==================================================

        movie_dirs = dir(well_path);
        movie_dirs = movie_dirs([movie_dirs.isdir]);

        movie_dirs = movie_dirs(~ismember( ...
            {movie_dirs.name}, ...
            {'.','..'}));


        %% =================================================
        % Identify valid movie folders
        % ==================================================

        valid_movie_dirs = movie_dirs([]);

        for m = 1:length(movie_dirs)

            candidate_path = fullfile( ...
                movie_dirs(m).folder, ...
                movie_dirs(m).name);

            movie_file = fullfile( ...
                candidate_path, ...
                "frames1.bin");

            metadata_file = fullfile( ...
                candidate_path, ...
                "output_data.mat");

            % Only treat folder as FOV if required files exist
            if isfile(movie_file) && isfile(metadata_file)

                valid_movie_dirs(end+1) = movie_dirs(m); %#ok<AGROW>

            end

        end


        fprintf( ...
            '    %s: found %d valid movie FOVs\n', ...
            well_folder_name, ...
            length(valid_movie_dirs));


        %% =================================================
        % Assign FOV number based on folder order
        % ==================================================

        for f = 1:length(valid_movie_dirs)

            fov_folder_name = ...
                string(valid_movie_dirs(f).name);

            fov_path = fullfile( ...
                valid_movie_dirs(f).folder, ...
                valid_movie_dirs(f).name);


            % FOV number is simply its order within the well
            fov_number = uint16(f);


            %% Add to output
            plate_col(end+1,1) = plate_name;

            well_col(end+1,1) = well_name;

            fov_col(end+1,1) = fov_number;

            fov_name_col(end+1,1) = ...
                fov_folder_name;

            fov_path_col(end+1,1) = ...
                string(fov_path);

        end

    end

end


%% =========================================================
% Build final table
% ==========================================================

fov_table = table( ...
    plate_col, ...
    well_col, ...
    fov_col, ...
    fov_name_col, ...
    fov_path_col, ...
    'VariableNames', { ...
        'plate', ...
        'well', ...
        'fov', ...
        'fov_name', ...
        'fov_path'});

end