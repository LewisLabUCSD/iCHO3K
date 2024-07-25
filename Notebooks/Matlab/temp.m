intervals = {'P0', 'P2'; 'P2', 'P4'; 'P4', 'P6'; 'P6', 'P8'; 'P8', 'P12'; 'P12', 'P14'};
preprocessing = 'raw'; % 'smoothened'

% Ensure results structures are initialized
results = [];
results_for_check = [];
cumulative_results = containers.Map;


% List of batch IDs
batch_ids = {'U1', 'U2', 'U3', 'U4', 'U5', 'U6', 'U7', 'U8'};
% Initialize result table
results = table();
results_for_check = table();

for i = 1:length(batch_ids)
    batch_id = batch_ids{i};
    if isKey(batch_dfs_meta, batch_id)
        % Extract the table from the map
        batch_data = batch_dfs_meta(batch_id);
        
        % Drop the first row
        batch_data(1, :) = [];
        
        % Update the map with the modified table
        batch_dfs_meta(batch_id) = batch_data;
        
        % Get the list of metabolites
        met_list = batch_data.Properties.VariableNames(3:end);

        % Calculate the specific rates for each interval
        for j = 1:size(intervals, 1)
            start = intervals{j, 1};
            end_ = intervals{j, 2};
            interval_key = [start '_' end_];  % Create a unique key for the interval

            end_row = batch_data(strcmp(batch_data.SampleID, end_), :);
            start_row = batch_data(strcmp(batch_data.SampleID, start), :);

            % Proceed only if both rows exist
            if ~isempty(end_row) && ~isempty(start_row)
                try
                    if strcmp(preprocessing, 'smoothened')
                        % Extract necessary values from the smoothened data
                        batch_data_smooth = batch_dfs(batch_id); % Assuming batch_dfs is a map with smoothened data
                        growth_rate = batch_data_smooth(strcmp(batch_data_smooth.SampleID, end_), :).GrowthRate;
                        VC_t = batch_data_smooth(strcmp(batch_data_smooth.SampleID, end_), :).fitted_VCD;
                        VC_i = batch_data_smooth(strcmp(batch_data_smooth.SampleID, start), :).fitted_VCD;
                        IVCD = [];

                        current_time = batch_data_smooth(strcmp(batch_data_smooth.SampleID, end_), :).Age_h_;
                        initial_time = batch_data_smooth(strcmp(batch_data_smooth.SampleID, start), :).Age_h_;
                    elseif strcmp(preprocessing, 'raw')
                        % Extract necessary values from the raw data
                        batch_data_raw = batch_dfs_raw(batch_id); % Assuming batch_dfs_raw is a map with raw data
                        growth_rate = batch_data_raw(strcmp(batch_data_raw.SampleID, end_), :).GrowthRate;
                        IVCD = batch_data_raw(strcmp(batch_data_raw.SampleID, end_), :).VCD;
                        VC_t = [];
                        VC_i = [];

                        current_time = batch_data_raw(strcmp(batch_data_raw.SampleID, end_), :).Age_h_;
                        initial_time = batch_data_raw(strcmp(batch_data_raw.SampleID, start), :).Age_h_;
                    end

                    % Save the growth rate in the cumulative results
                    new_entry = table({batch_id}, {interval_key}, {'SGR'}, growth_rate, ...
                                      'VariableNames', {'BatchID', 'Interval', 'Metabolite', 'SpecificRate'});
                    results = [results; new_entry];

                    % Loop through each metabolite to calculate specific rates
                    for k = 1:length(met_list)
                        met_name = met_list{k};
                        delta_met = end_row.(met_name) - start_row.(met_name);

                        % Calculate specific rate
                        specific_rate = calculate_specific_rate(delta_met, growth_rate, initial_time, current_time, VC_t, VC_i, IVCD);

                        % Save the results
                        new_entry = table({batch_id}, {interval_key}, {met_name}, specific_rate, ...
                                          'VariableNames', {'BatchID', 'Interval', 'Metabolite', 'SpecificRate'});
                        results = [results; new_entry];

                        % Save the check results
                        check_entry = table({num2str(delta_met)}, {num2str(growth_rate)}, {num2str(VC_t)}, ...
                                            {num2str(VC_i)}, {num2str(IVCD)}, {num2str(initial_time)}, {num2str(current_time)}, ...
                                            'VariableNames', {'DeltaMet', 'GrowthRate', 'VC_t', 'VC_i', 'IVCD', 'InitialTime', 'CurrentTime'});
                        results_for_check = [results_for_check; check_entry];
                    end
                catch ME
                    fprintf('Error processing interval %s-%s for batch %s: %s\n', start, end_, batch_id, ME.message);
                end
            end
        end
    end
end

% Display the cumulative results for verification
disp(results);
disp(results_for_check);
