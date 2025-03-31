% Calculate specific rate function
function SR = calculate_specific_rate(delta_met, growth_rate, initial_time, current_time, VC_t, VC_i, IVCD)
    try
        time_interval = current_time - initial_time;
        if ~isempty(VC_t) && ~isempty(VC_i) && isempty(IVCD)
            cell_delta = VC_t - VC_i;
        elseif isempty(VC_t) && isempty(VC_i) && ~isempty(IVCD)
            cell_delta = IVCD;
        else
            cell_delta = NaN; % Default to NaN if none of the conditions are met
        end

        % Zero division handling for `cell_delta`
        if cell_delta == 0
            SR = 0;
            return;
        end

        % Zero division handling for other parameters
        if time_interval == 0 || growth_rate == 0
            SR = NaN;
            return;
        end

        SR = delta_met / cell_delta * growth_rate / (200 / 1000);
    catch ME
        fprintf('Error calculating specific rate: %s\n', ME.message);
        SR = NaN;
    end
end