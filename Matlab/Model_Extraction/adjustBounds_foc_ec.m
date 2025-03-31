%% Function to Adjust Bounds if Infeasible
function [bounds_dict, time_field] = adjustBounds_foc_ec(bounds_dict, rxn, relaxation_factor, time)
    % Adjust the bounds for reaction rxn by relaxation_factor
    original_bounds = bounds_dict.(rxn);
    lb = original_bounds.(time)(1);
    ub = original_bounds.(time)(2);
    time_field = matlab.lang.makeValidName(time);
    % Adjust lower bound
    if lb < 0
        lb = lb * relaxation_factor; % Make lb more negative
    elseif lb > 0
        lb = lb / relaxation_factor; % Make lb smaller
    else
        % lb == 0, set to small negative value
        lb = -0.1;
    end

    % Adjust upper bound
    if ub > 0
        ub = ub * relaxation_factor; % Increase ub
    elseif ub < 0
        ub = ub / relaxation_factor; % Decrease ub
    else
        % ub == 0, set to small positive value
        ub = 0.1;
    end

    % Update the new_bounds variable with adjusted values
    new_bounds = bounds_dict.(rxn); % Start with existing fields
    new_bounds.(time_field) = [lb, ub]; % Add or update the time field

    % Store the adjusted bounds in bounds_dict
    bounds_dict.(rxn) = new_bounds;
    % Add time field to all other reactions if they exist
    all_rxns = fieldnames(bounds_dict);
    for i = 1:length(all_rxns)
        other_rxn = all_rxns{i};
            if ~strcmp(other_rxn, rxn)
            bounds_dict.(other_rxn).(time_field) = bounds_dict.(other_rxn).(time); % Add the time field to other reactions
            end
    end

    % Print adjusted bounds for verification
    fprintf('Adjusted bounds for reaction %s at time %s: Lower Bound = %f, Upper Bound = %f\n', rxn, time, lb, ub);
end