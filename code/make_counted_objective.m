function [counted_fobj, get_state] = make_counted_objective(base_fobj, max_fes)
%MAKE_COUNTED_OBJECTIVE Wrap an objective function with a strict FE budget.
%
%   [FCOUNTED, GET_STATE] = MAKE_COUNTED_OBJECTIVE(FBASE, MAX_FES)
%   returns a counted objective handle and a state-reader handle. Calls made
%   after MAX_FES return Inf without evaluating FBASE or increasing the count.

fe_count = 0;
best_value = inf;
best_position = [];
best_curve = nan(1, max_fes);

counted_fobj = @evaluate;
get_state = @state;

    function value = evaluate(position)
        if fe_count >= max_fes
            value = inf;
            return;
        end

        value = base_fobj(position);
        fe_count = fe_count + 1;

        if value < best_value
            best_value = value;
            best_position = position;
        end
        best_curve(fe_count) = best_value;
    end

    function output = state()
        output.fe_count = fe_count;
        output.best_value = best_value;
        output.best_position = best_position;
        output.best_curve = best_curve(1:fe_count);
    end
end
