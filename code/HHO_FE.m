function [Best_fit, Best_pos, Convergence_curve] = HHO_FE( ...
    Npop, MaxFEs, lb, ub, dim, fobj)
%HHO_FE Harris Hawks Optimization with an FE-based stopping condition.
%
% This implementation preserves the HHO update rules while expressing the
% decreasing escape-energy schedule as a fraction of the consumed objective
% evaluations. It is used only for equal-FE benchmarking because standard HHO
% iterations have a branch-dependent number of objective evaluations.

if isscalar(lb), lb = lb * ones(1,dim); end
if isscalar(ub), ub = ub * ones(1,dim); end

X = rand(Npop,dim).*(ub-lb)+lb;
fit = inf(Npop,1);
Best_fit = inf;
Best_pos = zeros(1,dim);
fe_count = 0;
Convergence_curve = [];

for i = 1:Npop
    fit(i) = fobj(X(i,:));
    fe_count = fe_count+1;
    if fit(i) < Best_fit
        Best_fit = fit(i);
        Best_pos = X(i,:);
    end
end

iteration = 0;
while fe_count < MaxFEs
    iteration = iteration+1;
    progress = fe_count/MaxFEs;
    E1 = 2*(1-progress);
    X_mean = mean(X,1);

    for i = 1:Npop
        if fe_count >= MaxFEs, break; end

        E0 = 2*rand-1;
        Esc = E1*E0;

        if abs(Esc) >= 1
            q = rand;
            X_rand = X(randi(Npop),:);
            if q >= 0.5
                candidate = X_rand-rand(1,dim).* ...
                    abs(X_rand-2*rand(1,dim).*Best_pos);
            else
                candidate = (Best_pos-X_mean)-rand(1,dim).* ...
                    (lb+rand(1,dim).*(ub-lb));
            end
            [X,fit,Best_fit,Best_pos,fe_count] = evaluate_candidate( ...
                candidate,i,X,fit,Best_fit,Best_pos,fe_count, ...
                MaxFEs,lb,ub,fobj);
            continue;
        end

        r = rand;
        delta = Best_pos-X(i,:);
        if r >= 0.5 && abs(Esc) >= 0.5
            candidate = delta-Esc.*abs(delta);
            [X,fit,Best_fit,Best_pos,fe_count] = evaluate_candidate( ...
                candidate,i,X,fit,Best_fit,Best_pos,fe_count, ...
                MaxFEs,lb,ub,fobj);
        elseif r >= 0.5
            candidate = Best_pos-Esc.*abs(delta);
            [X,fit,Best_fit,Best_pos,fe_count] = evaluate_candidate( ...
                candidate,i,X,fit,Best_fit,Best_pos,fe_count, ...
                MaxFEs,lb,ub,fobj);
        else
            Y = Best_pos-Esc.*abs(delta);
            [X,fit,Best_fit,Best_pos,fe_count] = evaluate_candidate( ...
                Y,i,X,fit,Best_fit,Best_pos,fe_count, ...
                MaxFEs,lb,ub,fobj);
            if fe_count < MaxFEs
                Z = Y+rand(1,dim).*levy_step(dim,1.5);
                [X,fit,Best_fit,Best_pos,fe_count] = evaluate_candidate( ...
                    Z,i,X,fit,Best_fit,Best_pos,fe_count, ...
                    MaxFEs,lb,ub,fobj);
            end
        end
    end
    Convergence_curve(iteration) = Best_fit; %#ok<AGROW>
end
end

function [X,fit,Best_fit,Best_pos,fe_count] = evaluate_candidate( ...
    candidate,i,X,fit,Best_fit,Best_pos,fe_count,MaxFEs,lb,ub,fobj)
if fe_count >= MaxFEs, return; end
candidate = max(lb,min(ub,candidate));
value = fobj(candidate);
fe_count = fe_count+1;
if value < fit(i)
    X(i,:) = candidate;
    fit(i) = value;
end
if value < Best_fit
    Best_fit = value;
    Best_pos = candidate;
end
end

function L = levy_step(dim,beta)
num = gamma(1+beta)*sin(pi*beta/2);
den = gamma((1+beta)/2)*beta*2^((beta-1)/2);
sigma = (num/den)^(1/beta);
u = randn(1,dim)*sigma;
v = randn(1,dim);
L = u./abs(v).^(1/beta);
end
