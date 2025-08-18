function [Ceq1_fit, Ceq1, Convergence_curve] = PFO(Particles_no, Max_iter, lb, ub, dim, fobj)
%-----------------------------------------------------------------------------------------------------------
% PFO: Plane Force System Optimizer 
%-----------------------------------------------------------------------------------------------------------

% Initialization
Ceq1 = zeros(1, dim); Ceq1_fit = inf; 
Ceq2 = zeros(1, dim); Ceq2_fit = inf; 

C = initialization(Particles_no, dim, ub, lb); 
Convergence_curve = zeros(1, Max_iter);
fitness = zeros(1, Particles_no);

% Initial fitness calculation
for i = 1:Particles_no
    C(i,:) = bound_reflect(C(i,:), lb, ub); % Rebound boundary 
    fit_val = fobj(C(i,:));
    % Stability check
    if isnan(fit_val) || isinf(fit_val)
        fprintf('Warning: Initial fitness NaN or Inf at particle %d, set to large number\n', i);
        fit_val = 1e10;
    end
    fitness(i) = fit_val;
    if fitness(i) < Ceq1_fit
        Ceq2_fit = Ceq1_fit; Ceq2 = Ceq1;
        Ceq1_fit = fitness(i); Ceq1 = C(i,:);
    elseif fitness(i) > Ceq1_fit && fitness(i) < Ceq2_fit
        Ceq2_fit = fitness(i); Ceq2 = C(i,:);
    end
end

for t = 1:Max_iter

    for i = 1:Particles_no
        P = randi([1,2]); % Branch selection

        old_pos = C(i,:);
        old_fit = fitness(i);

        if P == 1 % ------------------ Branch 1 ------------------
            r_branch1 = rand;

            if r_branch1 > 0.8 % Plane Concurrent Force System
                epsilon = 0.1; omega = 0.75;
                mu = pi * (Max_iter - t) / Max_iter;
                delta = omega * (1 + t / Max_iter);
                eta = epsilon + delta * (1 - epsilon) * mu;

                rand_indices = randperm(Particles_no, 2);
                Fj = C(rand_indices(1), :);
                Fk = C(rand_indices(2), :);

                S = eta * (Fj - Fk);

                r_m1 = rand(1, dim);
                l_m1 = rand(1, dim);
                m1 = 2 * r_m1 .* (exp(-l_m1 * t / Max_iter) - 1);

                new_pos = old_pos + (m1 .* (Ceq1 - old_pos) + S);

            elseif r_branch1 > 0.2 % Plane General Force System
                if rand > 0.5
                    u = 0.2; v = 0.9;
                    w = v - t * ((v - u) / Max_iter);
                    j_idx = randi(Particles_no);
                    Fj = C(j_idx, :);
                    theta = 2 * pi * rand;
                    phi = w * (Fj - old_pos) * exp(theta);
                    r1 = rand; r2 = rand;
                    new_pos = old_pos + phi * r1 .* sin(2 * pi * r2);
                else
                    Fm = mean(C, 1);
                    sigma_rot = cos((pi / 2) * (t / Max_iter)^2) .* (Ceq1 - Fm);
                    gamma_rot = (old_pos + 1e-6*randn(1,dim)) .* sigma_rot;
                    new_pos = gamma_rot;
                end
            else
                new_pos = old_pos; % r_branch1 <= 0.2
            end

        else % P == 2
            if rand > 0.5 % Plane Force Couple System
                alpha = 1 - (t / Max_iter);
                beta_couple = (Ceq1 - old_pos) * alpha;
                v_couple = rand;
                new_pos = v_couple * Ceq1 + (1 - v_couple) * beta_couple;
            else % Plane Parallel Force System
                h = 0.5;
                r_par = rand;
                theta_par = 2 * pi * rand;
                levy_term = (r_par * cos(theta_par)) .* Levy(dim) .* r_par .* (h * (Ceq2 - old_pos));
                new_pos = Ceq1 + levy_term + h * (Ceq1 - old_pos);
            end
        end

        % Rebound boundary
        new_pos = bound_reflect(new_pos, lb, ub);

        % Calculate the new fitness value
        new_fit = fobj(new_pos);
        % Stability check
        if isnan(new_fit) || isinf(new_fit)
            fprintf('Warning: Fitness NaN or Inf at iter %d, particle %d, set to large number\n', t, i);
            new_fit = 1e10;
        end

        if new_fit < old_fit
            C(i,:) = new_pos;
            fitness(i) = new_fit;

            % Update the global optimum
            if new_fit < Ceq1_fit
                Ceq2_fit = Ceq1_fit; Ceq2 = Ceq1;
                Ceq1_fit = new_fit; Ceq1 = new_pos;
            elseif new_fit > Ceq1_fit && new_fit < Ceq2_fit
                Ceq2_fit = new_fit; Ceq2 = new_pos;
            end
        end
 
    end

    % Record the convergence values and ensure there are no NaN/Inf values.
    if isnan(Ceq1_fit) || isinf(Ceq1_fit)
        fprintf('Warning: Ceq1_fit NaN or Inf at iter %d, set to large number\n', t);
        Ceq1_fit = 1e10;
    end
    Convergence_curve(t) = Ceq1_fit;
end
end

%% Initialization

function C = initialization(Particles_no, dim, ub, lb)
    Boundary_no = numel(ub);
    if Boundary_no == 1
        C = rand(Particles_no, dim) .* (ub - lb) + lb;
    else
        C = zeros(Particles_no, dim);
        for i = 1:dim
            C(:,i) = rand(Particles_no, 1) .* (ub(i) - lb(i)) + lb(i);
        end
    end
end

%% Rebound boundary (scalar/vector compatible)
function pos = bound_reflect(pos, lb, ub)
    dim = numel(pos);
    if numel(lb) == 1
        lb = repmat(lb, 1, dim);
    end
    if numel(ub) == 1
        ub = repmat(ub, 1, dim);
    end
    for d = 1:dim
        if pos(d) < lb(d)
            pos(d) = lb(d) + (lb(d) - pos(d));
            if pos(d) > ub(d)
                pos(d) = ub(d);
            end
        elseif pos(d) > ub(d)
            pos(d) = ub(d) - (pos(d) - ub(d));
            if pos(d) < lb(d)
                pos(d) = lb(d);
            end
        end
    end
end

%% Levy Distribution
function o = Levy(d)
    beta = 1.5;
    sigma_num = gamma(1 + beta) * sin(pi * beta / 2);
    sigma_den = gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2);
    sigma = (sigma_num / sigma_den)^(1 / beta);

    u = randn(1, d) * sigma;
    v = randn(1, d);
    step = u ./ abs(v).^(1 / beta);
    o = step;
end
