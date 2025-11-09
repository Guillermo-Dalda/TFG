% A straightforward implementation of the NIZAR Optimization Algorithm (NOA)
% https://link.springer.com/article/10.1007/s11227-023-05579-4
% N.C. Cruz, University of Granada, Spain, 2025

% INPUTS:
% func: Function handle of of the objective function (to minimize)
% bounds: nDim x 2 matrix with var. bounds. 1st Col.: lower bnds, 2nd uppr.
% popSize: Population size. A population is a nDim x popSize matrix
% maxCycles: Maximum number of iterations

function [sol,fval] = NIZAROptimizer(func, bounds, popSize, maxCycles)
    [pop, vals] = CreatePopulation(popSize, bounds, func); % A row per dim. A column per individual
    for c = 1:1:maxCycles
        % Step 4:
        [xBest, ~] = TakeBest(pop, vals);
        alfas = [rand(), (rand() - 3/8), (rand() - 1/4), (rand() - 3/8)];
        lambdas = randi([0 1], [1, 8]); % Row vec. of 8 rand. ints 0 or 1. Equivalent to paper's mod(a_n, 2)
        for i=1:1:popSize
            xi = pop(:, i);
            candidates = [1:i-1, i+1:popSize]; % Every individual apart from 'i' can be chosen
            jkm = randperm(popSize-1, 3); % Random values between [1, popSize-1]
            [xj, xk, xm] = deal( pop(:, candidates(jkm(1))), pop(:, candidates(jkm(2))), pop(:, candidates(jkm(3))));
            [B1, B2] = deal(rand(), rand()); % Two random decimal values in [0, 1]
            p1 = ConstructP1(lambdas, xm, xBest, B1); % COMMONLY REQUIRED VALUE: <PRELIMINARIES>
            % Diversity phase
            if lambdas(1) == 1 % Eq. (2.9) <Step 5>
                xi_prime = BuildS1(lambdas, p1, xi, xj, xk, B1, B2); % <Step 6><Eq. (2.10)> Modified candidate solution
            else % Otherwise, construct the effective points using Eqs. (2.14) and (2.15) <Step 5>
                t3 = ConstructT3(lambdas, xBest, xi); % COMMONLY REQUIRED VALUES: <PRELIMINARIES>
                t2 = ConstructT2(lambdas, alfas, t3, p1, B2); % ORDER MATTERS: T2 USES T3, SO T3 FIRST
                p2 = ConstructP2(lambdas, xm, t2, B1);
                p3 = ConstructP3(lambdas, alfas, p1, t2, t3, B2);
                vj = fi_3(alfas(1), xj, t3); % Eq. (2.12)
                vk = fi_3(alfas(1), xk, t3); % Eq. (2.12)
                dj = B1*((-1)^(candidates(jkm(1)))); % dj = B1 * (-1)^j ; dk = B2 * (-1)^k <Between Eq. (2.11) & Eq. (2.12)>
                dk = B2*((-1)^(candidates(jkm(2))));
                xi_prime = BuildS2(lambdas, dj, dk, p2, p3, vj, vk); % <Step 6><Eq. (2.11)>
            end
            % Overlap phase
            if (all( abs(xi_prime - xBest) < eps) ) || (all( abs(xi - xBest) < eps) ) || rand() <= 1/4 % <Step 7><Eq. (2.19)> % Comparing equality of floating point vectors with an epsilon
                xi_prime = RebuildS1(lambdas, xi_prime, xj, xk, xm, xBest);
            end
            xi_prime = FixLimits(xi_prime, xi, bounds); % <Step 8: Call to Eq. (2.21)>
            val_prime = func(xi_prime);
            if val_prime < vals(i) % Assuming minimization: The modified solution is better-> Replace
                pop(:, i) = xi_prime;
                vals(i) = val_prime;
            end % Discarding the modified version otherwise
        end
    end
    [sol, focusBest] = TakeBest(pop, vals);
    fval = vals(focusBest);
end

% Internal auxiliary functions:

function [pop, vals] = CreatePopulation(popSize, bounds, func)
    nDim = size(bounds, 1); % As many dims. as rows in bounds
    pop = rand(nDim, popSize); % Preliminary 0-1 matrix
    pop = (pop .* (bounds(:, 2) - bounds(:, 1))) + bounds(:, 1); % ind = [0, 1]*(max - min) + min
    vals = zeros(1, popSize); % A value per individual
    for i=1:1:popSize
        vals(i) = func(pop(:, i)); % An individual per column
    end
end

function [xBest, focusBest] = TakeBest(pop, vals)
    [~, focusBest] = min(vals);
    xBest = pop(:, focusBest); % Selecting the best individual (column)
end

function T = fi_1(alfa, x, r) % Eq. (2.2). Alfa is a threshold, x is a solution vector, r a random scalar
    if alfa <= 0.5
        T = x;
    else
        T = round(x + (r^2*(ones(length(x), 1)))); % As many rows in <ones> as rows in X, single column vector
    end
end

function T = fi_2(alfa, x, r) % Eq. (2.3). Alfa is a threshold, x is a solution vector, r a random scalar
    if alfa <= 0.5
        T = x;
    else
        T = round(x * (r^2)); % Scaling
    end
end

function T = fi_3(alfa, xi, xj) % Eq. (2.4).  Alfa is a threshold, xi is a solution vector and xj is another solution vector
    if alfa <= 0.5
        T = xi;
    else
        T = xj;
    end
end

function output = replace(input, target) % See Fig. 2 to 4
    nDim = length(target);
    selMask = (randi([0, 1], [1, nDim]) > 0); % Same probability of being replaced (directly turning double to logic)
    output = input;
    output(selMask) = target(selMask);
end

function output = scramble(input, target) % See Fig. 2 to 4 % WARNING: This function may violate bounds when variables have different limits
    nDim = length(target);
    selMask = (randi([0, 1], [1, nDim]) > 0); % Same probability of being replaced (directly turning double to logic)
    output = input;
    chosen = target(selMask);
    selMask = randperm(nDim, length(chosen)); % Define random injection points
    output(selMask) = chosen;
end

function output = distribute(input, target) % See Fig. 2 to 4 % WARNING: This function may violate bounds when variables have different limits
    nDim = length(target);
    chosen = target(randi(nDim));
    output = input;
    selMask = (randi([0, 1], [1, nDim]) > 0); % Same probability of being replaced (directly turning double to logic)
    output(selMask) = chosen;
end

function p1 = ConstructP1(lambdas, xm, xBest, B1)
    if lambdas(3) == 1 % Eq. (2.13)
        p1 = xm;
    else
        if lambdas(6) == 1 % Eq. (2.16) % T1 will only become p1, so p1
            p1 = 0.5*(xBest + xm); % Column vector expected
        else
            p1 = xm*B1 + (1-B1)*xBest; % Column vector expected
        end
    end
end

function p2 = ConstructP2(lambdas, xm, t2, B1)
    if lambdas(4) == 1 % Eq. (2.14)
        p2 = xm + (B1*ones(length(xm), 1)); % Using random scalr B1 to get a column vector with ones(D, 1), i.e., col. vec of ones.
    else
        p2 = t2;
    end
end

function p3 = ConstructP3(lambdas, alfas, p1, t2, t3, B2)
    if lambdas(5) == 1 % Eq. (2.15)
        p3 = fi_2( alfas(2), fi_3( alfas(3), distribute(t3, p1), t2 ) , B2 );
    else
        p3 = fi_2( alfas(2), fi_3( alfas(3), distribute(p1, t3), t2 ) , B2 );
    end
end

function t2 = ConstructT2(lambdas, alfas, t3, p1, B2)
    if lambdas(7) == 1 % Eq. (2.17)
        t2 = fi_1(alfas(4), replace(t3, p1) , B2);
    else
        t2 = fi_1(alfas(4), scramble(t3, p1) , B2);
    end
end

function t3 = ConstructT3(lambdas, xBest, xi)
    if lambdas(8) == 1 % Eq. (2.18)
        t3 = xBest;
    else
        t3 = xi;
    end
end

function s1 = BuildS1(lambdas, p1, xi, xj, xk, B1, B2)
    if lambdas(2) == 1 % Eq. (2.10)
        s1 = p1 + B1*(xi - xj) + B2*(xi - xk);
    else
        s1 = xi + B1*(p1 - xj) - B2*(p1 - xk);
    end
end

function s2 = BuildS2(lambdas, dj, dk, p2, p3, vj, vk)
    if lambdas(2) == 1 % Eq. (2.11)
        s2 = p2 + dj*(p3 - vj) + dk*(p3 - vk);
    else
        s2 = p3 + dj*(p2 - vj) - dk*(p2 - vk);
    end
end

function s1 = RebuildS1(lambdas, xi_prime, xj, xk, xm, xBest) % reupdate the new solution
    [B1, B2] = deal( rand(length(xj), 1)*2 - 1, rand(length(xj), 1)*2 - 1 ); % Explanation below Eq. (2.19) in the paper <RANGE BETWEEN -1 AND 1>
    if lambdas(3) == 1 % Refactor of p1 with a new kind of Beta during the overlap stage
        p1 = xm;
    else
        if lambdas(6) == 1 % <Injected Eq. (2.16)>
            p1 = 0.5 * (xBest + xm);
        else
            p1 = xm .* B1 + (ones(length(xj), 1) - B1) .* xBest;
        end
    end
    if lambdas(2) == 1 % Modified Eq. (10) to reupdate the new solution during the overlap phase
        s1 = p1 + B1 .* (xi_prime - xj) + B2 .* (xi_prime - xk);
    else
        s1 = xi_prime + B1 .* (p1 - xj) -  B2 .* (p1 - xk);
    end
end

function xi_prime = FixLimits(xi_prime, xi, bounds) % nDim x 2 matrix with var. bounds. 1st Col.: lower bnds, 2nd uppr.
    tooBig = xi_prime > bounds(:, 2);
    tooSmall = xi_prime < bounds(:, 1);
    xi_prime(tooBig) = xi(tooBig);
    xi_prime(tooSmall) = xi(tooSmall);
end
