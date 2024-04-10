function P = ipot_WD(a1, a2, C, beta)
    % Solve the optimal transport problem and return the OT matrix
    % The function solves the following optimization problem:
    %   gamma = argmin_gamma <gamma, C>_F
    %   s.t. gamma * 1 = a1
    %        gamma' * 1 = a2
    %        gamma >= 0
    % where:
    %   - C is the (ns, nt) metric cost matrix
    %   - a1 and a2 are source and target weights (sum to 1)
    % The algorithm uses proximal point method
    max_iter=2000;
    % L=1;
    % use_path = 1;
    % return_map = 1;
    % return_loss = 0;
    



    n = length(a1);
    v = ones(n, 1);
    u = ones(n, 1);
    P = ones(n,n) / n^2;

    K = exp(-C / beta);
    % if return_loss
    %     loss = zeros(1, max_iter);
    % end

    for outer_i = 1:max_iter
        Q = K .* P;

        % if ~use_path
        %     v = ones(n, 1);
        %     u = ones(n, 1);
        % end
            u = a1 ./ bsxfun(@times, Q ,v);
            v = a2 ./ bsxfun(@times, Q' ,v);
        PP = bsxfun(@times, Q ,v);
        P=bsxfun(@times, PP ,u);
        % if return_loss
        %     W = sum(P1(:) .* C(:));
        %     loss(outer_i) = W;
        % end
    end

    % if return_loss
    %     if return_map
    %         P = P1;
    %     else
    %         P = loss;
    %     end
    % else
    %     if return_map
    %         P = P1;
    %     else
    %         P = [];
    %     end
    % end
end

% function bar = geometricBar(weights, alldistribT)
%     % Return the weighted geometric mean of distributions
%     assert(length(weights) == size(alldistribT, 2));
%     bar = exp(log(alldistribT) * weights');
% end
% 
% function mean_dist = geometricMean(alldistribT)
%     % Return the geometric mean of distributions
%     mean_dist = exp(mean(log(alldistribT), 2));
% end
% 
% function barycenter = ipot_barycenter(A, M, beta, weights, numItermax)
%     % Compute the Wasserstein barycenter of distributions A
%     % The function solves the following optimization problem:
%     %   a = argmin_a sum_i W(a, a_i)
%     % where:
%     %   - W(a, a_i) is the Wasserstein distance
%     %   - a_i are training distributions in the columns of matrix A
%     %   - M is the cost matrix for OT
%     % The algorithm absorbs many tricks in Python package "POT"
% 
%     if nargin < 5
%         numItermax = 1000;
%     end
% 
%     if nargin < 4
%         weights = ones(1, size(A, 2)) / size(A, 2);
%     else
%         assert(length(weights) == size(A, 2));
%     end
% 
%     [n, k] = size(A);
%     cpt = 0;
% 
%     K = exp(-M / beta);
%     Pi = ones(n, n, k) / (n * n);
%     K = reshape(K, [n, n, 1]);
%     Q = K .* Pi;
% 
%     v = A ./ sum(Q, 1);
%     UKv = sum(Q .* reshape(v, [1, n, k]), 2);
%     u = geometricMean(UKv) ./ UKv';
% 
%     while cpt < numItermax
%         cpt = cpt + 1;
% 
%         Q = K .* Pi;
% 
%         for i = 1
%             u = (u' .* geometricBar(weights, UKv))' ./ UKv; % for numerical stability
%             v = A ./ sum(Q .* reshape(u, [n, 1, k]), 1);
%             Pi = reshape(u, [n, 1, k]) .* Q .* reshape(v, [1, n, k]);
%             Pi = Pi ./ sum(sum(Pi)); % for numerical stability
%             UKv = u' .* sum(Q .* reshape(v, [1, n, k]), 2);
%         end
%     end
% 
%     barycenter = geometricBar(weights, UKv);
% end
