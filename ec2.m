clear;
transition = [0.5, 0.5; 0.5, 0.5];
emission = [0.4, 0.1, 0.5; 0.1, 0.5, 0.4];
LA = 1;
NY = 2;
NULL = 3;
O = [NULL, LA, LA, NULL, NY, NULL, NY, NY, NY, NULL, NY, NY, NY, NY,...
    NY, NULL, NULL, LA, LA, NY];
pi = [0.5, 0.5];
iteration = 0;
while true
    iteration = iteration + 1;
    delta = zeros(2, size(O, 2));
    beta = zeros(2, size(O, 2));
    for i = 1 : 2
        delta(i, 1) = pi(i) * emission(i, NULL);
        beta(i, size(O, 2)) = 1;
    end 
    for t = 2 : size(O, 2)
        for i = 1 : 2
            delta(i, t) = emission(i, O(t)) * (delta(:, t-1)' * ...
                            transition(:, i));
        end
    end
    for t = size(O, 2)-1 : -1 : 1
        for i = 1 : 2
            beta(i, t) = beta(LA, t+1) * transition(i, LA) * emission(LA, O(t+1))...
                        + beta(NY, t+1) * transition(i, NY) * emission(NY, O(t+1));
        end 
    end
    gamma = zeros(2, size(O, 2));
    for t = 1 : size(O, 2)
        div = delta(:, t)' * beta(:, t);
        gamma(:, t) = delta(:, t) .* beta(:, t) ./ div;
    end
    epsilon = zeros(2, 2, size(O, 2)-1);
    for t = 1 : size(O, 2)-1
        div = delta(:, t)' * beta(:, t);
        for i = 1 : 2
            for j = 1 : 2
                epsilon(i, j, t) = delta(i, t) * transition(i, j) * beta(j, t+1) *...
                                    emission(j, O(t+1)) / div;
            end 
        end 
    end
    for i = 1 : 2
        pi(i) = gamma(i, 1);
        for j = 1 : 2
            transition(i, j) = sum(epsilon(i, j, :)) / sum(gamma(i, 1:end-1));
        end
        for j = 1 : 3
            gSum = 0;
            for t = 1 : size(O, 2)
                if O(t) == j
                    gSum = gSum + gamma(i, t);
                end
            end
            emission(i, j) = gSum / sum(gamma(i, :));
        end
    end
    P(iteration) = sum(delta(:, end));
    if (iteration > 1) && abs(P(iteration) - P(iteration-1)) < eps
        break;
    end
end
figure;
plot(P(1, :));
xlabel('Iteration');
ylabel('Evaluation');
display(pi);
display(transition);
display(emission);