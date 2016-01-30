clear;
data = csvread('tshirts-I.csv', 1, 1);
[m, ~] = size(data);
rng(0);
randIndex = randperm(size(data, 1));
% Customize number of clusters
K = 4;
% Initial center points
u = data(randIndex(1:K), :);
clusteredValues = [data ones(m, 1)];
% First figure required
figure;
hold on;
plot(u(:, 1), u(:, 2), 'bo');
plot(data(:, 1), data(:, 2), 'rx');
xlabel('weight');
ylabel('height');
title('Initial Seeds');
% Generate visible colors for K values
colorR = zeros(1, K);
colorG = zeros(1, K);
colorB = zeros(1, K);
for i = 1 : K
    colorVal = 2.1*i / K;
    if colorVal >= 0.7
        colorR(i) = 0.7;
    elseif colorVal > 0
        colorR(i) = colorVal;
    end
    colorVal = colorVal - 0.7;
    if colorVal >= 0.7
        colorG(i) = 0.7;
    elseif colorVal > 0
        colorG(i) = colorVal;
    end
    colorVal = colorVal - 0.7;
    if colorVal >= 0.7
        colorB(i) = 0.7;
    elseif colorVal > 0
        colorB(i) = colorVal;
    end
end
% Start clustering
iteration = 1;
while true
    preU = u;
    uCounts = ones(K, 1);
    for i = 1 : m
        % Calculate closest cluster center
        minIndex = 1;
        minDistance = norm(data(i, :) - u(1, :));
        for k = 2 : K
            distance = norm(data(i, :) - u(k, :));
            if distance < minDistance 
                minIndex = k;
                minDistance = distance;
            end
        end
        % Update cluster center
        clusteredValues(i, end) = minIndex;
        u(minIndex, :) = (u(minIndex, :) .* uCounts(minIndex, 1)+...
                                data(i, :)) ./ (uCounts(minIndex, 1)+1);
        uCounts(minIndex, 1) = uCounts(minIndex, 1) + 1;
    end
    if iteration == 1
        % Second figure required
        figure;
        hold on;
        for i = 1 : K
            plot(u(i, 1), u(i, 2), 'o', 'Color', [colorR(i), colorG(i), colorB(i)]);
            cluster = clusteredValues(:, end) == i;
            plot(clusteredValues(cluster, 1), clusteredValues(cluster, 2), 'x',...
                'Color', [colorR(i), colorG(i), colorB(i)]);
        end
        xlabel('weight');
        ylabel('height');
        title('Iteration 1');
    end
    iteration = iteration + 1;
    if norm(preU - u) < eps
        break;
    end
end
% Last figure after clustering
figure;
hold on;
for i = 1 : K
    plot(u(i, 1), u(i, 2), 'o', 'Color', [colorR(i), colorG(i), colorB(i)]);
    cluster = clusteredValues(:, end) == i;
    plot(clusteredValues(cluster, 1), clusteredValues(cluster, 2), 'x',...
        'Color', [colorR(i), colorG(i), colorB(i)]);
end
xlabel('weight');
ylabel('height');
title(['Iteration ' num2str(iteration)]);