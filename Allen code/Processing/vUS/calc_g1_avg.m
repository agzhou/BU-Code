load([savepath, 'fUSdata-', num2str(1), '.mat'], 'g1')

g1_sum = cell([size(g1)]);
for j = 1:3
    g1_sum{j} = zeros(size(g1{j}));
end

% for filenum = startFile + 1:endFile
for filenum = startFile:endFile
    tic
    load([savepath, 'fUSdata-', num2str(filenum), '.mat'], 'g1')


    for j = 1:3
        g1_sum{j} = g1_sum{j} + g1{j};
    end
    toc
end

g1_avg = cell(size(g1_sum));
for j = 1:3
    g1_avg{j} = g1_sum{j}./(endFile - startFile + 1);
end

save([savepath, 'g1_avg.mat'], 'g1_avg', '-v7.3')