% mat1 = ones(3, 3);
% % mat2 = ones(3, 3);
% mat2 = zeros(3, 3);
mat1 = [ones(5, 2), zeros(5, 3)];
% mat2 = mat1;
mat2 = [zeros(5, 2), ones(5, 3)];
% figure; imagesc(mat1)

h = 1;

testw = weightGWED(mat1, mat2, h);