%% Part 2.1
p = default_params();


% Test points, each row is one case
X  = [0.5 12 4   20;
      1.0  8 2   15;
      0.4 18 5.5 22];
TH = [2.0 0.0325 648000 4;
      2.3 0.035  612000 6;
      1.7 0.028  684000 2];

R_test = zeros(3,1);
g_test = zeros(3,4);
for k = 1:3
    [R_test(k), g_test(k,:)] = uav_range(X(k,:), TH(k,:), p);
end

% With help of Claude
fprintf('Task 2.1: UAV model test points\n');
fprintf('%4s %14s %14s %14s %14s %14s\n', 'Case', 'R [m]', 'g1', 'g2', 'g3', 'g4');
for k = 1:3
    fprintf('%4d %14.6g %14.6g %14.6g %14.6g %14.6g\n', k, R_test(k), g_test(k,:));
end

% Vectorization check, same theta in 5 rows should give 5 identical rows
[R_vec, g_vec] = uav_range(X(1,:), repmat(TH(1,:), 5, 1), p);
assert(isequal(size(R_vec), [5 1]) && isequal(size(g_vec), [5 4]), 'Wrong output size')
assert(all(R_vec == R_vec(1)), 'Rows of R do not match')