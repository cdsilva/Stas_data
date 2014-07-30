function pred_time = predict_times_PCA(train_data, train_times, test_data, nmodes)

[m, n] = size(train_data);
[m2, n2] = size(test_data);

mean_train_data = mean(train_data);

train_data = train_data - repmat(mean_train_data, m, 1);
test_data = test_data - repmat(mean_train_data, m2, 1);

[U, S, V] = svds(train_data, m);

d = diag(S).^2;

b = regress(train_times, [train_data*V(:, 1:nmodes) ones(m, 1)]);

pred_time = [test_data*V(:, 1:nmodes) ones(m2, 1)]*b;


