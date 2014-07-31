function pred_time = predict_times_PCA(train_data, train_times, test_data, nmodes)

[m, n] = size(train_data);
[m2, n2] = size(test_data);

mean_train_data = mean(train_data);

train_data = train_data - repmat(mean_train_data, m, 1);
test_data = test_data - repmat(mean_train_data, m2, 1);

[U, S, V] = svds(train_data, nmodes);

b = regress(train_times, [train_data*V ones(m, 1)]);

pred_time = [test_data*V ones(m2, 1)]*b;


