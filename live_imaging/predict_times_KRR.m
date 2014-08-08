function pred_time = predict_times_KRR(train_data, train_times, test_data, nmodes, lambda)

[m, n] = size(train_data);
[m2, n2] = size(test_data);

mean_train_data = mean(train_data);

train_data = train_data - repmat(mean_train_data, m, 1);
test_data = test_data - repmat(mean_train_data, m2, 1);

[U, S, V] = svds(train_data, nmodes);

xtrain = train_data * V;
xtest = test_data * V;

K = squareform(pdist(xtrain)).^2;
eps = median(K(:));
K = exp(-K/eps);

k = pdist2(xtrain, xtest).^2;
k = exp(-k/eps);

I = eye(m);

pred_time = (train_times' * inv(K + lambda * I) * k)';
