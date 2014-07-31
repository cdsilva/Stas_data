function shift_times = estimate_shift_times(PCA_data, time, nmodes)

nmovies = length(PCA_data);
bias = zeros(nmovies);
variance = zeros(nmovies);

%%
for train_movie=1:nmovies
        
    for i=1:nmovies
        pred_time = predict_times_PCA(PCA_data{train_movie}, time{train_movie}, PCA_data{i}, nmodes);
        
        bias(train_movie, i) = mean(pred_time - time{i});
        variance(train_movie, i) = mean((pred_time - time{i}).^2);        
    end
end

%% adjust movies

synch_matrix = (bias-bias')/2;
synch_scale = 4*(max(synch_matrix(:))-min(synch_matrix(:)));
synch_matrix = exp(sqrt(-1) * synch_matrix * (2*pi/synch_scale));

[V, ~] = eigs(synch_matrix, 1);

shift_times = atan2(imag(V), real(V)) * synch_scale/(2*pi);
shift_times = shift_times - max(shift_times);