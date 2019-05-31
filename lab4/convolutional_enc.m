function encoded_data = convolutional_enc(binary_data, impulse_response)

size_response = size(impulse_response);
num_channels = size_response(1);
% constraint_length = size_response(2);

encoded_data = zeros(1, num_channels*length(binary_data));

for train_no = 1:num_channels
   data_conv = conv(binary_data, impulse_response(train_no, 1:end));
   data_conv = data_conv(1:length(binary_data));
   encoded_data(train_no: num_channels:end) = mod(data_conv,2);
end
