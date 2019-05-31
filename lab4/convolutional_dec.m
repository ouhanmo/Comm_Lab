function decoded_data = convolutional_dec(binary_data, impulse_response)

size_response = size(impulse_response);
num_channels = size_response(1);
memory_length = size_response(2)-1;

size_raw = length(binary_data)/num_channels;

decoded_data = zeros(1, size_raw);

path = zeros(2^memory_length , size_raw);
distance = zeros(2^memory_length, size_raw +1);
distance(2:end, 1) = inf;

for transition_step = 1:size_raw
    received_bits = binary_data((transition_step-1)*num_channels+1 : transition_step*num_channels);
    
    for case_num = 0:2^memory_length-1
        memory_bits = de2bi(case_num,memory_length);
        predict_bits = zeros(2, num_channels);
        for channel_id = 1:num_channels
            predict_bits(1,channel_id) = mod(dot([memory_bits 0], impulse_response(channel_id, 1:end)),2);
            predict_bits(2,channel_id) = mod(dot([memory_bits 1], impulse_response(channel_id, 1:end)),2);
        end
        memory0 = bi2de([memory_bits(2:end) 0]);
        memory1 = bi2de([memory_bits(2:end) 1]);
        distance0 = distance(memory0 +1, transition_step) + sum(abs(received_bits-predict_bits(1,1:end)));
        distance1 = distance(memory1 +1, transition_step) + sum(abs(received_bits-predict_bits(2,1:end)));
        
        if distance0 <= distance1
            path(case_num + 1, transition_step) = memory0 +1 ;
            distance(case_num +1 ,transition_step +1) = distance0;
        else
            path(case_num + 1, transition_step) = memory1 +1 ;
            distance(case_num +1 ,transition_step +1) = distance1;
        end
    end
        
end

% Find Shortest Distance

[min_distance, min_position] = min(distance(1:end, end));


for ii = size_raw:-1:1
    last_memory = de2bi(min_position-1, memory_length);
    last_bit = last_memory(1);
    decoded_data(ii) = last_bit;
    min_position = path(min_position, ii);
end

if size_raw <100
    
distance
path
end



