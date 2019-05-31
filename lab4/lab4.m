close all
clear 

exper = 2;


if exper == 1
    num_bits = 1000;
    R = 1000;
    d = 1;
    power_avg = 10*d*d; 
    snrs = 0:3:15;
    num_bit_error = zeros(1, length(snrs));
    num_sym_error = zeros(1, length(snrs));
    figure;
    for ii = 1:length(snrs)
        snr = snrs(ii);
        for iteration = 1:R
            data_bits = randi([0 1], 1,num_bits);
            data_symbols = symbol_mapper(data_bits, 16, 1, "QAM", "Grey");
             % scatter(real(data_symbols), imag(data_symbols));
            noise = sqrt(power_avg/2)*(randn(1,length(data_symbols)) + j*(randn(1,length(data_symbols))));
            channel_out = data_symbols + 10^(-snr/20)*noise;
            received_bits = symbol_demapper(channel_out, 16, 1, "QAM", "Grey", "MD");
            num_bit_error(ii) = num_bit_error(ii) +  sum(abs(received_bits-data_bits));
            for symbol_it = 1:num_bits/4
                if sum(abs(received_bits((symbol_it-1)*4+1:(symbol_it)*4) - data_bits((symbol_it-1)*4+1:(symbol_it)*4))) > 0
                    num_sym_error(ii) = num_sym_error(ii) +1 ;
                end
            end
        end
    end
    
    figure;
    semilogy(snrs, num_bit_error/(num_bits*R), "bo-" ,"Linewidth",2);
    hold on
    semilogy(snrs, num_sym_error/(num_bits*R/4), "go-" ,"Linewidth",2);
    legend("BER","SER");
    xlabel("SNR(dB)")
    ylabel("Error Rates")
    title("BER vs SER R=" + R)
    figure;
    % semilogy(snrs, num_bit_error/(num_bits*R), "bo-" ,"Linewidth",2);
  
    semilogy(snrs, num_sym_error/(num_bits*R/4), "go-" ,"Linewidth",2);
    hold on
    grid on
    theory = zeros(1,length(snrs));
    upper  = zeros(1,length(snrs));
    lower  = zeros(1,length(snrs));
    
    for ii = 1:length(snrs)
        Eb = (3/4)*erfc(1/2*sqrt((3)/(2*15*10^(-snrs(ii)/10))));
        theory(ii) =  2*Eb-Eb^2;
        upper(ii)  = 2*Eb;
        lower(ii)  = 2/3*Eb;
    end
    
    semilogy(snrs,theory, "kx-", "Linewidth", 2);
    semilogy(snrs,upper, "rx-", "Linewidth", 2);
    semilogy(snrs,lower, "yx-", "Linewidth", 2);
    legend("Simulation","Theory", "Upper Bound", "Lower Bound");
    xlabel("SNR(dB)")
    ylabel("Error Rates")
    title("SERs")
end

if exper == 2
    data_bits = randi([0 1] , 1, 1000);
    impulse_response = [ 1 1 1 1; 1 1 0 1];

    encoded_bits = convolutional_enc(data_bits, impulse_response);
    decoded_bits = convolutional_dec(encoded_bits, impulse_response);
    num_errors = sum(abs(decoded_bits - data_bits))
end

if exper == 3
    [x,fs] = audioread("handel.ogg");
    x = x';
    x_quan = quantizer_L_level(x,max(x), 16);
    symbols = unique(x_quan);
    symbol_count = zeros(1,16);
    for ii = 1:length(x_quan)
        j = find(symbols==x_quan(ii));
        symbol_count(j) = symbol_count(j) + 1;
    end
    p = symbol_count/length(x_quan);
    dict = huffman_dict(symbols,p);
    x_sc = huffman_enc(x_quan(1:100000), dict);
    x_cc = convolutional_enc(x_sc, [ 1 1 1 1; 1 1 0 1]);
    x_symbols = symbol_mapper(x_cc, 4, 1, "PSK" , "Grey");
    
    snr = 10;
    
    y_symbols = x_symbols + 0.5*10^(-snr/20)*(randn(size(x_symbols)) + j*randn(size(x_symbols)));
    y_cc = symbol_demapper(y_symbols, 4, 1, "PSK" , "Grey", "MD");
    y_sc = convolutional_dec(y_cc,  [ 1 1 1 1; 1 1 0 1]);
    y_quan = huffman_dec(y_sc , dict);
    
    sound(y_quan , fs);
    
end


