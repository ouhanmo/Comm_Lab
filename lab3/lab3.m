close all
clear
clc

exper = 5;




if exper == 1
    binary_seq = [1 0 0 1 1 0 1 1 1 1 0 1 0 1 0 0 0 1 0 1 0 0 1 1];
    symbols = symbol_mapper(binary_seq, 16 , 1, "QAM", "Binary")
    symbols_grey = symbol_mapper(binary_seq, 16 , 1, "QAM", "Grey")
    
    

end

if exper == 2
    W = 50;
    T = 1/(2*W);
    oversampling_factor = 1000;
    T_os = 1/oversampling_factor; % time spacing after oversampling
    pulse_duration = 1; % 1 sec
    t_axis = (-pulse_duration/2 : T_os : pulse_duration/2 - T_os);
    
    
    
    rect = -sign(min(abs(t_axis),T/2)-T/2);
    figure;
    plot(t_axis,rect);
    hold on
    scatter(-0.5:T:0.5, zeros(size(-0.5:T:0.5)),"x")
    title("Rectangular")
    xlabel("Time(s)")
    ylabel("Magnitude")
    X_r = circshift(T_os*fft(circshift(rect, oversampling_factor/2)), oversampling_factor/2);
    figure;
    plot(-oversampling_factor/2:oversampling_factor/2-1,abs(X_r));
    title("Fourier Transform  - Rectangular")
    xlabel("Frequency(Hz)")
    ylabel("Magnitude")
    
    X_sum_r = zeros(1,2*W);
    for ii = 1:10
        X_sum_r = X_sum_r + X_r((ii-1)*100 +1 : ii*100);
    end 
    figure;
    plot(1:100,abs(X_sum_r));
    title("Nyquist Condition Sum Check - Rectangular")
    xlabel("Frequency(Hz)")
    ylabel("Magnitude")
    
    g = sinc(t_axis/T);
    figure;
    plot(t_axis,g);
     hold on
    scatter(-0.5:T:0.5, zeros(size(-0.5:T:0.5)),"x")
    title("Sinc")
    xlabel("Time(s)")
    ylabel("Magnitude")
    
    X_s = circshift(T_os*fft(circshift(g, oversampling_factor/2)), oversampling_factor/2);
    figure;
    plot(-oversampling_factor/2:oversampling_factor/2-1,abs(X_s));
    title("Fourier Transform  - Sinc")
    xlabel("Frequency(Hz)")
    ylabel("Magnitude")
    
    X_sum_s = zeros(1,2*W);
    for ii = 1:10
        X_sum_s = X_sum_s + X_s((ii-1)*100 +1 : ii*100);
    end 
    figure;
    plot(1:100,abs(X_sum_s));
    title("Nyquist Condition Sum Check - Sinc")
    xlabel("Frequency(Hz)")
    ylabel("Magnitude")
    
    
    beta = 0.5;
    g_b = zeros(1,length(t_axis));
    for ii = 1:length(t_axis)
        t = t_axis(ii);
        g_b(ii) = sinc(t/T)*cos(pi*beta*t/T)/(1-4*beta^2*t^2/T^2); 
    end
    
    figure;
    plot(t_axis,g_b);
    title("Raised Cosine")
     hold on
    scatter(-0.5:T:0.5, zeros(size(-0.5:T:0.5)),"x")
    xlabel("Time(s)")
    ylabel("Magnitude")
    
    X_b = circshift(T_os*fft(circshift(g_b, oversampling_factor/2)), oversampling_factor/2);
    figure;
    plot(-oversampling_factor/2:oversampling_factor/2-1,abs(X_b));
    title("Fourier Transform  - Raised Cosine")
    xlabel("Frequency(Hz)")
    ylabel("Magnitude")
   
    X_sum_b = zeros(1,2*W);
    for ii = 1:10
        X_sum_b = X_sum_b + X_b((ii-1)*100 +1 : ii*100);
    end 
    figure;
    plot(1:100,abs(X_sum_b));
    title("Nyquist Condition Sum Check - Raised Cosine")
    xlabel("Frequency(Hz)")
    ylabel("Magnitude")
    
    % f 
    data_seq = randi([0,1] ,1,20)
    symbols = symbol_mapper(data_seq, 4 , 2, "PSK", "Grey");
    y = pulse_shaper(symbols, "raised cosine", 50);
    
    figure;
    plot((1:length(y))*T_os,y);
    title("Real Pulse Shaper Output - Raised Cosine")
    xlabel("Time(s)")
    ylabel("Magnitude")
    figure;
    plot((1:length(y))*T_os,imag(y));
    title("Imaginary Pulse Shaper Output - Raised Cosine")
    xlabel("Time(s)")
    ylabel("Magnitude")
end


if exper ==3
     W = 50;
     T = 1/(2*W);
     oversample_rate = 1000;
     T_os = 1/oversample_rate;
     data_seq = randi([0,1] ,1,20)
     symbols = symbol_mapper(data_seq, 4 , 2, "PSK", "Grey");
     t_axis = 0:T_os:length(symbols)*T-T_os;
     base_signal = pulse_shaper(symbols, "raised cosine", W);
     carrier_freq = 100;
     
     x = sqrt(2)*real(exp(j*2*pi*carrier_freq*t_axis).*base_signal);
     figure;
     plot(t_axis,x);
     title("Converted Signal- Raised Cosine")
     xlabel("Time(s)")
     ylabel("Magnitude")
     
     SNR = 25;
     x_noise = x + 2*10^(-SNR/20)*randn(size(x));
     figure;
     plot(t_axis,x_noise);
     title("Converted Signal with Noise")
     xlabel("Time(s)")
     ylabel("Magnitude")
     
     dem_i_no_noise = cos(2*pi*carrier_freq*t_axis).*x;
     dem_q_no_noise = sin(2*pi*carrier_freq*t_axis).*x;
     
     dem_i = cos(2*pi*carrier_freq*t_axis).*x_noise;
     dem_q = sin(2*pi*carrier_freq*t_axis).*x_noise;
     
     dem_i_lp = conv(dem_i, sinc(t_axis/T));
     dem_i_lp = dem_i_lp(1:length(x));
     
     dem_q_lp = conv(dem_q, sinc(t_axis/T));
     dem_q_lp = dem_q_lp(1:length(x));
     
     figure;
     plot(t_axis,dem_i_lp);
     title("Filtered Output - In Phase ")
     xlabel("Time(s)")
     ylabel("Magnitude")
     figure;
     plot(t_axis,dem_q_lp);
     title("Filtered Output - Quadrature Phase ")
     xlabel("Time(s)")
     ylabel("Magnitude")
     
end

if exper == 4 
    num_bits = 2000;
    data_bits = randi([0,1],1,num_bits);
    symbols = symbol_mapper(data_bits, 4, 2, "PSK","Binary");
    for SNR = 0:10:20
        symbols_noise = symbols + (randn(size(symbols)) + j*randn(size(symbols)))*10^(-SNR/20);
        figure;
        scatter(real(symbols_noise), imag(symbols_noise));
        title("Scatter of Noisy Symbols SNR=" + SNR)
        xlabel("In Phase")
        ylabel("Quadrature Phase")
    end
    
    demod_bits = symbol_demapper(symbols_noise, 4, 2, "PSK", "Binary", "MD");
    error = sum(abs(demod_bits-data_bits))
    
    c1 = "PAM";
    c2 = "PSK";
    M1 = 16;
    M2 = 16;
    m1 = "Grey";
    m2 = "Grey";
    
    num_bits = 96000;
    data_bits = randi([0,1],1,num_bits);
    symbols1 = symbol_mapper(data_bits, M1, 2, c1, m1);
    symbols2 = symbol_mapper(data_bits,  M2, 2, c2, m2 );
    error_s1 = zeros(1,26);
    error_s2 = zeros(1,26);
    error1 = zeros(1,26);
    error2 = zeros(1,26);
    bps1 = log2(M1);
    bps2 = log2(M2);
    power1 = sqrt(sum(abs(symbols1).^2)/(num_bits/bps1));
    power2 = sqrt(sum(abs(symbols2).^2)/(num_bits/bps2));
    
    figure;
    scatter(real(symbols1), imag(symbols1));
    figure;
    scatter(real(symbols2), imag(symbols2));
    
    n1 = (randn(size(symbols1)) + j*randn(size(symbols1)))*sqrt(0.5);
    n2 = (randn(size(symbols2)) + j*randn(size(symbols2)))*sqrt(0.5);
    
    for SNR = 0:25
        symbols_noise1 =  symbols1 + n1 * power1*10^(-SNR/20);
        symbols_noise2 =  symbols2 + n2 * power2*10^(-SNR/20);
        demod_bits1 = symbol_demapper(symbols_noise1,  M1, 2, c1, m1, "MD");
        demod_bits2 = symbol_demapper(symbols_noise2,  M2, 2, c2, m2, "MD");
        error1(SNR+1) = sum(abs(demod_bits1 - data_bits))/num_bits;
        error2(SNR+1) = sum(abs(demod_bits2 - data_bits))/num_bits;
        
        
        for ii = 1:num_bits/bps1
            if demod_bits1((ii-1)*bps1+1:ii*bps1) == data_bits((ii-1)*bps1+1:ii*bps1)
                
            else
                error_s1(SNR+1)  = error_s1(SNR+1) + 1/(num_bits/bps1);
            end
        end
        for ii = 1:num_bits/bps2
            if demod_bits2((ii-1)*bps2+1:ii*bps2) == data_bits((ii-1)*bps2+1:ii*bps2)
            else
                error_s2(SNR+1)  = error_s2(SNR+1) + 1/(num_bits/bps2);
            end
        end
    end
    
    M3 = 16;
    c3 = "QAM";
    m3 = "Grey";
    
    symbols3 = symbol_mapper(data_bits, M3, 2, c3, m3);
    error_s3 = zeros(1,26);
    error3= zeros(1,26);
    bps3 = log2(M3);
    power3 =  sqrt(sum(abs(symbols3).^2)/(num_bits/bps3));
    figure;
    scatter(real(symbols3), imag(symbols3));
    n3 = (randn(size(symbols3)) + j*randn(size(symbols3)))*sqrt(0.5);
    
    for SNR = 0:25
        symbols_noise3 =  symbols3 + n3 * power3*10^(-SNR/20);
        demod_bits3 = symbol_demapper(symbols_noise3,  M3, 2, c3, m3, "MD");
        error3(SNR+1) = sum(abs(demod_bits3 - data_bits))/num_bits;
           
        for ii = 1:num_bits/bps3
            if demod_bits3((ii-1)*bps3+1:ii*bps3) == data_bits((ii-1)*bps3+1:ii*bps3)
                
            else
                error_s3(SNR+1)  = error_s3(SNR+1) + 1/(num_bits/bps3);
            end
        end
    end
    

    
    figure;
    semilogy(0:25, error1);
    hold on;
    semilogy(0:25, error2);
    semilogy(0:25, error3);
    title("Comparison: "+ M1 + "-" + c1 + " + " + m1 + " Code v.s. " + M2 + "-" + c2 + " + " + m2 + " Code")
    title("Comparison Case 5")
    xlabel("SNR(dB)")
    ylabel("Bit Error Rate (log scale)")
    legend(M1 + "-" + c1 + " + " + m1 ,  M2 + "-" + c2 + " + " + m2,  M3 + "-" + c3 + " + " + m3)
    
    figure;
    semilogy(0:25, error_s1);
    hold on;
    semilogy(0:25, error_s2);
    semilogy(0:25, error_s3);
    title("Comparison: "+ M1 + "-" + c1 + " + " + m1 + " Code v.s. " + M2 + "-" + c2 + " + " + m2 + " Code")
    title("Comparison Case 5")
    xlabel("SNR(dB)")
    ylabel("Symbol Error Rate (log scale)")
    legend(M1 + "-" + c1 + " + " + m1 ,  M2 + "-" + c2 + " + " + m2,  M3 + "-" + c3 + " + " + m3)
end

if exper ==5
    [x,fs] = audioread("handel.ogg");
    x = x';
    x_quan = quantizer_L_level(x,max(x), 16);
    x_symbols = unique(x_quan);
    symbol_count = zeros(1,16);
    for ii = 1:length(x_quan)
        j = find(x_symbols==x_quan(ii));
        symbol_count(j) = symbol_count(j) + 1;
    end
    p = symbol_count/length(x_quan);
    
    dict = huffman_dict(x_symbols,p);
    binary_data = pcm_enc(x_quan, 4);
    
    symbols = symbol_mapper(binary_data, 4, sqrt(2) ,"PSK", "Grey");
    
    SNR = 0;
    symbols_noise = symbols + randn(size(symbols))*10^(-SNR/20);
    binary_demod = symbol_demapper(symbols_noise, 4 , sqrt(2), "PSK", "Grey", "MD");
    
    x_dec = pcm_dec(binary_demod, x_symbols  , 4);
    
    sum(abs(binary_demod- binary_data));
    
%     figure;
%     plot(1/fs*(1:length(x)) , x_dec);
%     title("Recovered Signal with SNR = " + SNR)
%     xlabel("Time(s)")
%     ylabel("Magnitude")
    
end
