close all;
clear;
clc;

exper = 4;

if (exper == 1)
    % 1
    fs = 400;
    Ts = 1/fs;
    T  = Ts*[1:4/Ts];
    frequencies = [10,25,50,100];
    x = [];
    for ii = 0:3
        x = [x cos(T(1/Ts*ii+1:1/Ts*ii+1/Ts)*2*pi*frequencies(ii+1))];

    end

    ret_window1 = [ones(1,400),zeros(1,1200)];
    x1 = x.*ret_window1;

    x1_w = fft(x1);
    x1_w = circshift(x1_w,801);
    freq = (-799:800)*1/2*fs/800;
    figure;
    plot(freq,abs(x1_w));
    xlabel("frequency(Hz)")
    ylabel("CTFT")

    ret_window2 = [zeros(1,1000),ones(1,400),zeros(1,200)];
    x2 = x.*ret_window2;

    x2_w = fft(x2);
    x2_w = circshift(x2_w,801);
    figure;
    plot(freq,abs(x2_w));
    xlabel("frequency(Hz)")
    ylabel("CTFT")
    
    window_size = 50;
    [stft, f, t ] = STFT(x,ones(1,window_size),4,window_size,fs);
    figure;
    surf(t,f,abs(stft));
    xlabel("time(s)")
    ylabel("frequency(Hz)")
    zlabel("STFT")

    [x,fs] = audioread("handel.ogg");
    x = x';
    
    x = quantizer_L_level(x, max(x),4);
    [stft, f, t ] = STFT(x,ones(1,5000),2000,5000,fs);
    figure;
    surf(t,f,abs(stft));
    xlabel("time(s)")
    ylabel("frequency(Hz)")
    zlabel("STFT")
end 

if exper==2
    symbol = ['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h'];
    p = [0.3 0.2 0.2 0.1 0.05 0.05 0.05 0.05];
    %dict = huffmandict(symbol,p);
    dict_my = huffman_dict(symbol, p);
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
    
    figure;
    plot(symbols,p);
    xlabel("symbol");
    ylabel("probabilities");
    
    dict = huffman_dict(symbols,p);
    y = huffman_enc(x_quan, dict);
    x_dec = huffman_dec(y,dict);
    
    y_pcm = pcm_enc(x_quan,4);
    x_dec_pcm = pcm_dec(y_pcm, symbols, 4);
end

if exper == 4
    p = [0.3 0.2 0.2 0.1 0.05 0.05 0.05 0.05];
    entropy = 0;
    for ii = 1:length(p)
        entropy = entropy - p(ii)*log2(p(ii));
    end
   
    
    dict = huffman_dict(0:7, p);
    
    avg_length = 0;
    for codeword = 1:length(p)
        avg_length = avg_length + length(dict{codeword,2})*p(codeword);
    end
   
    
    R = 1000;
    L_sum = 0;
    for sample_num = 1:R
        seq = randsrc(1,100,[0:7 ; p]);
        seq_huf = huffman_enc(seq, dict);
        L_sum = L_sum +  length(seq_huf);
    end
    avg = L_sum/R
    
    R = 10000;
    L_sum = 0;
    for sample_num = 1:R
        seq = randsrc(1,1000,[0:7 ; p]);
        seq_huf = huffman_enc(seq, dict);
        L_sum = L_sum +  length(seq_huf);
    end
    avg = L_sum/R
    
end

