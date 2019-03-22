close all;
clc;
clear;

exper = 1;

% 1
[x,fs] = audioread("handel.ogg");
x = x';
time = 1/fs*(1:length(x));
if (exper == 1)
    sound(x,fs);
end
audiowrite('output_file.flac', x, fs);
%sound(x,fs/2);
%sound(x,fs*2);
figure;
plot(time,x);
xlabel("Time (s)");
ylabel("Signal");

close all;

% 2 
if ( exper ==  2 )   

    x_circ = circshift(x, 100000);
    %sound(x_circ,fs);
    figure;
    plot(time,x_circ);
    xlabel("Time (s)");
    ylabel("Signal");

    x_rev = fliplr(x);
    %sound(x_rev,fs);
    figure;
    plot(time,x_rev);
    xlabel("Time (s)")
    ylabel("Signal");

    x_forNback = [x  x_rev];
    time_2 = (1:length(x_forNback))/fs;
    %sound(x_forNback,fs);
    figure;
    plot(time_2,x_forNback);
    xlabel("Time (s)")
    ylabel("Signal");

    upsample_rate = 2;
    x_up = zeros(1,length(x)*upsample_rate);
    for ind = 1:length(x)
        x_up(ind*upsample_rate) = x(ind);
    end
    time_up = 1/upsample_rate/fs*(1:length(x)*upsample_rate);
    %sound(x_up,fs*upsample_rate);
    figure;
    plot(time_up, x_up);
    xlabel("Time (s)")
    ylabel("Signal")

    downsample_rate = 2;
    x_down = x(downsample_rate*(1:length(x)/downsample_rate));
    %sound(x_down,fs/downsample_rate);
    time_down = downsample_rate/fs*(1:length(x)/downsample_rate);
    figure;
    plot(time_down, x_down);
    xlabel("Time (s)")
    ylabel("Signal")

end


%3
if (exper == 3)


    for T = [0.1, 0.5, 1, 0.01, 0.05 ]
        x_clip = max(min(T,x),-T);
        figure;
        if (T == 0.5)
            sound(x_clip,fs);
        end 
        plot(time,x_clip);
        xlabel("Time (s)")
        ylabel("Signal")
    end

    x_sq = x.*x;
    figure;
    %sound(x_sq,fs);
    plot(time,x_sq);
    xlabel("Time (s)")
    ylabel("Signal")

    x_neg = -x;
    figure;
    %sound(x_neg,fs);
    plot(time,x_neg);
    xlabel("Time (s)")
    ylabel("Signal")

end 

%4 

if (exper == 4)
    x_quan4 = quantizer_L_level(x,max(x),4);
    figure;
    %sound(x_quan4,fs);
    plot(time,x_quan4);
    xlabel("Time (s)")
    ylabel("Signal")
    
    x_quan2 = quantizer_L_level(x,max(x),2);
    figure;
    %sound(x_quan2,fs);
    plot(time,x_quan2);
    xlabel("Time (s)")
    ylabel("Signal")
     
    x_quan8 = quantizer_L_level(x,max(x),8);
    figure;
    %sound(x_quan8,fs);
    plot(time,x_quan8);
    xlabel("Time (s)")
    ylabel("Signal")
    
    x_quan15 = quantizer_L_level(x,max(x),15);
    figure;
    %sound(x_quan15,fs);
    plot(time,x_quan15);
    xlabel("Time (s)")
    ylabel("Signal")

end

%5
if (exper == 5 )
    fc = 100;
    mod_signal = exp(time*2*pi*fc*i);
    x_mod = x.*mod_signal;
    figure;
    sound(real(x_mod),fs);
    plot(time(1:2000),abs(10*x_mod(1:2000)));
    hold on
    plot(time(1:2000),angle(x_mod(1:2000)));
    xlabel("Time (s)")
    ylabel("Signal")
end

if (exper == 6)
    var = 0.01; 
    x_noise = x + sqrt(var)*randn(1,length(x));
    sound(x_noise,fs);
    figure;
    plot(time,x_noise);
    xlabel("Time (s)")
    ylabel("Signal")
    
end

if (exper == 7)
    f_low  = 94;
    f_high = 142;
    W = 50;
    filter_design = [f_low/(0.5*fs),f_high/(0.5*fs)];
    filter_male = fir1(W, filter_design);
    
    figure;
    stem(0:W,filter_male);
    x_male = filter (filter_male,1,x);
    xlabel("n")
    ylabel("Impulse Response")
    
    sound(x_male,fs);
    figure;
    mag_resp = abs(freqz(filter_male,1,2000));
    plot(1/2000*fs/2*(1:200),mag_resp(1:200));
    xlabel("Frequency(Hz)")
    ylabel("Magnitude Response")
    
    figure;
    plot(time,x_male);
    xlabel("Time (s)")
    ylabel("Signal")
end 

if (exper == 8)
        clear;
        fs = 8192;
        fc = 1200;
        time = 1/fs*(1:200);
        x = cos(2*pi*fc*time);
        %sound(x,fs);
        figure;
        plot(time,x);
        xlabel("Time (s)")
        ylabel("Signal")

        % DFT
        X = zeros(1,length(x));
        for i1 = 1:length(x)
            for i2 = 1:length(x)
               X(i1) = X(i1) + x(i2)*exp(-j*2*pi*i1*i2/length(x));
            end
        end
        figure;
        stem(abs(X));
        xlabel("k")
        ylabel("X[k]")
    
        % DTFT
        z = 2*pi*(1:length(x))/length(x);
        z = circshift(z,length(x)/2);
        z(1:length(x)/2) = z(1:length(x)/2) -2*pi;
        X_dtft = circshift(X,length(x)/2-1);
        figure;
        plot(z,abs(X_dtft));
        xlabel("w")
        ylabel("X(exp(jw))")

        % CTFT
        Omega = fs*(z);
        X_ctft = X_dtft/fs;
        figure;
        plot(Omega/2/pi,abs(X_ctft));
        xlabel("frequency(Hz)")
        ylabel("X(jw)")
end 




