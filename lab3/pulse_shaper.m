function y = pulse_shaper(x, pulse_shape, W)
T = 1/(2*W);
W_os = 1000;
T_os = 1/W_os;

time = -T/2:T_os:T/2-T_os;

switch pulse_shape
    case "sinc"
        g = sinc(time/T);
                      
    case "raised cosine"
        beta = 0.5;
        g = zeros(1,length(time));
        for ii = 1:length(time)
        t = time(ii);
        g(ii) = sinc(t/T)*cos(pi*beta*t/T)/(1-4*beta^2*t^2/T^2); 
        end
    case "rect"
        g = ones(1,length(time));
end

y = zeros(1, length(time)*length(x));

G = T_os/T*fft(g);
p = fft(sqrt(abs(G)));


for ii = 0 : length(x)-1
    y(ii*length(time)+1:(ii+1)*length(time)) = p*x(ii+1);
end