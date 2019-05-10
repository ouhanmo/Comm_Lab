function symbol_sequence = symbol_mapper(binary_sequence, M, d, constellation, mapping)
%{
? Inputs:
¡V M: The number of points in the signal constellation and M = 2,4,8,16 for 'PAM'
and 'PSK', M = 4,16 for 'QAM'.
¡V d: The minimum distance among the constellation.
¡V constellation: 'PAM', 'PSK' or 'QAM'.
¡V mapping: 'Binary' or 'Gray'.
For Gray code mapping of PAM, please assign 00 ¡P ¡P ¡P 0 to the leftest constellation point.
Gray code mapping of PSK and QAM is shown in Figure 1 and Figure 2. If input M is
not the power of 2 or if input constellation is not defined, your function should be able
to throw error and display message.
%}

symbol_sequence = zeros(1,length(binary_sequence)/log2(M));

switch constellation
    case "PAM"
        if M == 2 || M == 4 || M == 8 ||  M == 16
            dict = d*(-M+1 : 2: M-1)
        else 
            except = MException("SymbolMapper:ConstellationError", "Number of points not valid");
            throw(except);
        end
    case "PSK"
         radius = d/2/sin(pi/M);
         switch  M
             case 2
                 dict = radius*([-1,1])
                 
             case {4,8} 
                 dict = radius*(exp(-j*2*pi/M*(0:M-1)-j*3*pi/4 ))
                 
             case 16  
                 dict = radius*(exp(j*2*pi/M*(0:M-1)))
             
             otherwise
                 except = MException("SymbolMapper:ConstellationError", "Number of points not valid");
                 throw(except);
        end
    case "QAM"
         switch M
             case {4,16}
                 re  =  d/2*(-sqrt(M)+1 : 2: sqrt(M)-1);
                 im  =  d/2*(-sqrt(M)+1 : 2: sqrt(M)-1);
                 for ii = 1:M
                     dict(ii) = re(floor((ii-1)/sqrt(M))+1) + j*im( mod((ii-1),sqrt(M))+1);
                 end
                dict
             otherwise
                 except = MException("SymbolMapper:ConstellationError", "Number of points not valid");
                 throw(except);
         end
    otherwise
        except = MException("SymbolMapper:ConstellationError", "No such constellation type");
        throw(except);

end



switch M
    case 2
        grey_seq = [1 2];
        
    case 4 
        grey_seq = [1 2 4 3];
        if constellation == "QAM"
            grey_seq = [1 2 3 4];
        end
        
    case 8
        grey_seq = [1 2 4 3 8 7 5 6];
    
    case 16
        switch constellation
            case "PAM"
                grey_seq = [1 2 4 3 8 7 5 6 16 15 13 14 9 10 12 11];
            case "PSK"
                grey_seq = [1 2 4 3 8 7 5 6 16 15 13 14 9 10 12 11];
            case "QAM"
                grey_seq = [1 2 4 3 5 6 8 7 13 14 16 15 9 10 12 11];

        end
end




switch mapping
    case "Binary"
        for ii = 1:length(binary_sequence)/log2(M)
            bits = binary_sequence(ii*log2(M)-log2(M)+1:ii*log2(M));
            symbol_sequence(ii) = dict(bi2de(bits)+1);
        end
    case "Grey"
        for ii = 1:length(binary_sequence)/log2(M)
            bits = binary_sequence(ii*log2(M)-log2(M)+1:ii*log2(M));
            symbol_sequence(ii) = dict(grey_seq(bi2de(bits)+1));
        end
        dict = dict(grey_seq)
    otherwise
        except = MException("SymbolMapper:MappingError", "No such mapping style");
        throw(except);    
end

return
%plot


figure;
scatter(real(dict), imag(dict));
axis([-2 2 -2 2])
switch M
    case 2
        labels = ["0" "1"]
    case 4
        labels = ["00" "01" "10" "11" ];
    case 8
        labels = ["000" "001" "010" "011" "100" "101" "110" "111" ];       
    case 16
        labels = ["0000" "0001" "0010" "0011" "0100" "0101" "0110" "0111" "1000" "1001" "1010" "1011" "1100" "1101" "1110" "1111" ];
end
text(real(dict), imag(dict),labels,'VerticalAlignment','bottom','HorizontalAlignment','right');
title(  M + "-" + constellation + " Constellation Map - " + mapping)
xlabel("In Phase")
ylabel("Quadrature Phase")
