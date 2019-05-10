function x_dec = pcm_dec(y, symbols, numBits)
x_dec = [];

for ii = 0:length(y)/numBits-1
    symbol_num = bi2de(y(ii*numBits+1:ii*numBits+numBits))+1;
    x_dec = [x_dec  symbols(symbol_num)];
    
end
