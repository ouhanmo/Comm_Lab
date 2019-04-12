function y =pcm_enc(x,numBits)
y = [];
symbols = unique(x);
dict = de2bi(0:2^numBits-1);

for ii = 1:length(x)
    y = [y dict(find(symbols==x(ii)),1:numBits) ];
end

