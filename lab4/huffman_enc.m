function y = huffman_enc(x,dict)
y = [];

set = [];
for ii =1:length(dict)
   set(ii) = dict{ii,1};
end

for ii = 1 : length(x)
    symbol_ind = find(set==x(ii));
    y = [y dict{symbol_ind,2}];
end 







