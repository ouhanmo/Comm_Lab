function  x_dec = huffman_dec(y,dict)
x_dec = [];

start = 1;
finish = 1;

while finish <=length(y)
    code_cand = y(start:finish);
    confirmed_symbol = 0;
    for symbol_id = 1 : length(dict)
        if isequal(dict{symbol_id, 2},code_cand)
            confirmed_symbol = symbol_id;
            break
        end
    end
    if confirmed_symbol > 0
        x_dec = [x_dec dict{confirmed_symbol, 1}];
        start = finish +1;
        finish = start;
    else
        finish = finish + 1;
    end
    
end
    