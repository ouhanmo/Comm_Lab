function binary_sequence = symbol_demapper(symbol_sequence, M, d, constellation, mapping, decision_rule)
bits_per_symbol = floor(log2(M));
binary_sequence = zeros(1,length(symbol_sequence)*bits_per_symbol);

possible_bits = de2bi(0:M-1);
possible_symbols = symbol_mapper(reshape(de2bi(0:M-1)',1,[]), M , d , constellation, mapping);

if decision_rule == "MD"
    for ii = 1:length(symbol_sequence)
        symbol = symbol_sequence(ii);
        for candidate = 1:M
            distance(candidate) = abs(possible_symbols(candidate) - symbol);
        end
        [~, md ] = min(distance);
        binary_sequence((ii-1)*bits_per_symbol +1 : ii*bits_per_symbol) = possible_bits(md,1:end);
    end
end