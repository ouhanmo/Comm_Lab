function dict = huffman_dict(symbols,p)
dict = cell(length(symbols),2);
merging = zeros(1,length(symbols));
merge_seq = zeros(1,length(symbols));

for iteration = 1:length(symbols)-1
    [merged_prob,merged_index] = min(p);
    p(merged_index) = 2;
    [~,merge_index] = min(p);
    p(merge_index) =  p(merge_index) + merged_prob;
    merging(merged_index) = merge_index;
    merge_seq(iteration) = merged_index;
end

[~,merge_seq(end)] = min(p);

for ii = length(symbols)-1:-1:1
    right = merge_seq(ii);
    left  = merging(right);
    dict{right, 2} = [dict{left,2} 1];
    dict{left, 2} = [dict{left,2} 0];
end 

for ii = 1:length(symbols)
    dict(ii,1) = {symbols(ii)};
end
