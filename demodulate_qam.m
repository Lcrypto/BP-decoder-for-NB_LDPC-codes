function [result] = demodulate_qam(h, cdata, sigma)
    result = [];
    M = h.M;
    data_to_compare = modulate(h, [0:1:M-1]);
    for i = 1:length(cdata)
        square_distances = (cdata(i) - data_to_compare).*conj(cdata(i) - data_to_compare); 
        llr_one_symbol = -square_distances/(2*sigma^2);
        llr_one_symbol = llr_one_symbol - llr_one_symbol(1);
        result = [result llr_one_symbol'];
    end
end