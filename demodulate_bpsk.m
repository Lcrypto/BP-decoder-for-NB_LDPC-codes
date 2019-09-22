function [data] = demodulate_bpsk(Q, cdata, sigma)
   data = [];
   
   hm = modem.pskmod('M', 2);
   data_to_compare = real(modulate(hm, de2bi([0:1:Q-1], log2(Q))));
   
   for i = 1:length(cdata)/log2(Q)
      input = cdata(log2(Q)*(i-1)+1:log2(Q)*i);
      temp = data_to_compare - repmat(input, Q, 1);
      temp = temp.^2;
      square_distances = sum(temp, 2);
      llr_one_symbol = -square_distances/(2*sigma^2);
      llr_one_symbol = llr_one_symbol - llr_one_symbol(1);
      data = [data llr_one_symbol];
   end
   
end