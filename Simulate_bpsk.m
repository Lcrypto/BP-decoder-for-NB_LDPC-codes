

%Copyright(c) 2013, Vasiliy Usatyuk
%All rights reserved.

%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met :
%*Redistributions of source code must retain the above copyright
%notice, this list of conditions and the following disclaimer.
%* Redistributions in binary form must reproduce the above copyright
%notice, this list of conditions and the following disclaimer in the
%documentation and / or other materials provided with the distribution.
%* Neither the name of the <organization> nor the
%names of its contributors may be used to endorse or promote products
%derived from this software without specific prior written permission.

%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED.IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
%DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
      
    %parpool(8);
    clear all;
    s = RandStream('mt19937ar','Seed','shuffle'); %use Mersenne twister 
    RandStream.setGlobalStream(s);
    ebn0=1;
    BER=0;  
    SER=0;
    FER=0;
    
    Codewords=0;
    FER_error=100; %number of error after which we stop simulation
    Rate=0.5; %code rate
    q = 64;                              % Field parameter 
    nbits = log2(q);                    % bits per symbol 
    
    %h = qc2sparse('6x12x198_2.txt');
    %[row,col,v] = find(h) ;
	%Number=size(v); 
    %for i=1:Number
   	%	h(row(i),col(i))=randi([1,q-1],1);
    %end
    
    %from alist format (for example generated by PEG)
      
    %h = qc2sparse('3x6x460_GF16.txt');
    %h = qc2sparse('4x8x50.txt');
    %h = qc2sparse('6x4x360ir_gF16.txt');
   %colmnelement=1;
   % [row, col]=size(h);
   %  for i=1:row
   %     row_in_optimal_degree=randi([1,4],1); % ramdom choice of optimal labelin for layer from table of optimal coefficient
   %     for j=1:col
   %         if h(i,j)~=0
   %             h(i,j)=optimal_row_labeling(row_in_optimal_degree,colmnelement);
   %             colmnelement=colmnelement+1;
   %         end      
   %      end
   %     colmnelement=1;
   % end
    
%     [row,col,v] = find(h) ;
% 	Number=size(v);
%     for i=1:Number
%    		h(row(i),col(i))=randi([1,q-1],1);
%     end
    
    h = alist2sparse('Dv2Dc4_G16_N160');
    [row,col,v] = find(h) ;
	Number=size(v);
    
	for i=1:Number
   		h(row(i),col(i))=randi([1,q-1],1);
    end
    %optimal_row_labeling array must have optimal edge labeling for CNs
    %h = qc2sparse('2x4x40_QC_proto_for_NB.txt');
   
    %{
    h = alist2sparse('Dv2Dc4_G16_N160.txt');
    colmnelement=1;
    [row, col]=size(h);
     for i=1:row
        row_in_optimal_degree=randi([1,4],1); % ramdom choice of optimal labelin for layer from table of optimal coefficient
        for j=1:col
            if h(i,j)~=0
                h(i,j)=optimal_row_labeling(row_in_optimal_degree,colmnelement);
                colmnelement=colmnelement+1;
            end      
         end
        colmnelement=1;
    end
    
    %}
    
    [H,G] = h2g(h,q);              % find systematic G and modify H (GE)
    
%     check h2g
%     R=gf(full(G),q)*gf(full(H'),q);
%     count=0;
%     for i=1:size(G,1)
%     for j=1:size(H,2)
%     if R(i,j)~=0
%     count=count+1
%     sprintf('Trouble with H2G');
%     end
%     end
%     end
    
    sparse2alist(q,H,'simulate.alist');
    init_ldpc = @(k) decode_soft(0, k);
    [ldpc, ~, ~] = init_ldpc('simulate.alist');

    [K, N] = size(G);
    Rate = K/N;% rate after linear dependency delete 
    %Rate = 0.57;
    sigma = sqrt( 0.5 / Rate/ power( 10, ebn0 / 10 ) );    % AWGN noise deviation 
    snr = -10*log(2*sigma*sigma)/log((10));
  
    
    while(FER<FER_error)
    Codewords=Codewords+1;              % decoded codewords
    x = randi(q, 1, K) - 1;    % random symbols 
    %x=zeros(1,size(G,1));             %zero codeword
    
    
    y = ldpc_encode(x,G,q);                 % encoding  
    yb=((de2bi(y,nbits)))';
    yb=yb(:);
    nm=modem.pskmod('M',2);
    zb=real(modulate(nm,yb));
    zb=zb+sigma*randn(size(zb));
    
    
    
    pp = demodulate_bpsk(q, zb', sigma);
    
    decode_ldpc = @(decoder, ldpc, in_llrs, iter)decode_soft(2, ldpc, pp, 30);
    [result, number_of_iter, est_cwd, out_llrs] = decode_ldpc(2, ldpc,pp , 30);
    z_hat=est_cwd';
     %if result ~=1
     %      sprintf('Syndrom is ok');
     %send
    x_hat = z_hat(size(G,2)+1-size(G,1):size(G,2)); 
    x_hat = x_hat';    
        if size(find(x~= x_hat),2) ~=0     
            FER=FER+1;
            SER=SER+size(find(x~= x_hat),2);% number of errors in symbol after decoding
        end
        
    b=(fliplr(de2bi(x,nbits)))';     % binary x
	b=b(:);    
	a=(fliplr(de2bi(x_hat,nbits)))'; %  binary result of encoding
	a=a(:);
	BER=BER+size(find(a~=b),1);     %compare bit representation except parity-check bits
	BERstat=BER/(size(H,2)*nbits*Codewords)
    FERstat=FER/Codewords
    Codewords
    FER
    ebn0=snr-10*log10(Rate)
    end
    