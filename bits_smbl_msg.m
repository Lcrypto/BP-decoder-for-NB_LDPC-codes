function [pb0, pb1, py] = bits_smbl_msg(b0,b1,y)
% Message passing bitween symbol and its bits
%
%   Use:     [pb0, pb1, py] = bits_smbl_msg(b0,b1,y)
%
%   b0,b1 (i,t) - proportional to probability that i-th bit is 0 and 1
%                 (may be not normalized) at time t (source time)
%   y(j,t)      - proportional to prob that y=j at time  t
%
%   Example let Y be a symbol with the value from 0 to 7.
%   let b1 b2 b3 be its binary representation
%   i.e.   b1*4 + b2*2 + b3 = Y
%   then given incoming messages P* ("apriori distributions" b0, b1 and y)
%   we compute output messages (pb0, pb1 and py) for each branch by multipling
%   the incoming distribution with the local function and
%   summarizing it for that branch

% Igor Kozintsev, igor.v.kozintsev@intel.com


    pb0 = zeros( size( b0 ) );
    pb1 = zeros( size( b1 ) );
    py =  zeros( size( y  ) );

    nbits = size(b0,1);
    T =     size(b0,2);

    for i = 1 : size(y,1)
        %convert to binary format
        x=[]; for j=1:nbits, x(j,1) = bitget(i-1,nbits+1-j); end;        
        x = x(:,ones(1,T));

        %message to y variable
        py(i,:) = prod(x.*b1 + (1-x).*b0);

        %accumulate messages for bits
        h = py(i,:) .* y(i,:);
        h = h(ones(1,nbits),:); 
        denom = (x.*b1 + (1-x).*b0); 
        denom(find(denom == 0))=realmin; 
        
        pb0 = pb0 + h.*(1-x)./denom;
        pb1 = pb1 + h.*x./denom;
      
    end

    %normalize
    ppy = sum(py); py = py ./(ppy(ones(1,size(y,1)),:));    
    ppb = pb0 +pb1;
    pb0 = pb0./ppb; pb1 = pb1./ppb;

