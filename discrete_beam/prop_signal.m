function [ P ] = prop_signal( n, m, G, K, sig )
%PROP_SIGNAL creates the proprioceptive signals
%  
    P = zeros(n-2,1);
    %for neurons at the HEAD, less proprioceptive coupling since fewer
    %anterior neurons
    for ii= 1:m
       %guassian kernel for proprioception - need new one for each size
       grid2 = -ii:ii;
       G2 = exp(-grid2.^2./(2*sig^2)); %Gaussian kernel of width 2m
       G2 = 2*G2/sum(G2); %normalize so half of G has area 1   
       for jj=0:ii-1
            P(ii) = P(ii) + K(ii-jj)*G2(ii-jj);
            %signal filtered with Gaussian kernel
       end
    end
    
    %for neurons m posterior from the head, m anterior neurons yield max
    %proprioceptive coupling
    for ii = 1+m:n-2
        for jj = 0:m-1
            P(ii) = P(ii) + K(ii-jj)*G(m-jj);
            %signal filtered with Gaussian kernel
        end
    end

end

