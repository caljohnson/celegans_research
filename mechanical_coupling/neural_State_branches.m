function [ y ] = neural_State_branches( x, thresh1, thresh2 )
%NEURAL_STATE_BRANCHES 
  persistent b;
  
  if( isempty(b) )
    b = 0;
  end
  
  if( x > thresh2 )
    y = 1;
    b = 1;
  elseif( x < thresh1 )
    y = 0;
    b = 0;
  else
    y=b;
  end
end

