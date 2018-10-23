function [Gp] = gperiodic_coeff_nop_linJVIE(cube)
%%
% The following is currently being updated due to the signs of 
% new basis functions
% L,M,N
if     strcmp(cube, 'L')
    
    Gp   = [ +1.0 -1.0 +1.0 +1.0 +1.0 +1.0 +1.0];
         
elseif strcmp(cube, 'M')
     
    Gp   = [ +1.0 +1.0 -1.0 +1.0 +1.0 +1.0 +1.0];
           
elseif strcmp(cube, 'N')
    
    Gp   = [ +1.0 +1.0 +1.0 +1.0 -1.0 +1.0 +1.0];
           
elseif strcmp(cube, 'LM')
    
    Gp   = [ +1.0 -1.0 -1.0 +1.0 +1.0 +1.0 +1.0];
        
elseif strcmp(cube, 'LN')
    
    Gp   = [ +1.0 -1.0 +1.0 +1.0 -1.0 +1.0 +1.0];
         
elseif strcmp(cube, 'MN')    
    
    Gp   = [ +1.0 +1.0 -1.0 +1.0 -1.0 +1.0 +1.0];

elseif strcmp(cube, 'LMN')
    
    Gp   = [ +1.0 -1.0 -1.0 +1.0 -1.0 +1.0 +1.0];
      
end
Gp = Gp';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%