
function [JOut_full_out]=lse_sparse_precon_multiply(JOut_full_in,Ae,nodeid_4_grnd,nodeid_4_injectcurr)

global A_inv LL UU PP QQ RR Sch_sparse slct_decomp_sch fl_cholmod
global fl_sparseprecon

if (fl_sparseprecon == 0)
    JOut_full_out = JOut_full_in;
    return
end

% ---------------------------------------------------------------------
% Sparse preconditioner [E F; G H]
% ---------------------------------------------------------------------

fl_fast_multiply = 1; % for single inversion of Schur complement
fl_volt_source = 2; % symmetric voltage source implementation

num_node=size(Ae,1);
num_curr=size(Ae,2);
JOut_full_out=zeros(num_node+num_curr,1);

switch slct_decomp_sch
    
    case 'lu_decomp'
        
        if (fl_volt_source == 1) % voltage source
            Ae_tmp=Ae; % Attention::: temporarily doubling the memory for Ae!!!
            Ae_tmp(nodeid_4_grnd,:)=0;
            Ae_tmp(nodeid_4_injectcurr,:)=0;
            
            if (fl_fast_multiply == 1)
                % for inverting system [A_inv B; C D][x;y]=[a;b]
                % Solve for y of (D - C A_inv B)y = b-C(A_inv)a;
                % Get potentials first. Then obtain currents via
                % (A_inv) x = a-By
                
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(num_curr+1:num_curr+num_node) - ( Ae_tmp * A_inv * JOut_full_in(1:num_curr) ) )))));
                
                JOut_full_out(1:num_curr) = (JOut_full_in(1:num_curr) - ((-Ae')*JOut_full_out(num_curr+1:num_curr+num_node))) .* diag(A_inv);
                
            else
                
                % block E contribution
                JOut_full_out(1:num_curr) = A_inv*JOut_full_in(1:num_curr)+A_inv*(-Ae')*...
                    QQ * (UU \ (LL \ (PP * (RR \ (Ae_tmp*A_inv*JOut_full_in(1:num_curr))))));
                
                % block F contribution
                JOut_full_out(1:num_curr) = JOut_full_out(1:num_curr)...
                    + A_inv * (Ae') * QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(num_curr+1:num_curr+num_node))))));
                
                % block G contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    -QQ * (UU \ (LL \ (PP * (RR \ (Ae_tmp*A_inv*JOut_full_in(1:num_curr))))));
                
                % block H contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = JOut_full_out(num_curr+1:num_curr+num_node)...
                    +QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(num_curr+1:num_curr+num_node))))));
            end
            
        else % current source or symmetic voltage source
            
            %d = S^-1*(b - Ar*Y^-1*a)
            %c = Y^-1*(a + Ar'*d)
            
            if (fl_fast_multiply == 1)
                
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(num_curr+1:num_curr+num_node) - ( Ae * A_inv * JOut_full_in(1:num_curr) ) )))));
                
                JOut_full_out(1:num_curr) = (JOut_full_in(1:num_curr) - ((-Ae')*JOut_full_out(num_curr+1:num_curr+num_node))) .* diag(A_inv);
                
                
            else
                
                % block E contribution
                JOut_full_out(1:num_curr) = A_inv*JOut_full_in(1:num_curr)+A_inv*(-Ae')*...
                    QQ * (UU \ (LL \ (PP * (RR \ (Ae*A_inv*JOut_full_in(1:num_curr))))));
                
                % block F contribution
                JOut_full_out(1:num_curr) = JOut_full_out(1:num_curr)...
                    + A_inv * (Ae') * QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(num_curr+1:num_curr+num_node))))));
                
                % block G contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    -QQ * (UU \ (LL \ (PP * (RR \ (Ae*A_inv*JOut_full_in(1:num_curr))))));
                
                % block H contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = JOut_full_out(num_curr+1:num_curr+num_node)...
                    +QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(num_curr+1:num_curr+num_node))))));
                
            end
            
        end
        
        
    case 'ldlt_decomp'
        
        if (fl_volt_source == 1) % voltage source
            
            disp('Change decomposition !!! ')
            disp('LDLT decomposition can not work with voltage source !!!')
            error('as the Schur complement is not symmetric for such case ... ')
            
        else % current source or symmetic voltage source
            
            if (fl_fast_multiply == 1)
                
                if (fl_cholmod == 1)

                    %JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    %    (RR * PP * (LL' \ (QQ \ (LL \ (PP' * RR * (JOut_full_in(num_curr+1:num_curr+num_node) - ( Ae * A_inv * JOut_full_in(1:num_curr) ) ))))));
                    
                    bb_dum = JOut_full_in(num_curr+1:num_curr+num_node) - (Ae * A_inv * JOut_full_in(1:num_curr));
                    
                    JOut_full_out(num_curr+QQ) = ldlsolve (LL,bb_dum(QQ));
                    
                    JOut_full_out(1:num_curr) = (JOut_full_in(1:num_curr) - ((-Ae')*JOut_full_out(num_curr+1:num_curr+num_node))) .* diag(A_inv);
                    
                    
                else
                    
                    JOut_full_out(num_curr+1:num_curr+num_node) = ...
                        (RR * PP * (LL' \ (QQ \ (LL \ (PP' * RR * (JOut_full_in(num_curr+1:num_curr+num_node) - ( Ae * A_inv * JOut_full_in(1:num_curr) ) ))))));
                    
                    JOut_full_out(1:num_curr) = (JOut_full_in(1:num_curr) - ((-Ae')*JOut_full_out(num_curr+1:num_curr+num_node))) .* diag(A_inv);
                end
            else
                
                % block E contribution
                JOut_full_out(1:num_curr) = A_inv*JOut_full_in(1:num_curr)+A_inv*(-Ae')*...
                    (RR * PP * (LL' \ (QQ \ (LL \ (PP' * RR * (Ae*A_inv*JOut_full_in(1:num_curr)))))));
                
                % block F contribution
                JOut_full_out(1:num_curr) = JOut_full_out(1:num_curr)...
                    + A_inv * (Ae') * (RR * PP * (LL' \ (QQ \ (LL \ (PP' * RR * (JOut_full_in(num_curr+1:num_curr+num_node)))))));
                
                % block G contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    - RR * PP * (LL' \ (QQ \ (LL \ (PP' * RR * (Ae*A_inv*JOut_full_in(1:num_curr))))));
                
                % block H contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = JOut_full_out(num_curr+1:num_curr+num_node)...
                    + (RR * PP * (LL' \ (QQ \ (LL \ (PP' * RR * (JOut_full_in(num_curr+1:num_curr+num_node)))))));
                
            end
            
        end
        
    case 'chol_decomp'
        
        if (fl_volt_source == 1)
            
            disp('Change decomposition !!! ')
            disp('Cholesky decomposition can not work with voltage source !!!')
            error('as the Schur complement is not symmetric for such case ... ')
            
        else % current source or symmetic voltage source
            
            if (fl_fast_multiply == 1)
                
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    (PP * (LL' \ (LL \ (PP' * (JOut_full_in(num_curr+1:num_curr+num_node) - ( Ae * A_inv * JOut_full_in(1:num_curr) ) )))));
                
                JOut_full_out(1:num_curr) = (JOut_full_in(1:num_curr) - ((-Ae')*JOut_full_out(num_curr+1:num_curr+num_node))) .* diag(A_inv);
                
            else
                
                % block E contribution
                JOut_full_out(1:num_curr) = A_inv*JOut_full_in(1:num_curr)+A_inv*(-Ae')*...
                    (PP * (LL' \ (LL \ (PP' * (Ae*A_inv*JOut_full_in(1:num_curr))))));
                
                % block F contribution
                JOut_full_out(1:num_curr) = JOut_full_out(1:num_curr)...
                    + A_inv * (Ae') * (PP * (LL' \ (LL \ (PP' * (JOut_full_in(num_curr+1:num_curr+num_node))))));
                
                % block G contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    - (PP * (LL' \ (LL \ (PP' * (Ae*A_inv*JOut_full_in(1:num_curr))))));
                
                % block H contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = JOut_full_out(num_curr+1:num_curr+num_node)...
                    + (PP * (LL' \ (LL \ (PP' * (JOut_full_in(num_curr+1:num_curr+num_node))))));
                
            end
        end
        
        
    case 'no_decomp'
        
        
        if (fl_volt_source == 1)
            
            Ae_tmp=Ae; % Attention::: temporarily doubling the memory for Ae!!!
            Ae_tmp(nodeid_4_grnd,:)=0;
            Ae_tmp(nodeid_4_injectcurr,:)=0;
            
            if (fl_fast_multiply == 1)
                
                disp('Fast solution !!!')
                
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    (Sch_sparse \ (JOut_full_in(num_curr+1:num_curr+num_node) - ( Ae_tmp * A_inv * JOut_full_in(1:num_curr) ) ));
                
                JOut_full_out(1:num_curr) = (JOut_full_in(1:num_curr) - ((-Ae')*JOut_full_out(num_curr+1:num_curr+num_node))) .* diag(A_inv);
                
                
            else
                
                disp('Slow solution !!!')
                
                % block E contribution
                JOut_full_out(1:num_curr) = A_inv*JOut_full_in(1:num_curr)+A_inv*(-Ae')*...
                    (Sch_sparse \ (Ae_tmp*A_inv*JOut_full_in(1:num_curr)));
                
                % block F contribution
                JOut_full_out(1:num_curr) = JOut_full_out(1:num_curr)...
                    + A_inv * (Ae') * (Sch_sparse \ (JOut_full_in(num_curr+1:num_curr+num_node)));
                
                % block G contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    - Sch_sparse \ (Ae_tmp*A_inv*JOut_full_in(1:num_curr));
                
                % block H contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = JOut_full_out(num_curr+1:num_curr+num_node)...
                    + (Sch_sparse \ (JOut_full_in(num_curr+1:num_curr+num_node)));
                
            end
            
        else % current source or symmetic voltage source
            
            if (fl_fast_multiply == 1)
                
                %disp('Fast solution !!!')
                
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    (Sch_sparse \ (JOut_full_in(num_curr+1:num_curr+num_node) - ( Ae * A_inv * JOut_full_in(1:num_curr) )));
                
                JOut_full_out(1:num_curr) = (JOut_full_in(1:num_curr) - ((-Ae')*JOut_full_out(num_curr+1:num_curr+num_node))) .* diag(A_inv);
                
            else
                
                %disp('Slow solution !!!')
                
                % block E contribution
                JOut_full_out(1:num_curr) = A_inv*JOut_full_in(1:num_curr)+A_inv*(-Ae')*...
                    (Sch_sparse \ (Ae*A_inv*JOut_full_in(1:num_curr)));
                
                % block F contribution
                JOut_full_out(1:num_curr) = JOut_full_out(1:num_curr)...
                    + A_inv * (Ae') * (Sch_sparse \ (JOut_full_in(num_curr+1:num_curr+num_node)));
                
                % block G contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    - Sch_sparse \ (Ae*A_inv*JOut_full_in(1:num_curr));
                
                % block H contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = JOut_full_out(num_curr+1:num_curr+num_node)...
                    + (Sch_sparse \ (JOut_full_in(num_curr+1:num_curr+num_node)));
                
            end
        end
        
    otherwise
        
        error('Invalid decomposition selection for Schur complement')
        
end

% moving the printing of dots to the 'lse_matvect_mult.m' function
%fprintf ('.') ;
