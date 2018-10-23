function [schur_out]=lse_sparse_schur_compl_multiply(schur_in, Ae, A_inv)

schur_out = Ae * A_inv * Ae' * schur_in;

fprintf ('+') ;
