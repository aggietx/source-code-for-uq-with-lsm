function b1=form_sigmatrix(Dim_UB,sigmak2)
%%% only for shallow water array
%%% form conjuate sigma matrix
b1 = zeros(Dim_UB, Dim_UB);

for i = 1: 2:  Dim_UB-1
        
        b1(i,i) = sigmak2(i);
        b1(i+1, i+1) = -1i* sigmak2(i);
        b1(i, i+1) = 1i * sigmak2(i);
        b1(i+1, i) = 1  * sigmak2(i);
   
end