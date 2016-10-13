function f=DIVISION(A,B)
f=A./B;
f(isnan(f))=0;
f(isinf(f))=0; 
