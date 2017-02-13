function f=DIVISION(A,B)
f=A./B;
f(isnan(double(f)))=0;
f(isinf(double(f)))=0; 
end