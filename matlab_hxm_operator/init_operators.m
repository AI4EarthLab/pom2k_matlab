function ierr=init_operators(m,n,k)
    global OP
    OP=Operator(m,n,k);
    ierr=true;
end