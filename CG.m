function [x  numiter]=CG(A,b,tol,maxiter)
n=length(b);%size of b.
xprev=zeros(1,n)';% the initial guess.
rprev=b-(A*xprev); %the initial resuidal.
dprev=rprev; %the inital direction.
    for i=1:maxiter
        %calculatig the step size.
        alpha=(norm(rprev,2))^2/(dprev'*A*dprev);
        % the next approximation solution.
        xnew=xprev+(alpha*dprev);
        %the next residual.
        rnew=rprev-(alpha*A*dprev);
        %stopping criteria.
         if norm(rnew,2)<tol
            numiter=i;
            x=xnew;
            break
         end
         %the step size for the next direction.
          beta=(norm(rnew,2))^2/(norm(rprev,2))^2 ;
          % the next direction.det
          dnew=rnew+(beta*dprev);
          % updating the variables.
          xprev=xnew;
          rprev=rnew;
          dprev=dnew;
    end
    %in case of divergence.
    if norm(rnew,2)>=tol
        numiter=-1;
        x=xnew;
    end
end
        