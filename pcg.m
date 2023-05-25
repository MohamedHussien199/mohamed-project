function[x numiter]=pcg(A,b,tol,maxiter)
n=length(b); %size of b.
xprev=zeros(1,n)'; %the inital guess.
L=ichol(A); % incomplete decompsition of A.
rprev=b-(A*xprev); %thre initial residual.
zprev=L\rprev; %preconditioned initial residual.
dprev=L'\zprev; % preconditioned initial direction.
for k=1:maxiter
    % stopping criteria.
    if norm(dprev,2)<tol
       x=xprev;
        numiter=k;
        break
    end
    u=A*dprev;
    %the step size for the next approximation solution.
    alpha=(norm(zprev,2))^2/(dprev'*u);
    % the next approximation solution.
     xnew=xprev+(alpha*dprev);
     % the next residual.
    rnew=rprev-(alpha*u);
    %the next preconditioned residual.
    znew=L\rnew;
    % stopping criteria.
    if (norm(znew,2))^2<tol
        if norm(rnew,2)<tol
            x=xnew;
            numiter=k;
            break
        end
    end
    %the step size for the next direction.
    beta=(norm(znew,2))^2/(norm(zprev,2))^2;
    %the next preconditioned direction.
    dnew=(L'\znew)+(beta*dprev);
    % updating all the variable.
    zprev=znew;
    dprev=dnew;
    xprev=xnew;
    rprev=rnew;
end
% in case of divergence.
if norm(rnew,2)>=tol
    x=xnew;
    numiter=-1;
end
end


