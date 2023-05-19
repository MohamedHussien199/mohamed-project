function [x numiter Rerr]=jacobi(A,b,tol)
n=length(b); % size of b.
x=zeros(1,n); % inital guess.
iter=0;  % for counting the number of iterations.
while true
    iter=iter+1;
    % finding the optimal solution.
    for i=1:n
        y(i)=b(i);
        for j=1:i-1
            y(i)=y(i)-A(i,j)*x(j);
        end
        for j=i+1:n
            y(i)=y(i)-A(i,j)*x(j);
        end
        y(i)=y(i)/A(i,i);
    end
    % stopping criteria.
    if norm(x-y,2)<=tol
    Rerr=norm(x-y,2)/norm(x,2);
    x=y;
    numiter=iter; 
    break
    end
    % in case of divergence.
    if norm(x-y,2)> tol
    Rerr=norm(x-y,2)/norm(x,2);
    x=y;
    numiter=-1;
    end
end

            