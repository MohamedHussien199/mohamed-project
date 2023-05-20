function [x numiter]=des(A,b,h,tol)
n=length(b);% size of b.
xprev=zeros(1,n)';% the initial guess.
w=h*A'*A; % to construct the iteration matrix.
w1=h*A'*b;% the iteration vector.
N=eye(n)-w; % the iteration matrix.
iter=0; % to count the number of iterations.
while true
    iter=iter+1;
    xnew=(N*xprev)+w1;
    % stopping criteria.
    if norm(xnew-xprev,2)<tol
        numiter=iter;
        x=xnew;
        break
    end
    xprev=xnew;
end
% in case of divergence.
if norm(xnew-xprev)>=tol
    numiter=-1;
    x=xnew;
end
end

