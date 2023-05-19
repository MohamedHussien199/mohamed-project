function [x numiter]=des(A,b,h,tol)
n=length(b);
xprev=zeros(1,n)';
w=h*A'*A;
w1=h*A'*b;
N=eye(n)-w;
iter=0;
while true
    iter=iter+1;
    xnew=(N*xprev)+w1;
    if norm(xnew-xprev,2)<tol
        numiter=iter;
        x=xnew
        break
    end
    xprev=xnew;
end
if norm(xnew-xprev)>tol
    numiter=-1;
    x=xnew;
end
end

