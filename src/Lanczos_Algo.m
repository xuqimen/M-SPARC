%LANCZOS Algorithm to find out the min & max eigen values of a real symmetric matrix A.
function [Lmin,Lmax]=Lanczos_Algo(A,tol)
    N=length(A);
    vkm1 = rand(N,1) ;      % Arbitrary initial guess for LANCZOS vector
    vkm1 = vkm1/norm(vkm1); % Normalize it, this is q1 in the algorithm
    vk = A*vkm1 ;           % The first v is A*q (from algorithm)
    a(1) = vkm1'*vk ;       % alpha=q'*v   
    vk = vk - a(1)*vkm1 ;   %update v=v-beta0*q0-alpha1*q1 (beta0, q0 = 0)
    b(1) = norm(vk) ;       %beta=b=norm(vk), from algorithm
    vk = vk/b(1) ;          %update vk as new qi (to be used in next iteration) by normalizing
   
    k=1;
    L(1,1)=0;
    M(1,1)=0;
    DL=1;
    DM=1;
    while (DL>tol) || (DM>tol)      %run iterations until eigen values converge to within tol
        vkp1 = A*vk ;
        a(k+1) = transpose(vk)*vkp1  ;
        vkp1 = vkp1 - a(k+1)*vk - b(k)*vkm1 ; %same as vk=A*vk-a(k+1)*vk-b(k)*vkm1;
        vkm1 = vk ;               %this will be q(n-1) for next iteration
        b(k+1) = norm(vkp1) ;  %=norm(of current iteration's updated vk)
        vk = vkp1/b(k+1) ;
        s=length(b);
        J = diag(a) + diag(b(1:s-1),1) + diag(b(1:s-1),-1) ;
        [Evec,Eval] = eig(J) ;
        L(1,k+1)=min(diag(Eval)); %eigs(J,1,'la');
        M(1,k+1)=max(diag(Eval));
        DL=abs(L(1,k+1)-L(1,k));
        DM=abs(M(1,k+1)-M(1,k));
        k=k+1;
    end;
fprintf('Lanczos_Algo takes %d steps\n',k-1);
Lmin=L(end);
Lmax=M(end);
end