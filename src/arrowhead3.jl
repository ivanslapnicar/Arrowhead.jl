#------------------------
#--------  Functions

#--------  random generation

function GenSymArrow(n::Integer,i::Integer)
    # generates symmetric n x n arrowhad matrix with arrow top at (i,i)
    SymArrow(rand(n-1),rand(n-1),rand(),i)
end

function GenSymDPR1(n::Integer)
    # generates symmetric n x n DPR1 matrix 
    SymDPR1(rand(n),rand(n),rand())
end

#-------- Inverses  

function inv{T}(A::SymArrow{T},i::Integer,tols::Vector{Float64})

    # COMPUTES: inverse of a SymArrow matrix A, inv(A-A.D[i]*I) which is again SymArrow
    # uses higher precision to compute top of the arrow element accurately, if
    # needed. 
    # tols=[tolb,tolz] are tolerances, usually [1e3, 10*n]
    # [0.0,0.0] forces DoubleDouble, [1e50,1e50] would never use it
    # RETURNS:  SymArrow(D,z,b,i), Kb, Kz, Qout
    # Kb - condition Kb, Kz - condition Kz, Qout = 1 / 0 - double was / was not used 

    n=length(A.D)
    D=Array{T}(n)
    z=Array{T}(n)
    wz=1/A.z[i]
    shift=A.D[i]

    for k=1:i-1
        D[k]=one(T)/(A.D[k]-shift)
        z[k]=-A.z[k]*D[k]*wz
    end
    for k=i+1:n
        D[k-1]=one(T)/(A.D[k]-shift)
        z[k-1]=-A.z[k]*D[k-1]*wz
    end

    # D=[1./(A.D[1:i-1]-shift),1./(A.D[i+1:end]-shift)]
    # wz=1/A.z[i]
    # z=-[A.z[1:i-1], A.z[i+1:end]].*D1*wz1
    a=A.a-A.D[i]

    # compute the sum in a plain loop
    P=zero(T)
    Q=zero(T)
    Qout=0

    nn=n-1
    for k=1:i-1
        D[k]>0.0 ? P+=A.z[k]^2*D[k] : Q+=A.z[k]^2*D[k]
    end
    for k=i:nn
        D[k]>0.0 ? P+=A.z[k+1]^2*D[k] : Q+=A.z[k+1]^2*D[k]
    end
    a<0 ? P=P-a : Q=Q-a

    Kb=(P-Q)/abs(P+Q)
    Kz=maximum(abs,A.z)*abs(wz)

    if Kb<tols[1] ||  Kz<tols[2]
        b=(P+Q)*wz*wz
    else  # recompute in Double or BigFloat
        if Kz<1.0/eps()^2 
            Type=Double
            Qout=1
        else
        # Example of a matrix where BigFloat is neeed, courtesy of Stan Eisenstat, is:
        # A=SymArrow([1+eps(), 1-eps(), 0,-1+2*eps(), -1-2*eps()],[2,2,eps()^4,1,1.0],6.0,5)
            Type=BigFloat
            Qout=50
        end
        shiftd=map(Type,A.D[i])
        Dd=[[Type(A.D[k])-shiftd for k=1:i-1]; 
            [Type(A.D[k])-shiftd for
                            k=i+1:length(A.D)]]
        wzd=Type(A.z[i])
        ad=Type(A.a)-shiftd

        Pd,Qd=map(Type,(0.0,0.0))
        
        for k=1:i-1
            map(Float64,Dd[k])>0.0 ? Pd+=Type(A.z[k])^2/Dd[k] :
            Qd+=Type(A.z[k])^2/Dd[k]
        end
        
        for k=i:nn
            map(Float64,Dd[k])>0.0 ? Pd+=Type(A.z[k+1])^2/Dd[k] :
            Qd+=Type(A.z[k+1])^2/Dd[k]
            # @show P,Q
        end 

        map(Float64,ad)<0 ?   Pd=Pd-ad : Qd=Qd-ad 

        bd=(Pd+Qd)/(wzd*wzd)
        b=map(Float64,bd)
# @show b, b-(bd.hi+bd.lo)
    end

    if i<A.i 
        # SymArrow([D[1:A.i-2],zero(T),D[A.i-1:end]],[z[1:A.i-2],wz,z[A.i-1:end]],b,i),Kb,Kz,Qout
        # This reduces memory allocation (4 vectors) and is 30-40% faster
        for k=n-1:-1:A.i-1
            D[k+1]=D[k]
            z[k+1]=z[k]
        end
        D[A.i-1]=zero(T)
        z[A.i-1]=wz   
        SymArrow(D,z,b,i),Kb,Kz,Qout
    else
        # SymArrow([D[1:A.i-1],zero(T),D[A.i:end]],[z[1:A.i-1],wz,z[A.i:end]],b,i+1),Kb,Kz,Qout
        for k=n-1:-1:A.i
            D[k+1]=D[k]
            z[k+1]=z[k]
        end
        D[A.i]=zero(T)
        z[A.i]=wz   
        SymArrow(D,z,b,i+1),Kb,Kz,Qout
    end

end # inv


function inv{T}(A::SymArrow{T}, shift::Float64, tolr::Float64)

    # COMPUTES: inverse of the shifted SymArrow A, inv(A-shift*I) which is SymDPR1
    # uses DoubleDouble to compute rho accurately, if needed. 
    # tolr is tolerance, usually 1e3,  0.0 forces Double, 1e50 would never use it
    # RETURNS: SymDPR1(D,u,rho), Krho, Qout
    # Krho - condition Krho, Qout = 1 / 0 - Double was / was not used 

    n=length(A.D)
    D=Array{T}(n+1)
    u=Array{T}(n+1)

    # Ds=A.D-shift

    # D=[1./Ds[1:A.i-1];0.0;1./Ds[A.i:end]];
    # u=[A.z[1:A.i-1]./Ds[1:A.i-1];-1.0;A.z[A.i:end]./Ds[A.i:end]];

    for k=1:A.i-1
        D[k]=one(T)/(A.D[k]-shift)
        u[k]=A.z[k]*D[k]
    end
    D[A.i]=zero(T)
    u[A.i]=-one(T)
    for k=A.i:n
        D[k+1]=one(T)/(A.D[k]-shift)
        u[k+1]=A.z[k]*D[k+1]
    end

    # compute rho and Krho
    #--- compute the sum in a plain loop

    a=A.a-shift
    P=zero(T)
    Q=zero(T)
    Qout=0

    for k=1:A.i-1
        D[k]>0.0 ? P=P+A.z[k]^2*D[k] : Q=Q+A.z[k]^2*D[k]
    end
    for k=A.i:n
        D[k+1]>0.0 ? P=P+A.z[k]^2*D[k+1] : Q=Q+A.z[k]^2*D[k+1]
    end

    a<0 ? P=P-a : Q=Q-a

    # Condition of rho
    Krho=(P-Q)/abs(P+Q)


    if Krho<tolr
        rho=-1.0/(P+Q)

    else  # recompute in Double
        Qout=1
        Pd,Qd=map(Double,(0.0,0.0))
        shiftd=Double(shift)
        ad=Double(A.a)-shiftd
        #   Dd=[Double(A.D[k])-shiftd for k=1:length(A.D) ] 
        #   Dd=map(Double,A.D)-shiftd
        for k=1:A.i-1
            D[k]>0.0 ? Pd=Pd+Double(A.z[k])^2/(Double(A.D[k])-shiftd) : Qd=Qd+Double(A.z[k])^2/(Double(A.D[k])-shiftd)
        end
        for k=A.i:n
            D[k+1]>0.0 ? Pd=Pd+Double(A.z[k])^2/(Double(A.D[k])-shiftd) : Qd=Qd+Double(A.z[k])^2/(Double(A.D[k])-shiftd)
        end

        # for k=1:length(A.D)
        #    Dd[k].hi>0.0 ? Pd=Pd+Double(A.z[k])^2/Dd[k] : Qd=Qd+Double(A.z[k])^2/Dd[k]
        # end
        ad.hi+ad.lo<0 ? Pd=Pd-ad : Qd=Qd-ad
        if Pd+Qd!=Double(0.0,0.0)
            Krho=Float64((Pd-Qd)/abs(Pd+Qd))
            r=Double(1.0)/(Pd+Qd)
            rho=-(r.hi+r.lo)
        else
        # Here we need quadruple working precision. We are using BigFloat.
        # Example of a matrix where this is neeed, courtesy of Stan Eisenstat, is:
        # A=SymArrow([1+eps(), 1-eps(), -1+2*eps(), -1-2*eps()],[2,2,1,1.0],6.0,5) 
            Qout=100
            Pd,Qd=map(BigFloat,(0.0,0.0))
            shiftd=BigFloat(shift)
            ad=BigFloat(A.a)-shiftd
            for k=1:A.i-1
                D[k]>0.0 ? Pd=Pd+BigFloat(A.z[k])^2/(BigFloat(A.D[k])-shiftd) : Qd=Qd+BigFloat(A.z[k])^2/(BigFloat(A.D[k])-shiftd)
            end
            for k=A.i:n
                D[k]>0.0 ? Pd=Pd+BigFloat(A.z[k])^2/(BigFloat(A.D[k])-shiftd) : Qd=Qd+BigFloat(A.z[k])^2/(BigFloat(A.D[k])-shiftd)
            end
            ad<0 ? Pd=Pd-ad : Qd=Qd-ad
            Krho=Float64((Pd-Qd)/abs(Pd+Qd))
            r=BigFloat(-1.0)/(Pd+Qd)
            rho=Float64(r)
        end
    end

    # returns the following
    SymDPR1(D,u,rho), Krho, Qout

end # inv

function inv{T}(A::SymArrow{T}, shift::Double)

    # COMPUTES: inverse of the shifted SymArrow A, inv(A-shift*I), which is a SymDPR1
    # here shift is Double so it uses Double to compute everything 
    # RETURNS: SymDPR1(D1,u1,rho), Qout
    # Qout = 1 on exit meaning Double was used 

    n=length(A.D)
    D=Array(Double,n+1)
    u=Array(Double,n+1)

    # D[1:n]=map(Double,A.D)-shift
    for k=1:n
        D[k]=Double(A.D[k])-shift
    end

    u[1:n]=map(Double,A.z)
    a=map(Double,A.a)-shift
    i=A.i

    oned=Double(1.0)
    zerod=Double(0.0)


    # Ds=A.D-shift
    # D=[1./Ds[1:A.i-1];0.0;1./Ds[A.i:end]];
    # u=[A.z[1:A.i-1]./Ds[1:A.i-1];-1.0;A.z[A.i:end]./Ds[A.i:end]];

    for k=1:i-1
        #   tmp=Double(A.D[k])-shift
        #    D[k]=oned/(Double(A.D[k])-shift)
        u[k]=u[k]/D[k]  
        D[k]=oned/D[k]
        #   u[k]=Double(A.z[k])/tmp # *D[k]
    end

    for k=i:n
        u[k+1]=u[k]/D[k]
        D[k+1]=oned/D[k+1]  
        # D[k+1]=oned/(Double(A.D[k])-shift)
        # u[k+1]=Double(A.z[k])*D[k+1]
    end

    D[i]=zerod
    u[i]=Double(-1.0)

    # D=[Double(1.0)./D1[1:A.i-1];Double(0.0);Double(1.0)./D1[A.i:end]];
    # u=[z1[1:A.i-1]./D1[1:A.i-1];Double(-1.0);z1[A.i:end]./D1[A.i:end]];

    # compute rho and Krho
    #--- compute the sum in a plain loop

    P,Q=zerod,zerod
    Qout=1

    for k=1:i-1
        # D[k].hi+D[k].lo > 0.0  # .lo cannot change sign, I think
        # D[k].hi > 0.0 ? P=P+Double(A.z[k])^2*D[k] : Q=Q+Double(A.z[k])^2*D[k]
        D[k].hi > 0.0 ? P=P+Double(A.z[k])*u[k] : Q=Q+Double(A.z[k])*u[k]
    end
    for k=i:n
        D[k+1].hi >0.0 ? P=P+Double(A.z[k])*u[k+1] : Q=Q+Double(A.z[k])*u[k+1]
    end

    a.hi+a.lo<0.0  ? P=P-a : Q=Q-a
    r=oned/(P+Q)

    # for k=1:length(A.D)
    #    D1[k].hi+D1[k].lo>0.0 ? P=P+z1[k]^2/D1[k] : Q=Q+z1[k]^2/D1[k]
    # end

    rho=-(r.hi+r.lo)

    # returns the following
    D1=Array{T}(n+1)
    u1=Array{T}(n+1)
    for k=1:n+1
        D1[k]=D[k].hi+D[k].lo
    end
    for k=1:n+1
        u1[k]=u[k].hi+u[k].lo
    end

    # SymDPR1(T[x.hi+x.lo for x=D],T[x.hi+x.lo for x=u],rho), Qout

    SymDPR1(D1,u1,rho), Qout

end # inv


function  eig{T}( A::SymArrow{T},k::Integer,tols::Vector{Float64})

    # COMPUTES: k-th eigenpair of an ordered irreducible SymArrow
    # A = [diag (D) z; z' alpha]
    # tols=[tolb,tolz,tolnu,tolrho,tollambda] = [1e3,10.0*n,1e3,1e3,1e3]
    # RETURNS: lambda, v, Sind, Kb, Kz, Knu, Krho, Qout
    # lambda - k-th eigenvalue in descending order
    # v - lambda's normalized eigenvector
    # Kb, Kz, Knu, Krho - condition numbers
    # Qout = 1 / 0 - Double was / was not used 

    # If Qout>0, quad precision was used
    # i was the shift index

    # Set the dimension
    n = length(A.D) + 1

    # Set all conditions initially to zero
    Kb,Kz,Knu,Krho=0.0,0.0,0.0,0.0
    Qout=0
    v=zeros(T,n)
    # Kz is former kappa_nu

    # Determine the shift sigma, the shift index i, and whether lambda 
    # is on the left or the right side of the nearest pole
    # Exterior eigenvalues (k = 1 or k = n):

    if k == 1 
        sigma,i,side = A.D[1],1,'R' 
    elseif k==n 
        sigma,i,side = A.D[n-1],n-1,'L'
    else
        # Interior eigenvalues (k in (2,...n-1) ):
        Dtemp = A.D-A.D[k]
        atemp = A.a-A.D[k]
        middle = Dtemp[k-1]/2.0
        Fmiddle = (atemp-middle)-sum(A.z.^2 ./(Dtemp-middle))
        # middle=(A.D[k-1]-A.D[k])/2.0
        # Fmiddle=0.0
        # for l=1:n-1
        #     Fmiddle+=A.z[l]^2/((A.D[l]-A.D[k])-middle)
        # end
        # Fmiddle=((A.a-A.D[k])-middle)-Fmiddle
        # @show Fmiddle
        sigma,i,side = Fmiddle < 0.0 ? (A.D[k],k,'R') : (A.D[k-1],k-1,'L')
    end

    # Compute the inverse of the shifted matrix, A_i^(-1), Kb and Kz

    Ainv,Kb,Kz,Qout = inv(A,i,tols[1:2])

    # Compute the eigenvalue of the inverse shifted matrix

    nu = bisect( Ainv,side )

    #  nu=fastzero([invD1; 0; invD2], [w1;wz;w2], b, side); # implement later
    #  [nu-nu1]/nu, [nu-nueva]/nu, pause  # and compare

    if abs(nu)==Inf
        # this is nonstandard
        # Deflation in aheig (nu=Inf)
        v[i]=1.0
        lambda=sigma
    else
        # standard case, full computation
        # nu1 is the F- or 1-norm of the inverse of the shifted matrix
        # nu10=maximum([sum(abs.(Ainv.z))+abs(Ainv.a), maximum(abs,Ainv.D)+abs.(Ainv.z))])
        nu1=0.0
        for l=1:n-1
            nu1=max(nu1,abs(Ainv.D[l])+abs(Ainv.z[l]))
        end
        nu1=max(sum(abs,Ainv.z)+abs(Ainv.a), nu1)

        Knu=nu1/abs(nu)

        while  Knu>tols[3]
            # Remedies according to Remark 3 - we shift between original
            # eigenvalues and compute DPR1 matrix
            # 1/nu1+sigma, 1/nu+sigma
            # println("Remedy 3 ")
            nu = side=='R' ? abs(nu) : -abs(nu)
            nu1=-sign(nu)*nu1
            sigma1=(nu1+nu)/(2.0*nu*nu1)+sigma

            if findfirst(A.D-sigma1,0.0)>0 # we came back with a pole
                # recompute sigmav more accurately with Dekker
                sigmav=(Double(nu1)+Double(nu))/(Double(2.0)*Double(nu)*Double(nu1))+Double(sigma)
                # Compute the inverse of the shifted arrowhead (DPR1)
                Ainv,Qout1=inv(A,sigmav) # Ainv is Float64
                nu1=bisect(Ainv,side) 
                mu1 = 1.0/nu1
                D0=map(Double,A.D)-sigmav
                D1=[D0[l].hi+D0[l].lo for l=1:length(D0)]
                if findfirst(D1-mu1,0.0)>0 
                    ind=find(D1-mu1==0.0)
                    v=zeros(n)
                    v[ind]=1.0
                else
                    v=[ A.z./(D1-mu1);-1.0]
                end
                # Shift the eigenvalue back in Double
                lam = Double(1.0)/Double(nu1)+sigmav
                # Return this
                lambda,v = lam.hi+lam.lo, v/norm(v)
                return lambda,v,i,Kb,Kz,Knu,Krho,Qout
            else
                # Compute the inverse of the shifted arrowhead (DPR1)
                Ainv, Krho,Qout1=inv(A,sigma1,tols[4]) # Ainv is Float64
                # Compute the eigenvalue by bisect for DPR1
	        # Note: instead of bisection could use dlaed4 (need a wrapper) but
	        # it is not faster. There norm(u)==1
                nu= bisect(Ainv,side)
                if side=='R' && A.D[1]>0.0 && nu<0.0
                    nu=bisect(Ainv,'L')
                end
                mu=1.0/nu 
                nu1=maximum(abs,Ainv.D)+abs(Ainv.r)*dot(Ainv.u,Ainv.u)
                Knu=nu1/abs(nu)
                sigma=sigma1
                # standard v
                # v=[ A.z./((A.D-sigma1)-mu1);-1.0]
                # Return this - shift the eigenvalue back and normalize the vector
                # lambda, v = mu1+sigma1, v/norm(v)
            end
            Qout=Qout+2*Qout1
        end

        # Accuracy is fine, compute the eigenvector
        mu = 1.0/nu
        # v=[ A.z./((A.D-sigma)-mu);-1.0]
        for l=1:n-1
            v[l] = A.z[l]/((A.D[l]-sigma)-mu)
        end 
        v[n]=-1.0
        lambda, v = mu+sigma, v/norm(v)

        # Remedy according to Remark 1 - we recompute the the eigenvalue
        # near zero from the inverse of the original matrix (a DPR1 matrix).  
        if (abs(A.D[i])+abs(1.0/nu))/abs(lambda)>tols[5]

            if k==1 && A.D[1]<0.0 || k==n && A.D[n-1]>0.0 || i<n-1 && side=='L' && sign(A.D[i])+sign(A.D[i+1])==0 || i>1
                side=='R' && sign(A.D[i])+sign(A.D[i-1])==0
                # println("Remedy 1 ")
                # Compute the inverse of the original arrowhead (DPR1)
                Ainv,Krho,Qout1 = inv(A,0.0,tols[4]) # Ainv is Float64
                Qout=Qout+4*Qout1
                if abs(Ainv.r)==Inf
                    lambda=0.0
                else                
                    # Here we do not need bisection. We compute the Rayleigh
                    # quotient by using already computed vectors which is
                    # componentwise accurate
                    nu1=sum(v.^2 .*Ainv.D)+Ainv.r*sum(v.*Ainv.u)^2;
                    lambda=1.0/nu1
                end
            end

        end
    end

    # Return this
    lambda,v,i,Kb,Kz,Knu,Krho,Qout

end # eig (k)


function eig{T}(A::SymArrow{T}, tols::Vector{Float64})

    # COMPUTES: all eigenvalues and eigenvectors of a real symmetric SymArrow 
    # A = [diag(D) z;z' alpha] (notice, here we assume A.i==n)
    # tols = [tolb,tolz,tolnu,tolrho,tollambda] = [1e3,10.0*n,1e3,1e3,1e3] or similar
    # RETURNS: U, E, Sind, Kb, Kz, Knu, Krho, Qout
    # U = eigenvectors, E = eigenvalues in decreasing order 
    # Sind[k] - shift index i for the k-th eigenvalue
    # Kb, Kz, Knu, Krho [k] - respective conditions for the k-th eigenvalue
    # Qout[k] = 1 / 0 - Double was / was not used when computing k-th eigenvalue 
    
    n=length(A.D)+1
    n0=n

    # Ordering the matrix
    is=sortperm(A.D,rev=true)
    D=A.D[is]
    z=A.z[is]
    U=eye(T,n,n)
    E=zeros(T,n)
    Kb=zeros(T,n); Kz=zeros(T,n); Knu=zeros(T,n); Krho=zeros(T,n)
    Qout=zeros(Int,n); Sind=zeros(Int,n) 

    # Quick return for 1x1
    if n==1
        return U,[A.a],Sind,Kb,Kz,Knu,Krho,Qout
    end

    #  test for deflation in z
    z0=find(z.==0)
    zx=find(z.!=0)
    
    if isempty(zx)  # nothing to do
        E=[A.D;A.a]
        isE=sortperm(E,rev=true)
        E=E[isE]
        U=U[:,isE]
        return U,E,Sind,Kb,Kz,Knu,Krho,Qout
    end
        
    if !isempty(z0)
        E[z0]=D[z0]
        D=D[zx]
        z=z[zx]
        if isempty(z)
            E[n]=A.a
            # return U,E,Sind,Kb,Kz,Knu,Krho,Qout
        else
            n=length(z)+1
        end
    end
    zx=[zx;n0]
    
    #  Test for deflation in D
    g=D[1:n-2]-D[2:n-1]
    # Can play with inexact deflation
    # g0=find(abs(g)<eps)
    # gx=find(abs(g)>=eps)
    # Only exact deflation !!
    g0=find(g.==0.0)
    gx=find(g.!=0.0)
    if !isempty(g0)
        # Deflation
        Dgx=D[gx]; zgx=z[gx]
        lg0=length(g0)
        R=Array(Tuple{Givens{Float64},Float64},lg0)
        for l=lg0:-1:1
            # This will be changed in v0.4
            R[l]=givens(z[g0[l]],z[g0[l]+1],zx[g0[l]],zx[g0[l]+1])
            z[g0[l]]=R[l][2]; z[g0[l]+1]=0.0
            # A_mul_Bc!(U,R) # multiply R'*U later
            E[zx[g0[l]+1]]=D[g0[l]+1]
        end    
        # remains
        gx=[0;gx]+1
        nn=length(gx)
        
        zxx=zx[[gx;n]]
        for k=1:nn+1
            E[zxx[k]],U[zxx,zxx[k]],Sind[zxx[k]],Kb[zxx[k]],Kz[zxx[k]],Knu[zxx[k]],Krho[zxx[k]],Qout[zxx[k]]=
            eig(SymArrow(D[gx],z[gx],A.a,nn+1),k,tols)
        end
        
        for l=1:lg0
            U=R[l][1]'*U
        end 
        
    else
        
        # No deflation in D
        for k=1:n
            E[zx[k]],U[zx,zx[k]],Sind[zx[k]],Kb[zx[k]],Kz[zx[k]],Knu[zx[k]],Krho[zx[k]],Qout[zx[k]]=
            eig(SymArrow(D,z,A.a,n),k,tols)
        end
    end
    # end
    # back premutation of vectors
    isi=sortperm(is)
    
    # must sort E once more 
    es=sortperm(E,rev=true)
    E=E[es]
    
    U=U[[isi[1:A.i-1];n0;isi[A.i:n0-1]],es]
    
    # Return this
    U,E,Sind[es],Kb[es],Kz[es],Knu[es],Krho[es],Qout[es]
    
end # eig (all)


function bisect{T}(A::SymArrow{T}, side::Char)
    # COMPUTES: the leftmost (for side='L') or the rightmost (for side='R') eigenvalue
    # of a SymArrow A = [diag (D) z; z'] by bisection.
    # RETURNS: the eigenvalue
    
    # Determine the starting interval for bisection, [left; right]
    # left, right = side == 'L' ? {minimum([A.D-abs(A.z),A.a-sum(abs(A.z))]), minimum(A.D)} : 
    #   {maximum(A.D),maximum([A.D+abs.(A.z),A.a+sum(abs,A.z)])}
    
    absAz = abs.(A.z)
    if side == 'L'
        left  = minimum(A.D - absAz)::T
        left  = min(left, A.a - sum(absAz))::T
        right = minimum(A.D)::T
    else
        left  = maximum(A.D)::T
        right = maximum(A.D + absAz)::T
        right = max(right, A.a + sum(absAz))::T
    end
    
    # Bisection
    middle = (left + right) / convert(T,2)
    z2 = A.z .^ 2
    count, n = 0, length(A.D)
    
    while (right-left) > 2.0 * eps() * max(abs(left), abs(right))
        # in @time 50% of the time was garbage collection. The following line
        # assigns new vector every time it is called, so it is much better in the
        # loop?? Original timing were 30 secs for n=4000, 2.17 sec for n=1000 
        
        # Fmiddle = A.a-middle-sum(z2./(A.D-middle))
        
        Fmiddle = zero(T)
        for k=1:n
            Fmiddle += z2[k] / (A.D[k] - middle)
        end
        Fmiddle = A.a - middle - Fmiddle
        
        if Fmiddle > zero(T)
            left = middle
        else
            right = middle
        end
        middle = (left + right) / convert(T,2)
    end
    # Return the eigenvalue
    right
    
end # bisect

    function bisect( A::SymDPR1, side::Char )
        # COMPUTES: the leftmost (for side='L') or the rightmost (for side='R') eigenvalue
        # of a SymDPR1 matrix A = diagm(A.D) + A.r*A.u*(A.u)' by bisection.
        # RETURNS: the eigenvalue

        n=length(A.D)
        # Determine the starting interval for bisection, [left; right]
        indD=sortperm(A.D,rev=true)

        if A.r>0.0
            left, right = side == 'L' ? (A.D[indD[n]], A.D[indD[n-1]]) :  (A.D[indD[1]], A.D[indD[1]]+A.r*dot(A.u,A.u))
        else # rho<=0
            left, right = side == 'L' ? (A.D[indD[n]]+A.r*dot(A.u,A.u), A.D[indD[n]]) : (A.D[indD[2]], A.D[indD[1]])
        end

        # Bisection
        middle = (left + right)/2.0
        u2=A.u.^2
        n = length(A.D)
        while (right-left) > 2.0*eps()*maximum([abs(left),abs(right)])

            # Fmiddle = 1.0+A.r*sum(u2./(A.D-middle))

            Fmiddle=0.0
            for k=1:n
                Fmiddle=Fmiddle+u2[k]/(A.D[k]-middle)
            end    
            Fmiddle=1.0+A.r*Fmiddle

            sign(A.r)*Fmiddle < 0.0 ? left = middle : right = middle
            middle = (left + right)/2.0
        end
        # Return the eigenvalue
        right

    end # bisect

    # this is just for convenience
    # tols=[1e3,1000.0,1e3,1e3,1e3]
