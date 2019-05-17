#-------- Inverses of DPR1

function inv(A::SymDPR1{T},i::Integer,τ::Vector{Float64}=[1e3,10.0*length(A.D)]) where T
Kρ
    # COMPUTES: inverse of a σed SymDPR1 matrix A=diagm(A.D)+A.r*A.u*A.u',
    # inv(A-A.D[i]*I) which is a SymArrow.
    # Uses higher precision to compute top of the arrow element accurately, if
    # needed.
    # τ=[tolb,tolz] are tolerances, usually [1e3, 10*n]
    # [0.0,0.0] forces DoubleDouble, [1e50,1e50] would never use it
    # RETURNS:  SymArrow(D,z,b,i), Kb, Kz, Qout
    # Kb - condition Kb, Kz - condition Kz, Qout = 1 / 0 - double was / was not used

    n=length(A.D)
    D=Array{T}(n-1)
    z=Array{T}(n-1)
    wz=one(T)/A.u[i]
    σ=A.D[i]

    for k=1:i-1
        D[k]=one(T)/(A.D[k]-σ)
        z[k]=-A.u[k]*D[k]*wz
    end
    for k=i+1:n
        D[k-1]=one(T)/(A.D[k]-σ)
        z[k-1]=-A.u[k]*D[k-1]*wz
    end

    # compute the sum in a plain loop
    P=zero(T)
    Q=zero(T)
    Qout=0

    for k=1:i-1
        D[k]>0.0 ? P+=A.u[k]^2*D[k] : Q+=A.u[k]^2*D[k]
    end
    for k=i:n-1
        D[k]>0.0 ? P+=A.u[k+1]^2*D[k] : Q+=A.u[k+1]^2*D[k]
    end
    A.r>0 ? P=P+one(T)/A.r : Q=Q+one(T)/A.r

    Kb=(P-Q)/abs(P+Q)

    # compute Kz

    Kz=(sum(abs,A.u)-abs(A.u[i]))*abs(wz)


    # Kz=max(abs,A.z)*abs(wz)

    if Kb<τ[1] ||  Kz<τ[2]
        b=(P+Q)*wz*wz
    else  # recompute in Double or BigFoat
        if Kz<1.0/eps()^2
            Type=Double
            Qout=1
        else
        # Example of a matrix where BigFloat is neeed, courtesy of Stan Eisenstat, is:
        # A=SymDPR1([1+3*eps(), 1-3*eps(), 0, -(1-eps()), -(1+eps())],[1,1,eps()^4, 3,3],1/16)
            Type=BigFloat
            Qout=50
        end
        σ₁=map(Type,A.D[i])
        Dd=[[Type(A.D[k])-σ₁ for k=1:i-1];
            [Type(A.D[k])-σ₁ for k=i+1:length(A.D)]]
        wzd=Type(A.u[i])

        Pd,Qd=map(Type,(0.0,0.0))

        for k=1:i-1
            convert(Float64,Dd[k])>0.0 ? Pd+=Type(A.u[k])^2/Dd[k] :
            Qd+=Type(A.u[k])^2/Dd[k]
        end

        for k=i+1:n
            convert(Float64,Dd[k-1])>0.0 ? Pd+=Type(A.u[k])^2/Dd[k-1] :
            Qd+=Type(A.u[k])^2/Dd[k-1]
            # @show P,Q
        end

        A.r > 0 ?   Pd+=Type(1.0)/Type(A.r)  : Qd+=Type(1.0)/Type(A.r)

        bd=(Pd+Qd)/(wzd*wzd)
        b=convert(Float64,bd)
    end
    # return this
    SymArrow(D,z,b,i),Kb,Kz,Qout
end # inv

function inv(A::SymDPR1{T}, σ::Float64, tolr::Float64=1.0e3) where T
    # COMPUTES: inverse of the σed SymDPR1 A = diagm(A.D)+A.r*A.u*A.u',
    # inv(A-σ*I) = D + ρ*u*u', σ!=A.D[i], which is again a SymDPR1
    # uses DoubleDouble to compute A.r accurately, if needed.
    # τ is tolerance, usually 1e3,  0.0 forces Double, 1e50 would never use it
    # RETURNS: SymDPR1(D,u,ρ), Kρ, Qout
    # Kρ - condition Kρ, Qout = 1 / 0 - Double was / was not used

    n=length(A.D)
    D=Array{T}(n)
    u=Array{T}(n)

    for k=1:n
        D[k]=one(T)/(A.D[k]-σ)
        u[k]=A.u[k]*D[k]
    end

    # compute gamma and Kgamma
    #--- compute the sum in a plain loop

    P=zero(T)
    Q=zero(T)
    Qout=0

    for k=1:n
        D[k]>0.0 ? P=P+A.u[k]^2*D[k] : Q=Q+A.u[k]^2*D[k]
    end

    A.r>0 ? P=P+one(T)/A.r : Q=Q+one(T)/A.r

    # Condition of ρ
    Kρ=(P-Q)/abs(P+Q)

    if Kρ < τ
        ρ=-one(T)/(P+Q)

    else  # recompute in Double
        Qout=1
        Pd,Qd=map(Double,(0.0,0.0))
        σ₁=Double(σ)

        for k=1:n
            D[k]>0.0 ? Pd=Pd+Double(A.u[k])^2/(Double(A.D[k])-σ₁) : Qd=Qd+Double(A.u[k])^2/(Double(A.D[k])-σ₁)
        end

        A.r > 0 ?   Pd+=Double(1.0)/Double(A.r)  : Qd+=Double(1.0)/Double(A.r)
        Kρ=Float64((Pd-Qd)/abs(Pd+Qd))

        if Kρ<τ/eps()
            r=Double(1.0)/(Pd+Qd)
            ρ=-(r.hi+r.lo)
        else
        # Here we need quadruple working precision. We are using BigFloat.
        # Example of a matrix where this is neeed, courtesy of Stan Eisenstat, is:
        # A=SymDPR1([1+3*eps(), 1-3*eps(), -(1-eps()), -(1+eps())],[1,1,3,3.0],1/16)
            Qout=100
            Pd,Qd=map(BigFloat,(0.0,0.0))
            σ₁=BigFloat(σ)
            for k=1:n
                D[k]>0.0 ? Pd=Pd+BigFloat(A.u[k])^2/(BigFloat(A.D[k])-σ₁) : Qd=Qd+BigFloat(A.u[k])^2/(BigFloat(A.D[k])-σ₁)
            end

            A.r > 0 ?   Pd+=BigFloat(1.0)/BigFloat(A.r)  : Qd+=BigFloat(1.0)/BigFloat(A.r)
            Kρ=Float64((Pd-Qd)/abs(Pd+Qd))
            @show Pd,Qd
            r=BigFloat(1.0)/(Pd+Qd)
            ρ=-Float64(r)
        end
    end
    # returns the following
    SymDPR1(D,u,ρ), Kρ, Qout

end # inv


function inv{T}(A::SymDPR1{T}, σ::Double)

    # COMPUTES: inverse of the σed SymDPR1 A = diagm(A.D)+A.r*A.u*A.u',
    # inv(A-σ*I) = D + ρ*u*u', σ!=A.D[i], which is again a SymDPR1
    # here σ is Double so it uses Double to compute everything
    # RETURNS: SymDPR1(D,u,ρ), Qout
    # Qout = 1 on exit meaning Double was used

    n=length(A.D)
    D=Array(Double,n)
    u=Array(Double,n)

    for k=1:n
        D[k]=Double(A.D[k])-σ
    end


    u=map(Double,A.u)

    oned=Double(1.0)
    zerod=Double(0.0)

    for k=1:n
        u[k]=u[k]/D[k]
        D[k]=oned/D[k]
    end


    # compute ρ
    # compute the sum in a plain loop

    P,Q=zerod,zerod
    Qout=1

    for k=1:n
        D[k].hi > 0.0 ? P=P+Double(A.u[k])*u[k] : Q=Q+Double(A.u[k])*u[k]
    end

    A.r > 0 ?   P+=Double(1.0)/Double(A.r)  : Q+=Double(1.0)/Double(A.r)

    r=oned/(P+Q)
    ρ=-(r.hi+r.lo)

    # returns the following
    # SymDPR1(T[x.hi+x.lo for x=D],T[x.hi+x.lo for x=u],ρ), Qout

    D1=Array{T}(n)
    u1=Array{T}(n)

    for k=1:n
        D1[k]=D[k].hi+D[k].lo
        u1[k]=u[k].hi+u[k].lo
    end

    SymDPR1(D1,u1,ρ), Qout

end # inv


function  eigen( A::SymDPR1{T},k::Integer,τ::Vector{Float64}) where T

    # COMPUTES: k-th eigenpair of an ordered irreducible SymDPR1
    # A = diagm(A.D)+A.r*A.u*A.u', A.r > 0
    # τ=[tolb,tolz,tolnu,tolrho,tollambda] = [1e3,10.0*n,1e3,1e3,1e3]
    # RETURNS: lambda, v, Sind, Kb, Kz, Knu, Kρ, Qout
    # lambda - k-th eigenvalue in descending order
    # v - lambda's normalized eigenvector
    # Kb, Kz, Knu, Kρ - condition numbers
    # Qout = 1 / 0 - Double was / was not used

    # Set the dimension
    n = length(A.D)

    # Set all conditions initially to zero
    Kb,Kz,Knu,Kρ=0.0,0.0,0.0,0.0
    Qout=0
    v=zeros(n)
    # Kz is former kappa_nu

    # Determine the σ sigma, the σ index i, and whether lambda
    # is on the left or the right side of the nearest pole
    # Exterior eigenvalues (k = 1 or k = n):

    if k == 1
        sigma,i,side = A.D[1],1,'R'
    else
        # Interior eigenvalues (k in (2,...n-1) ):
        Dtemp = A.D-A.D[k]
        middle = Dtemp[k-1]/2.0
        Fmiddle = 1.0+A.r*sum(A.u.^2./(Dtemp-middle))
        sigma,i,side = Fmiddle > 0.0 ? (A.D[k],k,'R') : (A.D[k-1],k-1,'L')
    end

    # Compute the inverse of the σed matrix, A_i^(-1), Kb and Kz

    Ainv,Kb,Kz,Qout = inv(A,i,τ[1:2])

    # Compute the eigenvalue of the inverse σed matrix

    nu = bisect( Ainv,side )

    #  nu=fastzero([invD1; 0; invD2], [w1;wz;w2], b, side); # implement later
    #  [nu-nu1]/nu, [nu-nueva]/nu, pause  # and compare

    if abs(nu)==Inf
        # this is nonstandard
        # Deflation in dpr1eig (nu=Inf)
        v[i]=1.0
        lambda=sigma
    else
        # standard case, full computation
        # nu1 is the F- or 1-norm of the inverse of the σed matrix
        # nu10=maximum([sum(abs,Ainv.z)+abs(Ainv.a), maximum(abs.(Ainv.D)+abs.(Ainv.z))])
        nu1=0.0
        for k=1:n-1
            nu1=max(nu1,abs(Ainv.D[k])+abs(Ainv.z[k]))
        end
        nu1=max(sum(abs,Ainv.z)+abs(Ainv.a), nu1)

        Knu=nu1/abs(nu)

        while Knu>τ[3]
            # Remedies according to Remark 3 - we σ between original
            # eigenvalues and compute DPR1 matrix
            # 1/nu1+sigma, 1/nu+sigma
            # println("Remedy 3 ")
            nu = side=='R' ? abs(nu) : -abs(nu)
            nu1=-sign(nu)*nu1
            sigma1=(nu1+nu)/(2.0*nu*nu1)+sigma

            if findfirst(A.D-sigma1,0.0)>0 # we came back with a pole
                # recompute sigmav more accurately according with Dekker
                sigmav=(Double(nu1)+Double(nu))/(Double(2.0)*Double(nu)*Double(nu1))+Double(sigma)
                # Compute the inverse of the σed arrowhead (DPR1)
                Ainv,Qout1=inv(A,sigmav) # Ainv is Float64
                nu1=bisect(Ainv,side)
                mu1 = 1.0/nu1
                D0=map(A.D,Double)-sigmav
                D1=D0.hi+D0.lo
                if findfirst(D1-mu1,0.0)>0
                    ind=find(D1-mu1==0.0);
                    v=zeros(n)
                    v[ind]=1.0;
                else
                    v=[ A.u./(D1-mu1)]
                end
                # σ the eigenvalue back in Double
                lam = Double(1.0)/Double(nu1)+sigmav
                # Return this
                lambda,v = lam.hi+lam.lo, v/norm(v)
                return lambda,v,i,Kb,Kz,Knu,Kρ,Qout
            else
                # Compute the inverse of the σed arrowhead (DPR1)
                Ainv, Kρ,Qout1=inv(A,sigma1,τ[4]) # Ainv is Float64
                # Compute the eigenvalue by bisect for DPR1
	        # Note: instead of bisection could use dlaed4 (need a wrapper) but
	        # it is not faster. There norm(u)==1
                nu= bisect(Ainv,side)
                # Check the code below for DPR1
                # if side=='R' && A.D[1]>0.0 && nu<0.0
                #     nu=bisect(Ainv,'L')
                # end
                mu=1.0/nu
                nu1=maximum(abs,Ainv.D)+abs(Ainv.r)*dot(Ainv.u,Ainv.u)
                Knu=nu1/abs(nu)
                sigma=sigma1
                # standard v
                # v=A.u./((A.D-sigma1)-mu1)
                # Return this - σ the eigenvalue back and normalize the vector
                # lambda, v = mu1+sigma1, v/norm(v)
            end
            Qout=Qout+2*Qout1
        end
        # Accuracy is fine, compute the eigenvector
        mu = 1.0/nu
        # v=[ A.u./((A.D-sigma)-mu)]
        for l=1:n
            v[l]= A.u[l]/((A.D[l]-sigma)-mu)
        end
        lambda, v = mu+sigma, v/norm(v)


        # Remedy according to Remark 1 - we recompute the the eigenvalue
        # near zero from the inverse of the original matrix (a DPR1 matrix).
        if (abs(A.D[i])+abs(1.0/nu))/abs(lambda)>τ[5]

            if (k==1 && A.D[1]<0.0 || side=='L' && sign(A.D[i])+sign(A.D[i+1])==0 || i>1 &&
                side=='R' && sign(A.D[i])+sign(A.D[i-1])==0)
                # println("Remedy 1 ")
                # Compute the inverse of the original arrowhead (DPR1)
                Ainv,Kρ,Qout1 = inv(A,0.0,τ[4]) # Ainv is Float64
                Qout=Qout+4*Qout1
                if abs(Ainv.r)==Inf
                    lambda=0.0
                else
                    # Here we do not need bisection. We compute the Rayleigh
                    # quotient by using already computed vectors which is
                    # componentwise accurate
                    nu1=sum(v.^2.*Ainv.D)+Ainv.r*sum(v.*Ainv.u)^2;
                    lambda=1.0/nu1
                end
            end

        end
    end

    # Return this
    lambda,v,i,Kb,Kz,Knu,Kρ,Qout

end # eig (k)


function eigen(A::SymDPR1{T}, τ::Vector{Float64}) where T

    # COMPUTES: all eigenvalues and eigenvectors of a real symmetric SymDPR1
    # A = Diagonal(A.D)+A.r*A.u*A.u'
    # τ = [tolb,tolz,tolnu,tolrho,tollambda] = [1e3,10.0*n,1e3,1e3,1e3] or similar
    # RETURNS: U, E, Sind, Kb, Kz, Knu, Kρ, Qout
    # U = eigenvectors, E = eigenvalues in decreasing order
    # Sind[k] - σ index i for the k-th eigenvalue
    # Kb, Kz, Knu, Kρ [k] - respective conditions for the k-th eigenvalue
    # Qout[k] = 1 / 0 - Double was / was not used when computing k-th eigenvalue

    n=length(A.D)
    n0=n

    # Checking if A.r > 0
    signr =  A.r > 0.0 ? 1.0 : -1.0

    # Ordering the matrix
    D=signr*A.D
    is=sortperm(D,rev=true)
    D=D[is]
    z=A.u[is]
    ρ=signr*A.r

    U=eye(n,n)
    E=zeros(n)
    Kb=zeros(n); Kz=zeros(n); Knu=zeros(n); Kρ=zeros(n)
    Qout=zeros(Int,n); Sind=zeros(Int,n)

    # Quick return for 1x1, this is trivial for SymArrow, not so trivial here :)
    if n==1
        U=1;
        if (D==0)&&((ρ==0)|(z==0))
            E=0
        else
            E=A.D[1]+A.r*A.u[1]^2
            # Higher accuracy if needed
            KD=(abs(A.D[1])+abs(A.r)*A.u[1]^2)/abs(E)
            if KD>τ[1]
                Ed=Double(A.D[1])+Double(A.r)*Double(A.u[1])^2
                E=Ed.hi+Ed.lo
            end
        end
        return U, [E], Sind,Kb,Kz,Knu,Kρ,Qout
    end

    #  test for deflation in z
    z0=find(z.==0)
    zx=find(z.!=0)

    if isempty(zx)  # nothing to do
        E=A.D
        isE=sortperm(E,rev=true)
        E=E[isE]
        U=U[:,isE]
        return U,E,Sind,Kb,Kz,Knu,Kρ,Qout
    end

    if !isempty(z0)
        E[z0]=D[z0]
        D=D[zx]
        z=z[zx]
        if !isempty(z)
            # return U,E,Sind,Kb,Kz,Knu,Kρ,Qout
            n=length(z)
        end
    end


    #  Test for deflation in D
    g=D[1:n-1]-D[2:n]
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
            R[l]=givens(z[g0[l]],z[g0[l]+1],zx[g0[l]],zx[g0[l]+1])
            z[g0[l]]=R[l][2]; z[g0[l]+1]=0.0
            # A_mul_Bc!(U,R) # multiply R'*U later
            E[zx[g0[l]+1]]=D[g0[l]+1]
        end
        # remains
        gx=[0;gx]+1
        nn=length(gx)

        zxx=zx[gx]
        for k=1:nn
            E[zxx[k]],U[zxx,zxx[k]],Sind[zxx[k]],Kb[zxx[k]],Kz[zxx[k]],Knu[zxx[k]],Kρ[zxx[k]],Qout[zxx[k]]=
            eig(SymDPR1(D[gx],z[gx],ρ),k,τ)
        end

        for l=1:lg0
            U=R[l][1]'*U
        end

    else

        # No deflation in D
        for k=1:n
            E[zx[k]],U[zx,zx[k]],Sind[zx[k]],Kb[zx[k]],Kz[zx[k]],Knu[zx[k]],Kρ[zx[k]],Qout[zx[k]]=
            eig(SymDPR1(D,z,ρ),k,τ)
        end
    end

    # back premutation of vectors
    isi=sortperm(is)

    # change the sign if A.r was negative
    # must sort E once more
    E=signr*E
    es=sortperm(E,rev=true)
    E=E[es]

    U=U[isi,es]

    # Return this
    E,U,Sind[es],Kb[es],Kz[es],Knu[es],Kρ[es],Qout[es]

end # eig (all)
