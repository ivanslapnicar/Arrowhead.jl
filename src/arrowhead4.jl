#-------- Inverses of DPR1

function inv(A::SymDPR1{T},i::Integer,τ::Vector{Float64}=[1e3,10.0*length(A.D)]) where T
    # COMPUTES: inverse of a shifted SymDPR1 matrix A=diagm(A.D)+A.r*A.u*A.u',
    # inv(A-A.D[i]*I) which is a SymArrow.
    # Uses higher precision to compute top of the arrow element accurately, if
    # needed.
    # τ=[tolb,tolz] are tolerances, usually [1e3, 10*n]
    # [0.0,0.0] forces DoubleDouble, [1e50,1e50] would never use it
    # RETURNS:  SymArrow(D,z,b,i), Kb, Kz, Qout
    # Kb - condition Kb, Kz - condition Kz, Qout = 1 / 0 - double was / was not used
    n=length(A.D)
    D=Array{T}(undef,n-1)
    z=Array{T}(undef,n-1)
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

function inv!(B::SymArrow{T},A::SymDPR1{T},i::Integer,τ::Vector{Float64}=[1e3,10.0*length(A.D)]) where T
    # In-place version!
    # COMPUTES: inverse of a shifted SymDPR1 matrix A=diagm(A.D)+A.r*A.u*A.u',
    # inv(A-A.D[i]*I) which is a SymArrow.
    # Uses higher precision to compute top of the arrow element accurately, if
    # needed.
    # τ=[tolb,tolz] are tolerances, usually [1e3, 10*n]
    # [0.0,0.0] forces DoubleDouble, [1e50,1e50] would never use it
    # RETURNS:  SymArrow(D,z,b,i), Kb, Kz, Qout
    # Kb - condition Kb, Kz - condition Kz, Qout = 1 / 0 - double was / was not used
    n=length(A.D)
    # D=Array{T}(undef,n-1)
    # z=Array{T}(undef,n-1)
    wz=one(T)/A.u[i]
    σ=A.D[i]
    for k=1:i-1
        B.D[k]=one(T)/(A.D[k]-σ)
        B.z[k]=-A.u[k]*B.D[k]*wz
    end
    for k=i+1:n
        B.D[k-1]=one(T)/(A.D[k]-σ)
        B.z[k-1]=-A.u[k]*B.D[k-1]*wz
    end

    # compute the sum in a plain loop
    P=zero(T)
    Q=zero(T)
    Qout=0
    for k=1:i-1
        B.D[k]>0.0 ? P+=A.u[k]^2*B.D[k] : Q+=A.u[k]^2*B.D[k]
    end
    for k=i:n-1
        B.D[k]>0.0 ? P+=A.u[k+1]^2*B.D[k] : Q+=A.u[k+1]^2*B.D[k]
    end
    A.r>0 ? P+=one(T)/A.r : Q+=one(T)/A.r
    Kb=(P-Q)/abs(P+Q)

    # compute Kz
    Kz=(sum(abs,A.u)-abs(A.u[i]))*abs(wz)

    # Kz=max(abs,A.z)*abs(wz)
    if Kb<τ[1] ||  Kz<τ[2]
        b=(P+Q)*wz*wz
        B.a=b
        B.i=i
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
        B.a=b
        B.i=i
    end
    # SymArrow(D,z,b,i),Kb,Kz,Qout
    return Kb,Kz,Qout
end # inv!

function inv(A::SymDPR1{T}, σ::Float64, τ::Float64=1.0e3) where T
    # COMPUTES: inverse of the σed SymDPR1 A = diagm(A.D)+A.r*A.u*A.u',
    # inv(A-σ*I) = D + ρ*u*u', σ!=A.D[i], which is again a SymDPR1
    # uses DoubleDouble to compute A.r accurately, if needed.
    # τ is tolerance, usually 1e3,  0.0 forces Double, 1e50 would never use it
    # RETURNS: SymDPR1(D,u,ρ), Kρ, Qout
    # Kρ - condition Kρ, Qout = 1 / 0 - Double was / was not used
    n=length(A.D)
    D=Array{T}(undef,n)
    u=Array{T}(undef,n)
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
        D[k]>0.0 ? P+=A.u[k]^2*D[k] : Q+=A.u[k]^2*D[k]
    end
    A.r>0 ? P+=one(T)/A.r : Q+=one(T)/A.r
    # Condition of ρ
    Kρ=(P-Q)/abs(P+Q)
    if Kρ < τ
        ρ=-one(T)/(P+Q)
    else  # recompute in Double
        Qout=1
        Pd,Qd=map(Double,(0.0,0.0))
        σ₁=Double(σ)
        for k=1:n
            D[k]>0.0 ? Pd+=Double(A.u[k])^2/(Double(A.D[k])-σ₁) : Qd+=Double(A.u[k])^2/(Double(A.D[k])-σ₁)
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
                D[k]>0.0 ? Pd+=BigFloat(A.u[k])^2/(BigFloat(A.D[k])-σ₁) : Qd+=BigFloat(A.u[k])^2/(BigFloat(A.D[k])-σ₁)
            end
            A.r > 0 ?   Pd+=BigFloat(1.0)/BigFloat(A.r)  : Qd+=BigFloat(1.0)/BigFloat(A.r)
            Kρ=Float64((Pd-Qd)/abs(Pd+Qd))
            r=BigFloat(1.0)/(Pd+Qd)
            ρ=-Float64(r)
        end
    end
    # returns the following
    SymDPR1(D,u,ρ), Kρ, Qout
end # inv

function inv(A::SymDPR1{T}, σ::Double) where T
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

function  eigen( A::SymDPR1{T}, Ainv::SymArrow{T}, k::Integer,
    τ::Vector{Float64}=[1e3,10.0*length(A.D),1e3,1e3,1e3]) where T
    # COMPUTES: k-th eigenpair of an ordered irreducible SymDPR1
    # A = Diagonal(A.D)+A.r*A.u*A.u', A.r > 0
    # τ=[tolb,tolz,tolnu,tolrho,tollambda] = [1e3,10.0*n,1e3,1e3,1e3]
    # RETURNS: λ, v, Info
    # where
    # λ - k-th eigenvalue in descending order
    # v - λ's normalized eigenvector
    # Info = Sind, Kb, Kz, Kν, Kρ, Qout
    # Sind = shift index i for the k-th eigenvalue
    # Kb, Kz, Kν, Kρ - condition numbers
    # Qout = 1 / 0 - Double was / was not used
    # Set the dimension
    n = length(A.D)
    # Set all conditions initially to zero
    Kb,Kz,Kν,Kρ=0.0,0.0,0.0,0.0
    Qout=0
    v=zeros(T,n)
    # Kz is former kappa_nu
    # Determine the shift σ, the shift index i, and whether λ
    # is on the left or the right side of the nearest pole
    # Exterior eigenvalues (k = 1 or k = n):
    if k == 1
        σ,i,side = A.D[1],1,'R'
    else
        # Interior eigenvalues (k in (2,...n-1) ):
        Dtemp = A.D.-A.D[k]
        middle = Dtemp[k-1]/2.0
        Fmiddle = 1.0+A.r*sum(A.u.^2 ./(Dtemp.-middle))
        σ,i,side = Fmiddle > 0.0 ? (A.D[k],k,'R') : (A.D[k-1],k-1,'L')
    end

    # Compute the inverse of the shifted matrix, A_i^(-1), Kb and Kz
    # Ainv,Kb,Kz,Qout = inv(A,i)
    Kb,Kz,Qout = inv!(Ainv,A,i)
    # Compute the eigenvalue of the inverse σed matrix
    ν = bisect( Ainv,side )

    #  ν=fastzero([invD1; 0; invD2], [w1;wz;w2], b, side); # implement later
    #  [ν-ν₁]/ν, [ν-νeva]/ν, pause  # and compare
    if abs(ν)==Inf
        # this is nonstandard
        # Deflation in dpr1eig (ν=Inf)
        v[i]=1.0
        λ=sigma
    else
        # standard case, full computation
        # ν₁ is the F- or 1-norm of the inverse of the σed matrix
        # ν₁0=maximum([sum(abs,Ainv.z)+abs(Ainv.a), maximum(abs.(Ainv.D)+abs.(Ainv.z))])
        ν₁=0.0
        for k=1:n-1
            ν₁=max(ν₁,abs(Ainv.D[k])+abs(Ainv.z[k]))
        end
        ν₁=max(sum(abs,Ainv.z)+abs(Ainv.a), ν₁)
        Kν=ν₁/abs(ν)
        while Kν>τ[3]
            # Remedies according to Remark 3 - we σ between original
            # eigenvalues and compute DPR1 matrix
            # 1/ν₁+sigma, 1/ν+sigma
            # println("Remedy 3 ")
            ν = side=='R' ? abs(ν) : -abs(ν)
            ν₁=-sign(ν)*ν₁
            σ₁=(ν₁+ν)/(2.0*ν*ν₁)+σ
            if count(isequal(σ₁),A.D)>0 # we came back with a pole
                # recompute sigmav more accurately according with Dekker
                σᵥ=(Double(ν₁)+Double(ν))/(Double(2.0)*Double(ν)*Double(ν₁))+Double(σ)
                # Compute the inverse of the shifted arrowhead (DPR1)
                Ainv,Qout₁=inv(A,σᵥ) # Ainv is Float64
                ν₁=bisect(Ainv,side)
                μ₁ = 1.0/ν₁
                D0=map(A.D,Double)-σᵥ
                D1=D0.hi+D0.lo
                if count(isequal(μ₁),D1)>0
                    ind=findall(isequal(μ₁),D1)
                    v=zeros(T,n)
                    v[ind]=1.0
                else
                    v.=A.u./(D1.-μ₁)
                end
                # σ the eigenvalue back in Double
                lam = Double(1.0)/Double(ν₁)+σᵥ

                # Return this
                λ = lam.hi+lam.lo
                normalize!(v)
                κ=Info(i,Kb,Kz,Kν,Kρ,Qout)
                return λ,v,κ
            else
                # Compute the inverse of the Shifted DPR1 (DPR1)
                Ainv₁, Kρ,Qout₁=inv(A,σ₁,τ[4]) # Ainv is Float64
                # Compute the eigenvalue by bisect for DPR1
	            # Note: instead of bisection could use dlaed4 (need a wrapper) but
	            # it is not faster. There norm(u)==1
                ν= bisect(Ainv₁,side)
                # Check the code below for DPR1
                # if side=='R' && A.D[1]>0.0 && ν<0.0
                #     ν=bisect(Ainv,'L')
                # end
                μ=1.0/ν
                ν₁=maximum(abs,Ainv₁.D)+abs(Ainv₁.r)*dot(Ainv₁.u,Ainv₁.u)
                Kν=ν₁/abs(ν)
                σ=σ₁
                # standard v
                # v=A.u./((A.D-σ₁)-μ₁)
                # Return this - σ the eigenvalue back and normalize the vector
                # λ, v = μ₁+σ₁, v/norm(v)
            end
            Qout=Qout+2*Qout₁
        end
        # Accuracy is fine, compute the eigenvector
        μ = 1.0/ν
        # v=[ A.u./((A.D-σ)-μ)]
        for l=1:n
            v[l]= A.u[l]/((A.D[l]-σ)-μ)
        end
        λ = μ + σ
        normalize!(v)

        # Remedy according to Remark 1 - we recompute the the eigenvalue
        # near zero from the inverse of the original matrix (a DPR1 matrix).
        if (abs(A.D[i])+abs(1.0/ν))/abs(λ)>τ[5]
            if (k==1 && A.D[1]<0.0 || side=='L' && sign(A.D[i])+sign(A.D[i+1])==0 || i>1 &&
                side=='R' && sign(A.D[i])+sign(A.D[i-1])==0)
                # println("Remedy 1 ")
                # Compute the inverse of the original arrowhead (DPR1)
                Ainv₁,Kρ,Qout₁ = inv(A,0.0,τ[4]) # Ainv is Float64
                Qout=Qout+4*Qout₁
                if abs(Ainv₁.r)==Inf
                    λ=0.0
                else
                    # Here we do not need bisection. We compute the Rayleigh
                    # quotient by using already computed vectors which is
                    # componentwise accurate
                    ν₁=sum(v.^2 .*Ainv₁.D)+Ainv₁.r*sum(v.*Ainv₁.u)^2;
                    λ=1.0/ν₁
                end
            end
        end
    end
    # Return this
    κ=Info(i,Kb,Kz,Kν,Kρ,Qout)
    λ,v,κ
end # eigen (k)

function eigen(A::SymDPR1{T}, τ::Vector{Float64}=[1e3,10.0*length(A.D),1e3,1e3,1e3]) where T
    # COMPUTES: all eigenvalues and eigenvectors of a real symmetric SymDPR1
    # A = Diagonal(A.D)+A.r*A.u*A.u'
    # τ = [tolb,tolz,tolnu,tolrho,tollambda] = [1e3,10.0*n,1e3,1e3,1e3] or similar
    # RETURNS: Eigen(Λ,U), κ
    # where
    # Λ = eigenvalues in decreasing order, U = eigenvectors,
    # κ[k]=Info(Sind[k], Kb[k], Kz[k], Kν[k], Kρ[k], Qout[k]
    # Sind[k] - shift index i for the k-th eigenvalue
    # Kb, Kz, Kν, Kρ [k] - respective conditions for the k-th eigenvalue
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
    # Eigenvecgtor matrix and eigenvalue vector
    U=Array{T}(I,n,n)
    Λ=zeros(T,n)
    κ=[Info(0, 0.0, 0.0, 0.0, 0.0, 0) for i=1:n]

    # Quick return for 1x1, this is trivial for SymArrow, not so trivial here :)
    if n==1
        if (D==0)&&((ρ==0)|(z==0))
            Λ=0
        else
            Λ=A.D[1]+A.r*A.u[1]^2
            # Higher accuracy if needed
            KD=(abs(A.D[1])+abs(A.r)*A.u[1]^2)/abs(Λ)
            if KD>τ[1]
                Λd=Double(A.D[1])+Double(A.r)*Double(A.u[1])^2
                Λ=Λd.hi+Λd.lo
            end
            # Qout=1
            κ[1]=Info(0, 0.0, 0.0, 0.0, 0.0, 1)
        end
        return Eigen([Λ], U), κ
    end

    #  test for deflation in z
    z0=findall(iszero,z)
    zx=findall(!iszero,z)
    if isempty(zx)  # nothing to do
        Λ=A.D
        isΛ=sortperm(Λ,rev=true)
        Λ=Λ[isΛ]
        U=view(U,:,isΛ) # [:,isΛ]
        return U,Λ,Sind,Kb,Kz,Kν,Kρ,Qout
    end
    if !isempty(z0)
        Λ[z0]=D[z0]
        D=D[zx]
        z=z[zx]
        if !isempty(z)
            # return U,Λ,Sind,Kb,Kz,Kν,Kρ,Qout
            n=length(z)
        end
    end

    #  Test for deflation in D
    g=D[1:n-1]-D[2:n]
    # Can play with inexact deflation
    # g0=find(abs(g)<eps)
    # gx=find(abs(g)>=eps)
    # Only exact deflation !!
    g0=findall(iszero,g)
    gx=findall(!iszero,g)
    if !isempty(g0)
        # Deflation
        Dgx=D[gx]; zgx=z[gx]
        lg0=length(g0)
        R=Array(Tuple{LinearAlgebra.Givens{Float64},Float64},lg0)
        for l=lg0:-1:1
            R[l]=givens(z[g0[l]],z[g0[l]+1],zx[g0[l]],zx[g0[l]+1])
            z[g0[l]]=R[l][2]; z[g0[l]+1]=0.0
            # A_mul_Bc!(U,R) # multiply R'*U later
            Λ[zx[g0[l]+1]]=D[g0[l]+1]
        end
        # remains
        gx=[0;gx].+1
        nn=length(gx)
        zxx=zx[gx]
        Axx=SymDPR1(D[gx],z[gx],ρ)
        # Bxx=SymArrow(Vector{T}(undef,nn-1),Vector{T}(undef,nn-1),one(T),nn)
        Bxx=[SymArrow(Vector{T}(undef,nn-1),Vector{T}(undef,nn-1),one(T),nn) for i=1:Threads.nthreads()]
        Threads.@threads for k=1:nn
            tid=Threads.threadid()
            Λ[zxx[k]],U[zxx,zxx[k]],κ[zxx[k]]=eigen(Axx,Bxx[tid],k)
        end
        for l=1:lg0
            # U=R[l][1]'*U
            lmul!(R[l][1],U)
        end
    else
        # No deflation in D
        Ax=SymDPR1(D,z,ρ)
        Bx=[SymArrow(Vector{T}(undef,n-1),Vector{T}(undef,n-1),one(T),n) for i=1:Threads.nthreads()]
        # Bx=SymArrow(Vector{T}(undef,n-1),Vector{T}(undef,n-1),one(T),n)
        Threads.@threads for k=1:n
            tid=Threads.threadid()
            Λ[zx[k]],U[zx,zx[k]],κ[zx[k]]=eigen(Ax,Bx[tid],k)
        end
    end

    # back premutation of vectors
    isi=sortperm(is)
    # change the sign if A.r was negative
    # must sort Λ once more
    Λ=signr*Λ
    es=sortperm(Λ,rev=true)
    Λ.=Λ[es]
    # U=U[isi,es]
    # Return this
    U=view(U,isi,es)
    κ.=κ[es]
    # Return this
    return Eigen(Λ,U),κ
end # eigen (all)
