#------------------------
#--------  Functions

#--------  random generation
# using DoubleDouble

function doubledot(x::Vector{Float64}, y::Vector{Float64})
    dsum=Double(0.0)
    for k=1:length(x)
        dsum+=Double(x[k])*Double(y[k])
    end
    dsum
end

function GenHalfArrow(n::Integer, p::Int)
    # generates n x n half arrowhad matrix with/out arrow top
    # This is one without arrowhead point
    p==0 ?
    HalfArrow(rand(n-1).-0.5,rand(n-1).-0.5) :
    # This is one with the point
    HalfArrow(rand(n-1).-0.5,rand(n).-0.5)
end

#-------- Inverses

function inv(A::HalfArrow{T},i::Integer,τ::Vector{Float64}=[1e3, 10*length(A.D)]) where T
    # COMPUTES: inverse of a SymArrow matrix B, where B=A'*A, and A =[A.D, A.z]
    # is the HalfArrow. Here inv(B-B.D[i]*I) is again SymArrow.
    # Uses higher precision to compute top of the arrow element accurately, if
    # needed.
    # τ=[tolb,tolz] are tolerances, usually [1e3, 10*n]
    # [0.0,0.0] forces DoubleDouble, [1e50,1e50] would never use it
    # RETURNS:  SymArrow(D,z,b,i), Kb, Kz, Qout
    # Kb - condition Kb, Kz - condition Kz, Qout = 1 / 0 - double was / was not used

    n=length(A.D)
    D=Array{T}(undef,n)
    z=Array{T}(undef,n)
    z1=Array{T}(undef,n)
    z1.=A.z[1:n].*A.D
    wz=one(T)/(z1[i])
    σ=A.D[i]
    for k=1:i-1
        D[k]=one(T)/((A.D[k]-σ)*(A.D[k]+σ))
        z[k]=-z1[k]*D[k]*wz
    end
    for k=i+1:n
        D[k-1]=one(T)/((A.D[k]-σ)*(A.D[k]+σ))
        z[k-1]=-z1[k]*D[k-1]*wz
    end

    # D=[1./(A.D[1:i-1]-σ),1./(A.D[i+1:end]-σ)]
    # wz=1/A.z[i]
    # z=-[A.z[1:i-1], A.z[i+1:end]].*D1*wz1
    # Maybe this needs more accuracy??
    # if A.z[k]^2>A.D[i]^2 for some k, then we should use (A.z[k]+A.D[i])*(A.z[k]-A.d[i]) +
    # the rest, which is sum of squares
    # a=dot(A.z,A.z)-A.D[i]^2 # this was original formula
    ad=doubledot([A.z;A.D[i]],[A.z;-A.D[i]])
    a=ad.hi+ad.lo
    # compute the sum in a plain loop
    P=zero(T)
    Q=zero(T)
    Qout=0
    nn=n-1
    for k=1:i-1
        D[k]>0.0 ? P+=z1[k]^2*D[k] : Q+=z1[k]^2*D[k]
    end
    for k=i:nn
        D[k]>0.0 ? P+=z1[k+1]^2*D[k] : Q+=z1[k+1]^2*D[k]
    end
    a<0 ? P-=a : Q-=a

    Kb=(P-Q)/abs(P+Q)
    Kz=maximum(abs,z1)*abs(wz)
    if Kb<τ[1] ||  Kz<τ[2]
        b=(P+Q)*wz*wz
    else  # recompute in Double
        Qout=1
        σ₁=map(Double,A.D[i])
        Dd=[Double{Float64}[(Double(A.D[k])+σ₁)*(Double(A.D[k])-σ₁) for k=1:i-1];
            Double{Float64}[(Double(A.D[k])+σ₁)*(Double(A.D[k])-σ₁) for
                            k=i+1:length(A.D)]]
        wzd=Double(A.z[i])*σ₁
        # ad=doubledot(A.z,A.z)-σ₁^2 # we already have this if
        # we use doubledot above
        Pd,Qd=map(Double,(0.0,0.0))

        for k=1:i-1
            Dd[k].hi>0.0 ? Pd+=(Double(A.z[k])*Double(A.D[k]))^2/Dd[k] :
            Qd+=(Double(A.z[k])*Double(A.D[k]))^2/Dd[k]
        end
        for k=i:nn
            Dd[k].hi>0.0 ? Pd+=(Double(A.z[k+1])*Double(A.D[k+1]))^2/Dd[k] :
            Qd+=(Double(A.z[k+1])*Double(A.D[k+1]))^2/Dd[k]
            # @show P,Q
        end
        ad.hi<0 ?   Pd-=ad : Qd-=ad
        bd=(Pd+Qd)/(wzd*wzd)
        b=bd.hi+bd.lo
    end
    # SymArrow([D[1:A.i-2],zero(T),D[A.i-1:end]],[z[1:A.i-2],wz,z[A.i-1:end]],b,i),Kb,Kz,Qout
    # This reduces memory allocation (4 vectors) and is 30-40% faster
    D[n]=zero(T)
    z[n]=wz

    # return this
    SymArrow(D,z,b,i),Kb,Kz,Qout
end # inv

function inv!(B::SymArrow{T},A::HalfArrow{T},i::Integer,τ::Vector{Float64}=[1e3, 10*length(A.D)]) where T
    # In-place version!
    # COMPUTES: inverse of a SymArrow matrix B, where B=A'*A, and A =[A.D, A.z]
    # is the HalfArrow. Here inv(B-B.D[i]*I) is again SymArrow.
    # Uses higher precision to compute top of the arrow element accurately, if
    # needed.
    # τ=[tolb,tolz] are tolerances, usually [1e3, 10*n]
    # [0.0,0.0] forces DoubleDouble, [1e50,1e50] would never use it
    # RETURNS:  SymArrow(D,z,b,i), Kb, Kz, Qout
    # Kb - condition Kb, Kz - condition Kz, Qout = 1 / 0 - double was / was not used
    n=length(A.D)
    # D=Array{T}(undef,n)
    # z=Array{T}(undef,n)
    z1=Array{T}(undef,n)
    z1.=A.z[1:n].*A.D
    wz=one(T)/(z1[i])
    σ=A.D[i]
    for k=1:i-1
        B.D[k]=one(T)/((A.D[k]-σ)*(A.D[k]+σ))
        B.z[k]=-z1[k]*B.D[k]*wz
    end
    for k=i+1:n
        B.D[k-1]=one(T)/((A.D[k]-σ)*(A.D[k]+σ))
        B.z[k-1]=-z1[k]*B.D[k-1]*wz
    end

    ad=doubledot([A.z;A.D[i]],[A.z;-A.D[i]])
    a=ad.hi+ad.lo
    # compute the sum in a plain loop
    P=zero(T)
    Q=zero(T)
    Qout=0
    nn=n-1
    for k=1:i-1
        B.D[k]>0.0 ? P+=z1[k]^2*B.D[k] : Q+=z1[k]^2*B.D[k]
    end
    for k=i:nn
        B.D[k]>0.0 ? P+=z1[k+1]^2*B.D[k] : Q+=z1[k+1]^2*B.D[k]
    end
    a<0 ? P-=a : Q-=a

    Kb=(P-Q)/abs(P+Q)
    Kz=maximum(abs,z1)*abs(wz)
    if Kb<τ[1] ||  Kz<τ[2]
        b=(P+Q)*wz*wz
    else  # recompute in Double
        Qout=1
        σ₁=map(Double,A.D[i])
        Dd=[Double{Float64}[(Double(A.D[k])+σ₁)*(Double(A.D[k])-σ₁) for k=1:i-1];
            Double{Float64}[(Double(A.D[k])+σ₁)*(Double(A.D[k])-σ₁) for
                            k=i+1:length(A.D)]]
        wzd=Double(A.z[i])*σ₁
        # ad=doubledot(A.z,A.z)-σ₁^2 # we already have this if
        # we use doubledot above
        Pd,Qd=map(Double,(0.0,0.0))
        for k=1:i-1
            Dd[k].hi>0.0 ? Pd+=(Double(A.z[k])*Double(A.D[k]))^2/Dd[k] :
            Qd+=(Double(A.z[k])*Double(A.D[k]))^2/Dd[k]
        end
        for k=i:nn
            Dd[k].hi>0.0 ? Pd+=(Double(A.z[k+1])*Double(A.D[k+1]))^2/Dd[k] :
            Qd+=(Double(A.z[k+1])*Double(A.D[k+1]))^2/Dd[k]
            # @show P,Q
        end
        ad.hi<0 ?   Pd-=ad : Qd-=ad
        bd=(Pd+Qd)/(wzd*wzd)
        b=bd.hi+bd.lo
    end
    # SymArrow([D[1:A.i-2],zero(T),D[A.i-1:end]],[z[1:A.i-2],wz,z[A.i-1:end]],b,i),Kb,Kz,Qout
    # This reduces memory allocation (4 vectors) and is 30-40% faster
    B.D[n]=zero(T)
    B.z[n]=wz
    B.a=b
    B.i=i
    # return this
    # SymArrow(D,z,b,i),Kb,Kz,Qout
    return Kb, Kz, Qout
end # inv!

function inv(A::HalfArrow{T}, σ::Float64, τ::Float64=1.0e3) where T
    # COMPUTES: inverse of a SymArrow matrix B, where B=A'*A, and A =[A.D, A.z]
    # is the HalfArrow. Here inv(B-σ^2*I) is a SymDPR1.
    # uses DoubleDouble to compute ρ accurately, if needed.
    # τ is tolerance, usually 1e3,  0.0 forces Double, 1e50 would never use it
    # RETURNS: SymDPR1(D,u,ρ), Kρ, Qout
    # Kρ - condition Kρ, Qout = 1 / 0 - Double was / was not used

    n=length(A.D)
    D=Array{T}(undef,n+1)
    u=Array{T}(undef,n+1)
    z1=Array{T}(undef,n)
    z1.=A.z[1:n].*A.D
    for k=1:n
        D[k]=one(T)/((A.D[k]-σ)*(A.D[k]+σ))
        u[k]=z1[k]*D[k]
    end
    D[n+1]=zero(T)
    u[n+1]=-one(T)

    # compute ρ and Kρ
    #--- compute the sum in a plain loop
    P=zero(T)
    Q=zero(T)
    Qout=0
    # a=dot(A.z,A.z)-σ^2 # try with more accuracy
    ad=doubledot([A.z;σ],[A.z;-σ])
    a=ad.hi+ad.lo
    for k=1:n
        D[k]>0.0 ? P+=z1[k]^2*D[k] : Q+=z1[k]^2*D[k]
    end
    a<0 ? P-=a : Q-=a

    # Condition of ρ
    Kρ=(P-Q)/abs(P+Q)

    if Kρ<τ
        ρ=-1.0/(P+Q)
    else  # recompute in Double
        Qout=1
        Pd,Qd=map(Double,(0.0,0.0))
        σ₁=Double(σ)
        Dd=Double{Float64}[(Double(A.D[k])+σ₁)*(Double(A.D[k])-σ₁) for k=1:n]
        for k=1:n
            Dd[k].hi>0.0 ? Pd+=(Double(A.z[k])*Double(A.D[k]))^2/Dd[k] :
            Qd+=(Double(A.z[k])*Double(A.D[k]))^2/Dd[k]
        end
        ad.hi+ad.lo<0 ? Pd-=ad : Qd-=ad
        r=Double(1.0)/(Pd+Qd)
        ρ=-(r.hi+r.lo)
    end

    # returns the following
    SymDPR1(D,u,ρ), Kρ, Qout
end # inv

function inv(A::HalfArrow{T}, σ::Double) where T
    # COMPUTES: inverse of a SymArrow matrix B, where B=A'*A, and A =[A.D, A.z]
    # is the HalfArrow. Here inv(B-σ^2*I) is SymDPR1.
    # Here shift is Double so it uses Double to compute everything
    # RETURNS: SymDPR1(D1,u1,ρ), Qout
    # Qout = 1 on exit meaning Double was used

    n=length(A.D)
    D=Array{Double}(undef,n+1)
    u=Array{Double}(undef,n+1)
    oned=Double(1.0)
    zerod=Double(0.0)
    for k=1:n
        D[k]=(Double(A.D[k])-σ)*(Double(A.D[k])+σ)
        u[k]=Double(A.z[k])*Double(A.D[k])/D[k]
    end
    a=doubledot(A.z,A.z)-σ^2
    D[n+1]=zerod
    u[n+1]=Double(-1.0)

    # compute ρ and Kρ
    #--- compute the sum in a plain loop
    P,Q=zerod,zerod
    Qout=1
    for k=1:n
        D[k].hi>0.0 ? P+=(Double(A.z[k])*Double(A.D[k]))^2/D[k] :
        Q+=(Double(A.z[k])*Double(A.D[k]))^2/D[k]
        D[k]=oned/D[k]
    end
    a.hi+a.lo<0.0  ? P-=a : Q-=a
    r=oned/(P+Q)
    ρ=-(r.hi+r.lo)

    # returns the following
    D₁=Array{T}(undef,n+1)
    u₁=Array{T}(undef,n+1)
    for k=1:n+1
        D₁[k]=D[k].hi+D[k].lo
    end
    for k=1:n+1
        u₁[k]=u[k].hi+u[k].lo
    end
    SymDPR1(D₁,u₁,ρ), Qout
end # inv


function  svd(A::HalfArrow{T},Ainv::SymArrow{T},k::Integer,
    τ::Vector{Float64}=[1e3,10.0*length(A.D),1e3,1e3,1e3]) where T
    # COMPUTES: k-th singular value triple of an ordered irreducible HalfArrow
    # A = [Diagonal(A.D) A.z] with A.D > 0
    # τ=[tolb,tolz,tolnu,tolrho,tollambda] = [1e3,10.0*n,1e3,1e3,1e3]
    # RETURNS: λ, u, v, Info
    # where
    # λ - k-th singular value in descending order
    # u - λ's normalized left singular vector
    # v - λ's normalized right singular vector
    # Info = Sind, Kb, Kz, Kν, Kρ, Qout
    # Sind = shift index i for the k-th singular value
    # Kb, Kz, Kν, Kρ - condition numbers
    # Qout = 1 / 0 - Double was / was not used

    # Set the dimension
    n = length(A.D) + 1    # Set all conditions initially to zero
    Kb,Kz,Kν,Kρ=0.0,0.0,0.0,0.0
    Qout=0
    u=zeros(length(A.z))
    v=zeros(n)
    # z1=Array{Float64}(undef,n)
    z1=A.z[1:n-1].*A.D

    # Kz is former kappa_nu
    # Determine the shift σ, the shift index i, and whether λ
    # is on the left or the right side of the nearest pole
    # Exterior eigenvalues (k = 1 or k = n):

    if k == 1
        σ,i,side = A.D[1],1,'R'
    elseif k==n
        σ,i,side = A.D[n-1],n-1,'L'
    else
        # Interior eigenvalues (k in (2,...n-1) ):
        Dtemp = (A.D.-A.D[k]).*(A.D.+A.D[k])
        atemp = dot(A.z,A.z)-A.D[k]^2
        middle = Dtemp[k-1]/2.0
        Fmiddle = (atemp-middle)-sum(z1.^2 ./(Dtemp.-middle))
        σ,i,side = Fmiddle < 0.0 ? (A.D[k],k,'R') : (A.D[k-1],k-1,'L')
    end

    # Compute the inverse of the shifted matrix, (A'*A)_i^(-1), Kb and Kz
    # Ainv,Kb,Kz,Qout = inv(A,i)
    Kb,Kz,Qout = inv!(Ainv,A,i)

    # Compute the eigenvalue of the inverse shifted matrix
    ν = bisect( Ainv,side )

    #  ν=fastzero([invD1; 0; invD2], [w1;wz;w2], b, side); # implement later
    #  [ν-ν₁]/ν, [ν-nueva]/ν, pause  # and compare
    if abs(ν)==Inf
        # this is nonstandard
        # Deflation in aheig (ν=Inf)
        u[i]=1.0
        v[i]=1.0
        λ=σ
    else
        # standard case, full computation
        # ν₁ is the F- or 1-norm of the inverse of the shifted matrix
        # ν₁0=maximum([sum(abs,Ainv.z)+abs(Ainv.a), maximum(abs,Ainv.D)+abs.(Ainv.z))])
        ν₁=0.0
        for k=1:n-1
            ν₁=max(ν₁,abs(Ainv.D[k])+abs(Ainv.z[k]))
        end
        ν₁=max(sum(abs,Ainv.z)+abs(Ainv.a), ν₁)
        Kν=ν₁/abs(ν)
        if Kν<τ[3]
            # Accuracy is fine, compute the eigenvector
            μ = 1.0/ν
            #         v=[ z1./((A.D-σ).*(A.D+σ)-μ);-1.0]
            for k=1:n-1
	        v[k] = z1[k]/((A.D[k]-σ)*(A.D[k]+σ)-μ)
            end
            v[n]=-1.0
            normalize!(v)
            λ=sqrt(μ+σ^2) # this may have errors
            u[1:n-1]=λ*v[1:n-1]./A.D
            if length(A.z)==n
                u[n]=A.z[n]*v[n]/λ
            end
        else
            # Remedies according to Remark 3 - we shift between original
            # eigenvalues and compute DPR1 matrix
            # 1/ν₁+σ, 1/ν+σ
            ν = side=='R' ? abs(ν) : -abs(ν)
            ν₁=-sign(ν)*ν₁
            σ₁=(ν₁+ν)/(2.0*ν*ν₁)+σ^2

            if count(isequal(sqrt(σ₁)),abs.(A.D))>0 # we came back with a pole
                # recompute σᵥ more accurately according with dekker
                σᵥ=(Double(ν₁)+Double(ν))/(Double(2.0)*Double(ν)*Double(ν₁))+Double(σ)^2
                # Compute the inverse of the shifted arrowhead (DPR1)
                Ainv,Qout₁=inv(A,sqrt(σᵥ)) # Ainv is Float64, here it was sqrt(σᵥ)
                ν₁=bisect(Ainv,side)
                μ₁ = 1.0/ν₁
                D0=map(A.D,Double)
                D0=A.D.*A.D-σᵥ
                D1=D0.hi+D0.lo
                # Shift the eigenvalue back in Double
                lam = Double(1.0)/Double(ν₁)+σᵥ
                λ=sqrt(lam.hi+lam.lo)

                if count(isequal(μ₁),D1)>0
                    ind=findall(isequal(μ₁),D1)
                    u[ind]=1.0
                    v[ind]=1.0
                else
                    for k=1:n-1
	                v[k] = z1[k]/(D1[k]-μ₁)
                    end
                    v[n]=-1.0
                    normalize!(v)
                    u[1:n-1]=λ*v[1:n-1]./A.D
                    if length(A.z)==n
                        u[n]=A.z[n]*v[n]/λ
                    end
                    # v=[ A.z./(D1-μ₁);-1.0]
                end
            else
                # Compute the inverse of the shifted arrowhead (DPR1)
                σₛ=sqrt(σ₁)
                Ainv, Kρ,Qout₁=inv(A,σₛ,τ[4]) # Ainv is Float64
                # Compute the eigenvalue by bisect for DPR1
	        # Note: instead of bisection could use dlaed4 (need a wrapper) but
	        # it is not faster. There norm(u)==1
                ν₁= bisect(Ainv,side)
                μ₁=1.0/ν₁
                for k=1:n-1
	            v[k] = z1[k]/((A.D[k]-σₛ)*(A.D[k]+σₛ)-μ₁)
                end
                v[n]=-1.0
                normalize!(v)
                λ=sqrt(μ₁+σ₁)
                u[1:n-1]=λ*v[1:n-1]./A.D
                if length(A.z)==n
                    u[n]=A.z[n]*v[n]/λ
                end
            end
            Qout==1 && (Qout=Qout+2)
        end

        # Remedy according to Remark 1 - we recompute the the eigenvalue
        # near zero from the inverse of the original matrix (a DPR1 matrix).
        # Here this can only happen for k==n
        if k==n && (abs(A.D[i])^2+abs(1.0/ν))/abs(λ^2)>τ[5]
            # if k==1 && A.D[1]<0.0 || k==n && A.D[n-1]>0.0 || side=='L' && sign(A.D[i])+sign(A.D[i+1])==0 || side=='R' && sign(A.D[i])+sign(A.D[i-1])==0
            # Compute the inverse of the original arrowhead (DPR1)
            Ainv,Kρ,Qout₁ = inv(A,0.0,τ[4]) # Ainv is Float64
            Qout==1 && (Qout=Qout+4)
            if abs(Ainv.r)==Inf
                λ=0.0
            else
                # Here we do not need bisection. We compute the Rayleigh
                # quotient by using already computed vectors which are
                # componentwise accurate
                ν₁=sum(v.^2 .*Ainv.D)+Ainv.r*sum(v.*Ainv.u)^2;
                λ=sqrt(1.0/ν₁)
            end
        end
    end
    # Return this
    κ=Info(i,Kb,Kz,Kν,Kρ,Qout)
    return λ,u,v,κ
end # svd (k)

function svd(A::HalfArrow{T}, τ::Vector{Float64}=[1e3,10.0*length(A.D),1e3,1e3,1e3]) where T
    # COMPUTES: all singular values and singular vectors of a real HalfArrow
    # A = [diagm(A.D) A.z]
    # τ=[tolb,tolz,tolnu,tolrho,tollambda] = [1e3,10.0*n,1e3,1e3,1e3]
    # RETURNS: SVD(U, Σ, V), κ
    # where
    # U = left singular vectors
    # Σ = singular values in decreasing order
    # V = right singular vectors
    # κ[k]=Info(Sind[k], Kb[k], Kz[k], Kν[k], Kρ[k], Qout[k])
    # Sind[k] - shift index i for the k-th eigenvalue
    # Kb, Kz, Kν, Kρ [k] - respective conditions for the k-th singular value
    # Qout[k] = 1 / 0 - Double was / was not used when computing k-th singular value

    nd=length(A.D)
    n=length(A.z)
    n0=n
    # Ordering the matrix
    signD=sign.(A.D)
    D=abs.(A.D)
    is=sortperm(D,rev=true)
    D=D[is]
    signD=signD[is]
    z=A.z[1:n]
    z=z[is]
    if n==nd
        U=Array{T}(I,n,n)
        V=Array{T}(I,n+1,n)
    else
        U=Array{T}(I,n,n)
        V=Array{T}(I,n,n)
    end
    Σ=zeros(n)
    κ=[Info(0, 0.0, 0.0, 0.0, 0.0, 0) for i=1:n]

    #  test for deflation in z
    z0=findall(iszero,z)
    zx=findall(!iszero,z)
    if isempty(zx)  # nothing to do
        if n==nd
            Σ=[A.D]
            isΣ=sortperm(Σ,rev=true)
            Σ=Σ[isΣ]
            V[:,1:n]=view(V,:,isΣ) # V[:,isΣ]
            U=view(U,invperm(isΣ),:) # U[invperm(isΣ),:]
        else
            Σ=[A.D;z[n]]
            isΣ=sortperm(Σ,rev=true)
            Σ=Σ[isΣ]
            V=view(V,:,isΣ) # V[:,isΣ]
            U=view(U,invperm(isΣ),:) # U[invperm(isΣ),:]
        end
        return SVD(U,Σ,V), κ
    end

    if !isempty(z0)
        Σ[z0]=D[z0]
        for k=1:length(z0)
            V[z0[k],z0[k]]=signD[z0[k]]
        end
        D=D[zx]
        # signD=signD[zx]
        z=z[zx]
        if isempty(z)
            if n>nd
                Σ[n]=abs(A.z[n])
                V[n,n]=sign(A.z[n])
            end
            # return U,Σ,V,Sind,Kb,Kz,Kν,Kρ,Qout
        end
    else
        if n0==nd
            zxv=[zx;n0+1]
        else
            zx=zxv=[zx;n0]
            z=[z;A.z[n0]]
        end
        n=length(z)

        #  Test for deflation in D
        g=D[1:n-2]-D[2:n-1]
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
            R=Array{Tuple{LinearAlgebra.Givens{Float64},Float64}}(undef,lg0)
            for l=lg0:-1:1
                R[l]=givens(z[g0[l]],z[g0[l]+1],zx[g0[l]],zx[g0[l]+1])
                z[g0[l]]=R[l][2]; z[g0[l]+1]=0.0
                # A_mul_Bc!(U,R) # multiply R'*U later
                Σ[zx[g0[l]+1]]=D[g0[l]+1]
            end
            # remains
            gx=[0;gx].+1
            nn=length(gx)
            zxx=zx[[gx;n]]
            if n0==nd
                zxxv=[zxx,n0+1]
            else
                zxx=zxxv=[zxx,n0]
                # z=[z;A.z[n0]]
            end
            Axx=HalfArrow(D[gx],z[gx])
            Bxx=[SymArrow(Vector{T}(undef,nn),Vector{T}(undef,nn),one(T),nn+1) for i=1:Threads.nthreads()]
            # Threads.@threads
            for k=1:nn+1
                tid=Threads.threadid()
                Σ[zxx[k]],U[zxx,zxx[k]],V[zxxv,zxx[k]],κ[zxx[k]]=svd(Axx,Bxx[tid],k)
            end
            for l=1:lg0
                # U=R[l][1]'*U
                lmul!(R[l][1],U)
            end
        else
            # No deflation in D
            Ax=HalfArrow(D,z)
            Bx=[SymArrow(Vector{T}(undef,n),Vector{T}(undef,n),one(T),n+1) for i=1:Threads.nthreads()]
            # Threads.@threads
            for k=1:n
                tid=Threads.threadid()
                Σ[zx[k]],U[zx,zx[k]],V[zxv,zx[k]],κ[zx[k]]=svd(Ax,Bx[tid],k)
                V[zxv[1:end-1],zx[k]]=V[zxv[1:end-1],zx[k]].*signD[zxv[1:end-1]]
            end
        end
    end
    # back premutation of vectors
    isi=sortperm(is)
    # must sort Σ once more
    es=sortperm(Σ,rev=true)
    Σ.=Σ[es]
    if n0==nd
        U=view(U,isi,es)
        V=view(V,[isi;n0+1],es)
    else
        isiv=[isi;n0]
        U=view(U,isiv,es)
        V=view(V,isiv,es)
    end
    κ.=κ[es]
    # Return this
    return SVD(U,Σ,V),κ
end # svd (all)
