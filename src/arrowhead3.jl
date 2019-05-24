#------------------------
#--------  Functions

#--------  random generation

using Random
function GenSymArrow(n::Integer,i::Integer)
    # generates symmetric n x n arrowhad matrix with arrow top at (i,i)
    SymArrow(rand(n-1),rand(n-1),rand(),i)
end

function GenSymDPR1(n::Integer)
    # generates symmetric n x n DPR1 matrix
    SymDPR1(rand(n),rand(n),rand())
end

#-------- Inverses
function inv(A::SymArrow{T},i::Integer,τ::Vector{Float64}=[1e3;10.0*length(A.D)]) where T
    # COMPUTES: inverse of a SymArrow matrix A, inv(A-A.D[i]*I) which is again SymArrow
    # uses higher precision to compute top of the arrow element accurately, if
    # needed.
    # τ=[tolb,tolz] are tolerances, usually [1e3, 10*n]
    # [0.0,0.0] forces DoubleDouble, [1e50,1e50] would never use it
    # RETURNS:  SymArrow(D,z,b,i), Kb, Kz, Qout
    # Kb - condition Kb, Kz - condition Kz, Qout = 1 / 0 - double was / was not used

    n=length(A.D)
    D=Array{T}(undef,n)
    z=Array{T}(undef,n)
    wz=one(T)/A.z[i]
    σ=A.D[i]
    for k=1:i-1
        D[k]=one(T)/(A.D[k]-σ)
        z[k]=-A.z[k]*D[k]*wz
    end
    for k=i+1:n
        D[k-1]=one(T)/(A.D[k]-σ)
        z[k-1]=-A.z[k]*D[k-1]*wz
    end
    # D=[1./(A.D[1:i-1]-shift),1./(A.D[i+1:end]-σ)]
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
    a<0 ? P-=a : Q-=a
    Kb=(P-Q)/abs(P+Q)
    Kz=maximum(abs,A.z)*abs(wz)
    if Kb<τ[1] ||  Kz<τ[2]
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
        σ₁=map(Type,A.D[i])
        Dd=[[Type(A.D[k])-σ₁ for k=1:i-1];
            [Type(A.D[k])-σ₁ for k=i+1:length(A.D)]]
        wzd=Type(A.z[i])
        ad=Type(A.a)-σ₁

        Pd,Qd=map(Type,(0.0,0.0))
        for k=1:i-1
            convert(Float64,Dd[k])>0.0 ? Pd+=Type(A.z[k])^2/Dd[k] :
            Qd+=Type(A.z[k])^2/Dd[k]
        end
        for k=i:nn
            convert(Float64,Dd[k])>0.0 ? Pd+=Type(A.z[k+1])^2/Dd[k] :
            Qd+=Type(A.z[k+1])^2/Dd[k]
        end
        convert(Float64,ad)<0 ?   Pd=Pd-ad : Qd=Qd-ad
        bd=(Pd+Qd)/(wzd*wzd)
        b=convert(Float64,bd)
    end

    if i<A.i
        # This reduces memory allocation (4 vectors) and is 30-40% faster
        for k=n-1:-1:A.i-1
            D[k+1]=D[k]
            z[k+1]=z[k]
        end
        D[A.i-1]=zero(T)
        z[A.i-1]=wz
        SymArrow(D,z,b,i),Kb,Kz,Qout
    else
        for k=n-1:-1:A.i
            D[k+1]=D[k]
            z[k+1]=z[k]
        end
        D[A.i]=zero(T)
        z[A.i]=wz
        SymArrow(D,z,b,i+1),Kb,Kz,Qout
    end
end # inv

function inv!(B::SymArrow{T},A::SymArrow{T},i::Integer,τ::Vector{Float64}=[1e3;10.0*length(A.D)]) where T
    # In-place version!
    # COMPUTES: inverse of a SymArrow matrix A, inv(A-A.D[i]*I) which is again SymArrow
    # uses higher precision to compute top of the arrow element accurately, if
    # needed.
    # τ=[tolb,tolz] are tolerances, usually [1e3, 10*n]
    # [0.0,0.0] forces DoubleDouble, [1e50,1e50] would never use it
    # RETURNS:  B=SymArrow(D,z,b,i), Kb, Kz, Qout
    # Kb - condition Kb, Kz - condition Kz, Qout = 1 / 0 - double was / was not used

    n=length(A.D)
    # D=Array{T}(undef,n)
    # z=Array{T}(undef,n)
    wz=one(T)/A.z[i]
    σ=A.D[i]
    for k=1:i-1
        B.D[k]=one(T)/(A.D[k]-σ)
        B.z[k]=-A.z[k]*B.D[k]*wz
    end
    for k=i+1:n
        B.D[k-1]=one(T)/(A.D[k]-σ)
        B.z[k-1]=-A.z[k]*B.D[k-1]*wz
    end

    # D=[1./(A.D[1:i-1]-shift),1./(A.D[i+1:end]-σ)]
    # wz=1/A.z[i]
    # z=-[A.z[1:i-1], A.z[i+1:end]].*D1*wz1
    a=A.a-A.D[i]

    # compute the sum in a plain loop
    P=zero(T)
    Q=zero(T)
    Qout=0
    nn=n-1
    for k=1:i-1
        B.D[k]>0.0 ? P+=A.z[k]^2*B.D[k] : Q+=A.z[k]^2*B.D[k]
    end
    for k=i:nn
        B.D[k]>0.0 ? P+=A.z[k+1]^2*B.D[k] : Q+=A.z[k+1]^2*B.D[k]
    end
    a<0 ? P-=a : Q-=a

    Kb=(P-Q)/abs(P+Q)
    Kz=maximum(abs,A.z)*abs(wz)

    if Kb<τ[1] ||  Kz<τ[2]
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
        σ₁=map(Type,A.D[i])
        Dd=[[Type(A.D[k])-σ₁ for k=1:i-1];
            [Type(A.D[k])-σ₁ for k=i+1:length(A.D)]]
        wzd=Type(A.z[i])
        ad=Type(A.a)-σ₁

        Pd,Qd=map(Type,(0.0,0.0))
        for k=1:i-1
            convert(Float64,Dd[k])>0.0 ? Pd+=Type(A.z[k])^2/Dd[k] :
            Qd+=Type(A.z[k])^2/Dd[k]
        end
        for k=i:nn
            convert(Float64,Dd[k])>0.0 ? Pd+=Type(A.z[k+1])^2/Dd[k] :
            Qd+=Type(A.z[k+1])^2/Dd[k]
        end
        convert(Float64,ad)<0 ?   Pd-=ad : Qd-=ad
        bd=(Pd+Qd)/(wzd*wzd)
        b=convert(Float64,bd)
    end

    if i<A.i
        # This reduces memory allocation (4 vectors) and is 30-40% faster
        for k=n-1:-1:A.i-1
            B.D[k+1]=B.D[k]
            B.z[k+1]=B.z[k]
        end
        B.D[A.i-1]=zero(T)
        B.z[A.i-1]=wz
        B.i=i
        B.a=b
        # SymArrow(D,z,b,i),Kb,Kz,Qout
        return Kb,Kz,Qout
    else
        for k=n-1:-1:A.i
            B.D[k+1]=B.D[k]
            B.z[k+1]=B.z[k]
        end
        B.D[A.i]=zero(T)
        B.z[A.i]=wz
        B.i=i+1
        B.a=b
        # SymArrow(D,z,b,i+1),Kb,Kz,Qout
        return Kb,Kz,Qout
    end
end # inv!

function inv(A::SymArrow{T}, σ::Float64, τ::Float64=1.0e3) where T
    # COMPUTES: inverse of the shifted SymArrow A, inv(A-σ*I) which is SymDPR1
    # uses DoubleDouble to compute rho accurately, if needed.
    # τ is tolerance, usually 1e3,  0.0 forces Double, 1e50 would never use it
    # RETURNS: SymDPR1(D,u,ρ), Kρ, Qout
    # Kρ - condition Kρ, Qout = 1 / 0 - Double was / was not used

    n=length(A.D)
    D=Array{T}(undef,n+1)
    u=Array{T}(undef,n+1)
    for k=1:A.i-1
        D[k]=one(T)/(A.D[k]-σ)
        u[k]=A.z[k]*D[k]
    end
    D[A.i]=zero(T)
    u[A.i]=-one(T)
    for k=A.i:n
        D[k+1]=one(T)/(A.D[k]-σ)
        u[k+1]=A.z[k]*D[k+1]
    end

    # compute rho and Krho
    #--- compute the sum in a plain loop
    a=A.a-σ
    P=zero(T)
    Q=zero(T)
    Qout=0
    for k=1:A.i-1
        D[k]>0.0 ? P+=A.z[k]^2*D[k] : Q+=A.z[k]^2*D[k]
    end
    for k=A.i:n
        D[k+1]>0.0 ? P+=A.z[k]^2*D[k+1] : Q+=A.z[k]^2*D[k+1]
    end
    a<0 ? P-=a : Q-=a

    # Condition of rho
    Kρ=(P-Q)/abs(P+Q)
    if Kρ<τ
        ρ=-1.0/(P+Q)
    else  # recompute in Double
        Qout=1
        Pd,Qd=map(Double,(0.0,0.0))
        σ₁=Double(σ)
        ad=Double(A.a)-σ₁
        for k=1:A.i-1
            D[k]>0.0 ? Pd+=Double(A.z[k])^2/(Double(A.D[k])-σ₁) : Qd+=Double(A.z[k])^2/(Double(A.D[k])-σ₁)
        end
        for k=A.i:n
            D[k+1]>0.0 ? Pd+=Double(A.z[k])^2/(Double(A.D[k])-σ₁) : Qd+=Double(A.z[k])^2/(Double(A.D[k])-σ₁)
        end
        ad.hi+ad.lo<0 ? Pd-=ad : Qd-=ad
        if Pd+Qd!=Double(0.0,0.0)
            Kρ=Float64((Pd-Qd)/abs(Pd+Qd))
            r=Double(1.0)/(Pd+Qd)
            ρ=-(r.hi+r.lo)
        else
        # Here we need quadruple working precision. We are using BigFloat.
        # Example of a matrix where this is neeed, courtesy of Stan Eisenstat, is:
        # A=SymArrow([1+eps(), 1-eps(), -1+2*eps(), -1-2*eps()],[2,2,1,1.0],6.0,5)
            Qout=100
            Pd,Qd=map(BigFloat,(0.0,0.0))
            σd=BigFloat(σ)
            ad=BigFloat(A.a)-σ₁
            for k=1:A.i-1
                D[k]>0.0 ? Pd+=BigFloat(A.z[k])^2/(BigFloat(A.D[k])-σ₁) : Qd+=BigFloat(A.z[k])^2/(BigFloat(A.D[k])-σ₁)
            end
            for k=A.i:n
                D[k]>0.0 ? Pd+=BigFloat(A.z[k])^2/(BigFloat(A.D[k])-σ₁) : Qd+=BigFloat(A.z[k])^2/(BigFloat(A.D[k])-σ₁)
            end
            ad<0 ? Pd-=ad : Qd-=ad
            Kρ=Float64((Pd-Qd)/abs(Pd+Qd))
            r=BigFloat(-1.0)/(Pd+Qd)
            ρ=Float64(r)
        end
    end

    # returns the following
    SymDPR1(D,u,ρ), Kρ, Qout
end # inv

function inv(A::SymArrow{T}, σ::Double) where T
    # COMPUTES: inverse of the shifted SymArrow A, inv(A-shift*I), which is a SymDPR1
    # here shift is Double so it uses Double to compute everything
    # RETURNS: SymDPR1(D₁,u₁,ρ), Qout
    # Qout = 1 on exit meaning Double was used

    n=length(A.D)
    D=Array{Double}(undef,n+1)
    u=Array{Double}(undef,n+1)
    for k=1:n
        D[k]=Double(A.D[k])-σ
    end
    u[1:n]=map(Double,A.z)
    a=map(Double,A.a)-σ
    i=A.i
    oned=Double(1.0)
    zerod=Double(0.0)
    for k=1:i-1
        u[k]=u[k]/D[k]
        D[k]=oned/D[k]
    end
    for k=i:n
        u[k+1]=u[k]/D[k]
        D[k+1]=oned/D[k+1]
    end
    D[i]=zerod
    u[i]=Double(-1.0)

    # compute rho and Kρ
    P,Q=zerod,zerod
    Qout=1
    for k=1:i-1
        D[k].hi > 0.0 ? P+=Double(A.z[k])*u[k] : Q+=Double(A.z[k])*u[k]
    end
    for k=i:n
        D[k+1].hi >0.0 ? P+=Double(A.z[k])*u[k+1] : Q+=Double(A.z[k])*u[k+1]
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

mutable struct Info
    Sind::Int
    Kb::Float64
    Kz::Float64
    Kν::Float64
    Kρ::Float64
    Qout::Int
end


function  eigen( A::SymArrow{T}, Ainv::SymArrow{T}, k::Integer,
    τ::Vector{Float64}=[1e3,10.0*length(A.D),1e3,1e3,1e3]) where T
    # COMPUTES: k-th eigenpair of an ordered irreducible SymArrow
    # A = [Diagonal(D) z; z' alpha]
    # τ=[tolb,tolz,tolnu,tolrho,tollambda] = [1e3,10.0*n,1e3,1e3,1e3]
    # RETURNS: λ,v, Info
    # where
    # λ - k-th eigenvalue in descending order
    # v - λ's normalized eigenvector
    # Info = Sind, Kb, Kz, Kν, Kρ, Qout
    # Sind = shift index i for the k-th eigenvalue
    # Kb, Kz, Kν, Kρ - condition numbers
    # Qout = 1 / 0 - Double was / was not used

    # Set the dimension
    n = length(A.D) + 1
    # Set all conditions initially to zero
    Kb,Kz,Kν,Kρ=0.0,0.0,0.0,0.0
    Qout=0
    v=zeros(T,n)
    # Kz is former κ_ν
    # Determine the shift σ, the shift index i, and whether λ
    # is on the left or the right side of the nearest pole
    # Exterior eigenvalues (k = 1 or k = n):
    if k == 1
        σ,i,side = A.D[1],1,'R'
    elseif k==n
        σ,i,side = A.D[n-1],n-1,'L'
    else
        # Interior eigenvalues (k in (2,...n-1) ):
        Dtemp = A.D.-A.D[k]
        atemp = A.a-A.D[k]
        middle = Dtemp[k-1]/2.0
        Fmiddle = (atemp-middle)-sum(A.z.^2 ./(Dtemp.-middle))
        σ,i,side = Fmiddle < 0.0 ? (A.D[k],k,'R') : (A.D[k-1],k-1,'L')
    end

    # Compute the inverse of the shifted matrix, A_i^(-1), Kb and Kz
    # Ainv,Kb,Kz,Qout = inv(A,i)
    Kb,Kz,Qout = inv!(Ainv,A,i)

    # Compute the eigenvalue of the inverse shifted matrix
    ν = bisect( Ainv,side )

    if abs(ν)==Inf
        # this is nonstandard
        # Deflation in aheig (ν=Inf)
        v[i]=1.0
        λ=σ
    else
        # standard case, full computation
        # ν₁ is the F- or 1-norm of the inverse of the shifted matrix
        # ν₁=maximum([sum(abs.(Ainv.z))+abs(Ainv.a), maximum(abs,Ainv.D)+abs.(Ainv.z))])
        ν₁=0.0
        for l=1:n-1
            ν₁=max(ν₁,abs(Ainv.D[l])+abs(Ainv.z[l]))
        end
        ν₁=max(sum(abs,Ainv.z)+abs(Ainv.a), ν₁)
        Kν=ν₁/abs(ν)
        while  Kν>τ[3]
            # Remedies according to Remark 3 - we shift between original
            # eigenvalues and compute DPR1 matrix
            # 1/ν₁+σ, 1/ν+σ
            # println("Remedy 3 ")
            ν = side=='R' ? abs(ν) : -abs(ν)
            ν₁=-sign(ν)*ν₁
            σ₁=(ν₁+ν)/(2.0*ν*ν₁)+σ
            if count(isequal(σ₁),A.D)>0  # we came back with a pole
                # recompute σᵥ more accurately with Dekker
                σᵥ=(Double(ν₁)+Double(ν))/(Double(2.0)*Double(ν)*Double(ν₁))+Double(σ)
                # Compute the inverse of the shifted arrowhead (DPR1)
                Ainv,Qout₁=inv(A,σᵥ) # Ainv is Float64
                ν₁=bisect(Ainv,side)
                μ₁ = 1.0/ν₁
                D0=map(Double,A.D)-σᵥ
                D1=[D0[l].hi+D0[l].lo for l=1:length(D0)]
                if count(isequal(μ₁),D1)>0
                    ind=findall(isequal(μ₁),D1)
                    v=zeros(T,n)
                    v[ind]=1.0
                else
                    v=[ A.z./(D1.-μ₁);-1.0]
                end
                # Shift the eigenvalue back in Double
                lam = Double(1.0)/Double(ν₁)+σᵥ
                # Return this
                λ=lam.hi+lam.lo
                normalize!(v)
                κ=Info(i,Kb,Kz,Kν,Kρ,Qout)
                return λ,v,κ
            else
                # Compute the inverse of the shifted arrowhead (DPR1)
                Ainv₁, Kρ,Qout₁=inv(A,σ₁,τ[4]) # Ainv is Float64
                # Compute the eigenvalue by bisect for DPR1
	            # Note: instead of bisection could use dlaed4 (need a wrapper) but
	            # it is not faster. There norm(u)==1
                ν= bisect(Ainv₁,side)
                if side=='R' && A.D[1]>0.0 && ν<0.0
                    ν=bisect(Ainv₁,'L')
                end
                μ=1.0/ν
                ν₁=maximum(abs,Ainv₁.D)+abs(Ainv₁.r)*dot(Ainv₁.u,Ainv₁.u)
                Kν=ν₁/abs(ν)
                σ=σ₁
            end
            Qout=Qout+2*Qout₁
        end

        # Accuracy is fine, compute the eigenvector
        μ = 1.0/ν
        for l=1:n-1
            v[l] = A.z[l]/((A.D[l]-σ)-μ)
        end
        v[n]=-1.0
        λ=μ+σ
        normalize!(v)

        # Remedy according to Remark 1 - we recompute the the eigenvalue
        # near zero from the inverse of the original matrix (a DPR1 matrix).
        if (abs(A.D[i])+abs(1.0/ν))/abs(λ)>τ[5]
            if k==1 && A.D[1]<0.0 || k==n && A.D[n-1]>0.0 || i<n-1 && side=='L' &&
                sign(A.D[i])+sign(A.D[i+1])==0 || i>1 &&
                side=='R' && sign(A.D[i])+sign(A.D[i-1])==0
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
end # eigen(k)

function eigen(A::SymArrow{T}, τ::Vector{Float64}=[1e3,10.0*length(A.D),1e3,1e3,1e3]) where T
    # COMPUTES: all eigenvalues and eigenvectors of a real symmetric SymArrow
    # A = [diag(D) z;z' alpha] (notice, here we assume A.i==n)
    # τ = [tolb,tolz,tolnu,tolrho,tollambda] = [1e3,10.0*n,1e3,1e3,1e3] or similar
    # RETURNS: Eigen(Λ,U), κ
    # where
    # Λ = eigenvalues in decreasing order, U = eigenvectors,
    # κ[k]=Info(Sind[k], Kb[k], Kz[k], Kν[k], Kρ[k], Qout[k])
    # Sind[k] - shift index i for the k-th eigenvalue
    # Kb, Kz, Kν, Kρ [k] - respective conditions for the k-th eigenvalue
    # Qout[k] = 1 / 0 - Double was / was not used when computing k-th eigenvalue

    n=length(A.D)+1
    n0=n
    # Ordering the matrix
    is=sortperm(A.D,rev=true)
    D=A.D[is]
    z=A.z[is]
    U=Array{T}(I,n,n)
    Λ=zeros(T,n)
    κ=[Info(0, 0.0, 0.0, 0.0, 0.0, 0) for i=1:n]
    # Quick return for 1x1
    if n==1
        return Eigen([A.a],U),κ
    end

    #  test for deflation in z
    z0=findall(iszero,z)
    zx=findall(!iszero,z)
    if isempty(zx)  # nothing to do
        Λ=[A.D;A.a]
        isΛ=sortperm(Λ,rev=true)
        Λ=Λ[isΛ]
        U=view(U,:,isΛ) # U[:,isΛ]
        return Eigen(Λ,U),κ
    end
    if !isempty(z0)
        Λ[z0]=D[z0]
        D=D[zx]
        z=z[zx]
        if isempty(z)
            Λ[n]=A.a
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
            Λ[zx[g0[l]+1]]=D[g0[l]+1]
        end
        # remains
        gx=[0;gx].+1
        nn=length(gx)
        zxx=zx[[gx;n]]
        Axx=SymArrow(D[gx],z[gx],A.a,nn+1)
        Bxx=[deepcopy(Axx) for i=1:Threads.nthreads()]
        Threads.@threads for k=1:nn+1
            tid=Threads.threadid()
            Λ[zxx[k]],U[zxx,zxx[k]],κ[zxx[k]] = eigen(Axx,Bxx[tid],k)
        end
        for l=1:lg0
            # U=R[l][1]'*U
            lmul!(R[l][1],U)
        end
    else
        # No deflation in D
        Ax=SymArrow(D,z,A.a,n)
        Bx=[deepcopy(Ax) for i=1:Threads.nthreads()]
        Threads.@threads for k=1:n
            tid=Threads.threadid()
            Λ[zx[k]],U[zx,zx[k]],κ[zx[k]] = eigen(Ax,Bx[tid],k)
        end
    end
    # end
    # back premutation of vectors
    isi=sortperm(is)
    # must sort Λ once more
    es=sortperm(Λ,rev=true)
    Λ.=Λ[es]
    U=view(U,[isi[1:A.i-1];n0;isi[A.i:n0-1]],es)
    # U.=U[[isi[1:A.i-1];n0;isi[A.i:n0-1]],es]
    κ.=κ[es]
    # Return this
    return Eigen(Λ,U),κ
end # eigen (all)

function bisect(A::SymArrow{T}, side::Char) where T
    # COMPUTES: the leftmost (for side='L') or the rightmost (for side='R') eigenvalue
    # of a SymArrow A = [diag (D) z; z'] by bisection.
    # RETURNS: the eigenvalue

    # Determine the starting interval for bisection, [left; right]
    # left, right = side == 'L' ? {minimum([A.D-abs(A.z),A.a-sum(abs(A.z))]), minimum(A.D)} :
    #   {maximum(A.D),maximum([A.D+abs.(A.z),A.a+sum(abs,A.z)])}
    z = abs.(A.z)
    if side == 'L'
        left  = minimum(A.D .- z)
        left  = min(left, A.a - sum(z))
        right = minimum(A.D)
    else
        left  = maximum(A.D)
        right = maximum(A.D .+ z)
        right = max(right, A.a + sum(z))
    end

    # Bisection
    middle = (left + right) / convert(T,2)
    z.^= 2
    count, n = 0, length(A.D)
    while (right-left) > 2.0 * eps() * max(abs(left), abs(right))
        # in @time 50% of the time was garbage collection. The following line
        # assigns new vector every time it is called, so it is much better in the
        # loop?? Original timing were 30 secs for n=4000, 2.17 sec for n=1000

        # Fmiddle = A.a-middle-sum(z./(A.D-middle))
        Fmiddle = zero(T)
        for k=1:n
            Fmiddle += z[k] / (A.D[k] - middle)
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
