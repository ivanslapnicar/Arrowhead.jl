function tdc(T::SymTridiagonal{S}) where S
    # Tridiagonal divide-and-conquer using arrowhead matrices
    # COMPUTES: eigenvectors U and eigenvalues Λ of a SymTridiagonal matrix T,
    # T = U*Diagonal(Λ)*U'
    # RETURNS: Eigen(Λ, U)
    n=length(T.dv)
    # set the tolerances for eigen
    E=Eigen(Vector{Float64}(undef,1),Array{Float64}(undef,1,1))
    τ=[1e2,1e2,1e2,1e2,1e2]
    if n==1
        E=Eigen([T.dv[1]],ones(1,1))
    elseif n==2
        E,info=eigen(SymArrow([T.dv[1]],[T.ev[1]],T.dv[2],2),τ)
    else
        k=div(n,2)
        T₁=SymTridiagonal(T.dv[1:k],T.ev[1:k-1])
        T₂=SymTridiagonal(T.dv[k+2:end],T.ev[k+2:end])
        E₁=tdc(T₁)
        E₂=tdc(T₂)
        # form arrowhead
        D=[E₁.values;E₂.values]
        z=[E₁.vectors'[:,end]*T.ev[k]; E₂.vectors'[:,1]*T.ev[k+1]]
        a=T.dv[k+1]
        A=SymArrow(D,z,a,k+1)
        E,info=eigen(A,τ)
        E.vectors[1:k,:].=E₁.vectors*E.vectors[1:k,:]
        E.vectors[k+2:end,:].=E₂.vectors*E.vectors[k+2:end,:]
    end
    return E
end

# roots
function rootsah(pol::Union{Poly{Float32},Poly{Float64},Poly{Int32},Poly{Int64}}, D::Vector{Float64})
    # COMPUTES: the roots of polynomials with all distinct real roots.
    # The computation is forward stable. The program uses SymArrow (arrowhead) companion matrix and
    # corresponding eig routine
    # D are barycentric coordinates - elements od D must interpolate the roots of P,
    # for example
    #              D=roots(polyder(pol))
    # RETURNS: roots E

    # Type is Float64
    T = Float64

    # τ = [1e2,1e2,1e2,1e2,1e2] or similar is the vector of tolerances for eig
    τ=[1e2,1e2,1e2,1e2,1e2]
    Dm=map(T,D)
    Dm=sort(Dm,rev=true)
    D=sort(D,rev=true)
    p=map(T,[pol[i] for i=0:1:length(pol)-1])
    p=p[end:-1:1]

    # Compute z of the arrowhead
    # First we compute values s=p(D) we use Horner sheme with Double
    n=length(p)-1
    pD=map(Double,p)
    DD=map(Double,Dm)
    oneD=Double(one(T))
    s=Array{Double{T}}(undef,n-1)
    for i=1:n-1
        s[i]=oneD
    end
    s=pD[1]*s
    for i=2:n+1
        r=s.*DD
        s=r.+pD[i]
    end

    # Compute t's
    t=Array{Double{T}}(undef,n-1)
    for j=1:n-1
        h=oneD
        for i=1:j-1
            g=DD[j]-DD[i]
            h=h*g
        end
        for i=j+1:n-1
            g=DD[j]-DD[i]
            h=h*g
        end
        t[j]=h
    end
    # Compute αD
    αD=pD[2]/pD[1]
    for i=1:n-1
        αD+=DD[i]
    end
    # αD=-αD
    αD=Double(-αD.hi,-αD.lo)
    #  Compute z
    zD=Array{Double{T}}(undef,n-1)
    for i=1:n-1
        #         zD[i]=sqrt((-s[i])/(t[i]*pD[1]))
        zD[i]=sqrt(Double(-s[i].hi,-s[i].lo)/(t[i]*pD[1]))
    end

    z=Array{Float64}(undef,n-1)
    for i=1:n-1
        z[i]=Float64(zD[i].lo+zD[i].hi)
    end
    α=Float64(αD.hi+αD.lo)
    A=SymArrow(D,z,α,n)
    E=zeros(Float64,n)
    Qout=zeros(Int64,n)
    for k=1:n
        E[k], Qout[k] = eigen( A,zD,αD,k,τ )
    end
    E, Qout
end

#------------------------
#--------  Functions

#-------- Inverses

function inv(A::SymArrow{T},zD::Vector{Double{T}}, αD::Double{T},i::Int64,τ::Vector{Float64}=[1e3, 10*length(A.D)]) where T
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
    elseif Kb==Inf
       b=0.0
    else  # recompute in Double
        Qout=1
        AD=map(Double{T},A.D)
        σd=AD[i]
        Dd=[Double{T}[AD[k]-σd for k=1:i-1];
            Double{T}[AD[k]-σd for
                            k=i+1:length(A.D)]]
        wzd=zD[i] # Double(A.z[i])
        ad=αD-σd
        Pd,Qd=map(Double{T},(0.0,0.0))
        for k=1:i-1
            Dd[k].hi>0.0 ? Pd+=zD[k]^2/Dd[k] :
            Qd+=zD[k]^2/Dd[k]
        end
        for k=i:nn
            Dd[k].hi>0.0 ? Pd+=zD[k+1]^2/Dd[k] :
            Qd+=zD[k+1]^2/Dd[k]
        end
        ad.hi<0 ?   Pd-=ad : Qd-=ad

        bd=(Pd+Qd)/(wzd*wzd)
        b=Float64(bd.hi+bd.lo)
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

function inv(A::SymArrow{T}, zD::Vector{Double{T}}, αD::Double{T}, σ::Float64, τ::Float64=1.0e3) where T
    # COMPUTES: inverse of the shifted SymArrow A, inv(A-σ*I) which is SymDPR1
    # uses DoubleDouble to compute ρ accurately, if needed.
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

    # compute ρ and Kρ
    #--- compute the sum in a plain loop
    a=A.a-σ
    P=zero(T)
    Q=zero(T)
    Qout=0
    for k=1:A.i-1
        D[k]>0.0 ? P+=A.z[k]^2*D[k] : Q+=A.z[k]^2*D[k]
    end
    for k=A.i:n
        D[k+1]>0.0 ? P+=A.z[k]^2*D[k+1] : Q+=!=()A.z[k]^2*D[k+1]
    end
    a<0 ? P-=a : Q-=a

    # Condition of ρ
    Kρ=(P-Q)/abs(P+Q)
    if Kρ<τ
        ρ=-1.0/(P+Q)
    else  # recompute in Double
        Qout=1
        Pd,Qd=map(Double{T},(0.0,0.0))
        AD=map(Double{T},A.D)
        σd=map(Double{T},σ)
        ad=αD-σd
        for k=1:A.i-1
            D[k]>0.0 ? Pd+=zD[k]^2/(AD[k]-σd) : Qd+=zD[k]^2/(AD[k]-σd)
        end
        for k=A.i:n
            D[k+1]>0.0 ? Pd+=zD[k]^2/(AD[k]-σd) : Qd+=zD[k]^2/(A.D[k]-σd)
        end
        ad.hi+ad.lo<0 ? Pd=Pd-ad : Qd=Qd-ad
        r=map(Double{T},1.0)/(Pd+Qd)
        ρ=Float64(-(r.hi+r.lo))
    end

    # returns the following
    SymDPR1(D,u,ρ), Kρ, Qout
end # inv

function inv(A::SymArrow{T}, zD::Vector{Double{T}}, αD::Double{T}, σ::Double{Float64}) where T
    # COMPUTES: inverse of the σed SymArrow A, inv(A-σ*I), which is a SymDPR1
    # here σ is Double so it uses Double to compute everything
    # RETURNS: SymDPR1(D1,u1,ρ), Qout
    # Qout = 1 on exit meaning Double was used

    n=length(A.D)
    D=Array{Double{T}}(undef,n+1)
    u=Array{Double{T}}(undef,n+1)
    σd=map(Double{T},σ)
    for k=1:n
        D[k]=map(Double{T},A.D[k])-σd
    end
    u[1:n]=zD
    a=αD-σd
    i=A.i

    oned=Double(one(T))
    zerod=Double(zero(T))
    for k=1:i-1
        u[k]=u[k]/D[k]
        D[k]=oned/D[k]
    end
    for k=i:n
        u[k+1]=u[k]/D[k]
        D[k+1]=oned/D[k+1]
    end
    D[i]=zerod
    u[i]=-oned

    # compute ρ and Kρ
    #--- compute the sum in a plain loop
    P,Q=zerod,zerod
    Qout=1
    for k=1:i-1
        D[k].hi > 0.0 ? P+=zD[k]*u[k] : Q+=zD[k]*u[k]
    end
    for k=i:n
        D[k+1].hi >0.0 ? P+=zD[k]*u[k+1] : Q+=zD[k]*u[k+1]
    end
    a.hi+a.lo<0.0  ? P-=a : Q-=a
    r=oned/(P+Q)
    ρ=Float64(-(r.hi+r.lo))

    # returns the following
    D₁=Array{T}(undef,n+1)
    u₁=Array{T}(undef,n+1)
    for k=1:n+1
        D₁[k]=Float64(D[k].hi+D[k].lo)
    end
    for k=1:n+1
        u₁[k]=Float64(u[k].hi+u[k].lo)
    end
    SymDPR1(D₁,u₁,ρ), Qout
end # inv

function  eigen( A::SymArrow{T},zD::Vector{Double{T}},αD::Double{T},k::Int64,
    τ::Vector{Float64}=[1e3,10.0*length(A.D),1e3,1e3,1e3]) where T
    # COMPUTES: k-th eigenvalue of an ordered irreducible SymArrow
    # A = [Diagonal(D) z; z' α]
    # Specially designed to be used in the polynomial rootfinder rootsah !!!!
    # τ=[tolb,tolz,tolnu,tolrho,tollambda] = [1e3,10.0*n,1e3,1e3,1e3]
    # RETURNS: λ
    # λ - k-th eigenvalue
    # Kb, Kz, Kν, Kρ - condition numbers
    # Qout = 1 / 0 - Double was / was not used
    # If Qout>0, quad precision was used
    # i was the shift index

    # Set the dimension
    n = length(A.D) + 1
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
    elseif k==n
        σ,i,side = A.D[n-1],n-1,'L'
    else
        # Interior eigenvalues (k in (2,...n-1) ):
        # We need to compute this in Double
        Dd=map(Double,A.D)
        middle=(Dd[k-1]-Dd[k])/Double(2.0)
        Fmiddle=Double(zero(T))
        for l=1:n-1
            Fmiddle+=zD[l]^2/((Dd[l]-Dd[k])-middle)
        end
        Fmiddle=((αD-Dd[k])-middle)-Fmiddle
        σ,i,side = Fmiddle.hi+Fmiddle.lo < 0.0 ? (A.D[k],k,'R') : (A.D[k-1],k-1,'L')
    end

    # Compute the inverse of the shifted matrix, A_i^(-1), Kb and Kz
    Ainv,Kb,Kz,Qout = inv(A,zD,αD,i)

    # Compute the eigenvalue of the inverse shifted matrix
    ν = bisect( Ainv,side )

    #  ν=fastzero([invD1; 0; invD2], [w1;wz;w2], b, side); # implement later
    #  [ν-ν₁]/ν, [ν-nueva]/ν, pause  # and compare
    if abs(ν)==Inf
        # this is nonstandard
        # Deflation in aheig (ν=Inf)
        v[i]=1.0
        λ=σ
    else
        # standard case, full computation
        # ν₁ is the F- or 1-norm of the inverse of the shifted matrix
        # ν₁0=maximum([sum(abs,Ainv.z)+abs(Ainv.a), maximum(abs.(Ainv.D)+abs.(Ainv.z))])
        ν₁=0.0
        for l=1:n-1
            ν₁=max(ν₁,abs(Ainv.D[l])+abs(Ainv.z[l]))
        end
        ν₁=max(sum(abs,Ainv.z)+abs(Ainv.a), ν₁)
        Kν=ν₁/abs(ν)
        if Kν<τ[3]
            # Accuracy is fine, compute the eigenvector
            μ = 1.0/ν
            # v=[ A.z./((A.D-σ)-μ);-1.0]
            for l=1:n-1
            	v[l] = A.z[l]/((A.D[l]-σ)-μ)
            end
            v[n]=-1.0
            λ=μ+σ
            normalize!(v)
        else
            # Remedies according to Remark 3 - we shift between original
            # eigenvalues and compute DPR1 matrix
            # 1/ν₁+σ, 1/ν+σ
            # println(k, " Remedy 3 ")
            ν = side=='R' ? abs(ν) : -abs(ν)
            ν₁=-sign(ν)*ν₁
            σ₁=(ν₁+ν)/(2.0*ν*ν₁)+σ
            if count(isequal(σ₁),A.D)>0 # we came back with a pole
                # recompute σᵥ more accurately according with dekker
                σᵥ=(Double(ν₁)+Double(ν))/(Double(2.0)*Double(ν)*Double(ν₁))+Double(σ)
                # Compute the inverse of the shifted arrowhead (DPR1)
                Ainv,Qout₁=inv(A,σᵥ) # Ainv is Float64
                ν₁=bisect(Ainv,side)
                μ₁ = 1.0/ν₁
                D0=map(Double,A.D)-σᵥ
                D1=[D0[l].hi+D0[l].lo for l=1:length(D0)]
                if count(isequal(μ₁),D1)>0
                    ind=findall(isequal(μ₁),D1)
                    v=zeros(n)
                    v[ind]=1.0
                else
                    v=[ A.z./(D1-μ₁);-1.0]
                end
                # Shift the eigenvalue back in Double
                lam = Double(1.0)/Double(ν₁)+σᵥ
                # Return this
                λ = lam.hi+lam.lo
                normalize!(v)
            else
                # Compute the inverse of the shifted arrowhead (DPR1)
                Ainv, Kρ,Qout₁=inv(A,σ₁,τ[4]) # Ainv is Float64
                # Compute the eigenvalue by bisect for DPR1
	            # Note: instead of bisection could use dlaed4 (need a wrapper) but
	            # it is not faster. There norm(u)==1
                ν₁= bisect(Ainv,side)
                μ₁=1.0/ν₁
                # standard v
                v=[ A.z./((A.D-σ₁)-μ₁);-1.0]
                # Return this - shift the eigenvalue back and normalize the vector
                λ = μ₁+σ₁
                normalize!(v)
            end
            Qout==1 && (Qout=Qout+2)
        end

        # Remedy according to Remark 1 - we recompute the the eigenvalue
        # near zero from the inverse of the original matrix (a DPR1 matrix).
        if (abs(A.D[i])+abs(1.0/ν))/abs(λ)>τ[5]
            if k==1 && A.D[1]<0.0 || k==n && A.D[n-1]>0.0 || i<n-1 && side=='L' && sign(A.D[i])+sign(A.D[i+1])==0 ||
                side=='R' && sign(A.D[i])+sign(A.D[i-1])==0
                # println(k, " Remedy 1 ")
                # Compute the inverse of the original arrowhead (DPR1)
                Ainv,Kρ,Qout₁ = inv(A,0.0,τ[4]) # Ainv is Float64
                Qout==1 && (Qout=Qout+4)
                if abs(Ainv.r)==Inf
                    λ=0.0
                else
                    # Here we do not need bisection. We compute the Rayleigh
                    # quotient by using already computed vectors which is
                    # componentwise accurate
                    ν₁=sum(v.^2 .*Ainv.D)+Ainv.r*sum(v.*Ainv.u)^2;
                    λ=1.0/ν₁
                end
            end
        end
    end

    # Return this
    λ, Qout
end # eigen (k)
