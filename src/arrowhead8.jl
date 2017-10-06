# rootsah for BigInt and BigFloat
function rootsah(pol::Union{Poly{BigInt},Poly{BigFloat}}, D::Vector{Float64})

    # COMPUTES: the roots of polynomials with all distinct real roots.
    # The computation is forward stable. The program uses SymArrow (arrowhead) companion matrix and
    # corresponding eig routine
    # D are barycentric coordinates - elements od D must interpolate the roots of P,
    # for example
    #              D=roots(polyder(pol))
    # RETURNS: roots E

    T = BigFloat
    
    # tols = [1e2,1e2,1e2,1e2,1e2] or similar is the vector of tolerances for eig
    tols=[1e2,1e2,1e2,1e2,1e2]

    Dm=map(T,D)
    Dm=sort(Dm,rev=true)
    D=sort(D,rev=true)
    
    p=map(T,[pol[i] for i=0:1:length(pol)-1])
    p=p[end:-1:1]

    # Compute z of the arrowhead
    # First we compute values s=p(D) we use Horner sheme with double the working precision
    setprecision(512)
    
    n=length(p)-1
    # pD=map(Double,p)
    # DD=map(Double,Dm)
    # oneD=Double(one(T))

    s=Array{BigFloat}(n-1)
    for i=1:n-1
        s[i]=one(T)
    end
    s=p[1]*s
    for i=2:n+1
        r=s.*Dm
        s=r+p[i]
    end

    # Compute t's
    t=Array{BigFloat}(n-1)
    for j=1:n-1
        h=one(T)
        for i=1:j-1
            g=Dm[j]-Dm[i]
            h=h*g
        end
        for i=j+1:n-1
            g=Dm[j]-Dm[i]
            h=h*g
        end
        t[j]=h
    end
    # Compute alphaD 
    alphaD=p[2]/p[1]
    for i=1:n-1
        alphaD+=Dm[i]
    end
    alphaD=-alphaD
    #  Compute z
    zD=Array{BigFloat}(n-1) 
    for i=1:n-1
        zD[i]=sqrt((-s[i])/(t[i]*p[1]))
    end
    
    z=Array{Float64}(n-1)
    for i=1:n-1
        z[i]=Float64(zD[i])
    end
    alpha=Float64(alphaD)
    A=SymArrow(D,z,alpha,n)
    
    E=zeros(Float64,n)
    Qout=zeros(Int64,n)
    
    for k=1:n
        E[k], Qout[k] = eig( A,zD,alphaD,k,tols )
    end

    # reset the bigfloat_precision
    setprecision(256)

    E, Qout
end

#------------------------
#--------  Functions

#-------- Inverses  

function inv{T}(A::SymArrow{T},zD::Vector{BigFloat}, alphaD::BigFloat,i::Int64,tols::Vector{Float64})

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
    else  # recompute in double the working precision
        Qout=1
        AD=map(BigFloat,A.D)
        shiftd=AD[i]
        #shiftd=map(Double,A.D[i])
        Dd=[[AD[k]-shiftd for k=1:i-1]; 
            [AD[k]-shiftd for k=i+1:length(A.D)]]
        wzd=zD[i] #
        ad=alphaD-shiftd

        Pd,Qd=map(BigFloat,(0.0,0.0))
        
        for k=1:i-1
            Dd[k] > 0.0 ? Pd+=zD[k]^2/Dd[k] :
            Qd+=zD[k]^2/Dd[k]
        end
        
        for k=i:nn
            Dd[k] > 0.0 ? Pd+=zD[k+1]^2/Dd[k] :
            Qd+=zD[k+1]^2/Dd[k]
        end 

        ad<0 ?   Pd=Pd-ad : Qd=Qd-ad 

        bd=(Pd+Qd)/(wzd*wzd)
        b=Float64(bd)
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

function inv{T}(A::SymArrow{T}, zD::Vector{BigFloat}, alphaD::BigFloat, shift::Float64, tolr::Float64)

    # COMPUTES: inverse of the shifted SymArrow A, inv(A-shift*I) which is SymDPR1
    # uses DoubleDouble to compute rho accurately, if needed. 
    # tolr is tolerance, usually 1e3,  0.0 forces Double, 1e50 would never use it
    # RETURNS: SymDPR1(D,u,rho), Krho, Qout
    # Krho - condition Krho, Qout = 1 / 0 - Double was / was not used 

    n=length(A.D)
    D=Array(T,n+1)
    u=Array(T,n+1)

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

    else  # recompute in double the working precision
        Qout=1
        Pd,Qd=map(BigFloat,(0.0,0.0))
        AD=map(BigFloat,A.D)
        shiftd=map(BigFloat,shift)
        ad=alphaD-shiftd
        #   Dd=[Double(A.D[k])-shiftd for k=1:length(A.D) ] 
        #   Dd=map(Double,A.D)-shiftd
        for k=1:A.i-1
            D[k]>0.0 ? Pd=Pd+zD[k]^2/(AD[k]-shiftd) : Qd=Qd+zD[k]^2/(AD[k]-shiftd)
        end
        for k=A.i:n
            D[k+1]>0.0 ? Pd=Pd+zD[k]^2/(AD[k]-shiftd) : Qd=Qd+zD[k]^2/(A.D[k]-shiftd)
        end

        # for k=1:length(A.D)
        #    Dd[k].hi>0.0 ? Pd=Pd+Double(A.z[k])^2/Dd[k] : Qd=Qd+Double(A.z[k])^2/Dd[k]
        # end
        ad<0 ? Pd=Pd-ad : Qd=Qd-ad
        r=map(BigFloat,1.0)/(Pd+Qd)
        rho=Float64(-r)
    end

    # returns the following
    SymDPR1(D,u,rho), Krho, Qout

end # inv

function inv{T}(A::SymArrow{T}, zD::Vector{BigFloat}, alphaD::BigFloat, shift::BigFloat)

    # COMPUTES: inverse of the shifted SymArrow A, inv(A-shift*I), which is a SymDPR1
    # here shift is BigFloat so it uses double the working precision to compute everything 
    # RETURNS: SymDPR1(D1,u1,rho), Qout
    # Qout = 1 on exit meaning Double was used

    n=length(A.D)
    D=Array(BigFloat,n+1)
    u=Array(BigFloat,n+1)

    # D[1:n]=map(Double,A.D)-shift
    for k=1:n
        D[k]=map(BigFloat,A.D[k])-shift
    end

    u[1:n]=zD
    a=alphaD-shift
    i=A.i

    oned=one(BigFloat)
    zerod=zero(BigFloat)

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

    # compute rho and Krho
    #--- compute the sum in a plain loop

    P,Q=zerod,zerod
    Qout=1

    for k=1:i-1
        # D[k].hi+D[k].lo > 0.0  # .lo cannot change sign, I think
        # D[k].hi > 0.0 ? P=P+Double(A.z[k])^2*D[k] : Q=Q+Double(A.z[k])^2*D[k]
        D[k] > 0.0 ? P=P+zD[k]*u[k] : Q=Q+zD[k]*u[k]
    end
    for k=i:n
        D[k+1] > 0.0 ? P=P+zD[k]*u[k+1] : Q=Q+zD[k]*u[k+1]
    end

    a<0.0  ? P=P-a : Q=Q-a
    r=oned/(P+Q)


    rho=Float64(-r)

    # returns the following
    D1=Array(T,n+1)
    u1=Array(T,n+1)
    for k=1:n+1
        D1[k]=Float64(D[k])
    end
    for k=1:n+1
        u1[k]=Float64(u[k])
    end

    SymDPR1(D1,u1,rho), Qout

end # inv


function  eig{T}( A::SymArrow{T},zD::Vector{BigFloat},alphaD::BigFloat,k::Int64,tols::Vector{Float64})

    # COMPUTES: k-th eigenvalue of an ordered irreducible SymArrow
    # A = [diag (D) z; z' alpha]
    # Specially designed to be used in the polynomial rootfinder rootsah !!!!
    # tols=[tolb,tolz,tolnu,tolrho,tollambda] = [1e3,10.0*n,1e3,1e3,1e3]
    # RETURNS: lambda
    # lambda - k-th eigenvalue
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
        # We need to compute Fmiddle carefully
        
        Dtemp = A.D-A.D[k]
        middle = Dtemp[k-1]/2.0

        P=zero(Float64)
        Q=zero(Float64)
        P=P+sum(A.z[1:k-1].^2./(Dtemp[1:k-1]-middle))
        Q=Q+sum(A.z[k:n-1].^2./(Dtemp[k:n-1]-middle))

        P1=map(BigFloat,P)
        Q1=map(BigFloat,Q)
        atemp = (alphaD-map(BigFloat,A.D[k]))-map(BigFloat,middle)
        atemp<0 ? P1=P1-atemp : Q1=Q1-atemp

        # Fmiddle = (atemp-middle)-sum(zD.^2./(Dtemp-middle))
        Fmiddle=-(P1+Q1)
        sigma,i,side = Fmiddle < zero(BigFloat) ? (A.D[k],k,'R') : (A.D[k-1],k-1,'L')
    end

    # Compute the inverse of the shifted matrix, A_i^(-1), Kb and Kz

    Ainv,Kb,Kz,Qout = inv(A,zD,alphaD,i,tols[1:2])

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
        # nu10=maximum([sum(abs(Ainv.z))+abs(Ainv.a), maximum(abs(Ainv.D)+abs(Ainv.z))])
        nu1=0.0
        for l=1:n-1
            nu1=max(nu1,abs(Ainv.D[l])+abs(Ainv.z[l]))
        end
        nu1=max(sum(abs,Ainv.z)+abs(Ainv.a), nu1)

        Knu=nu1/abs(nu)
        if Knu<tols[3]
            # Accuracy is fine, compute the eigenvector
            mu = 1.0/nu
            # v=[ A.z./((A.D-sigma)-mu);-1.0]
            for l=1:n-1
            	v[l] = A.z[l]/((A.D[l]-sigma)-mu)
            end 
            v[n]=-1.0
            lambda, v = mu+sigma, v/norm(v)
        else
            # Remedies according to Remark 3 - we shift between original
            # eigenvalues and compute DPR1 matrix
            # 1/nu1+sigma, 1/nu+sigma
            # println("Remedy 3 ")
            nu = side=='R' ? abs(nu) : -abs(nu)
            nu1=-sign(nu)*nu1
            sigma1=(nu1+nu)/(2.0*nu*nu1)+sigma

            if findfirst(A.D-sigma1,0.0)>0 # we came back with a pole
                # recompute sigmav more accurately according with dekker
                setprecision(512)
                sigmav=(nu1+nu)/(BigFloat(2.0)*nu*nu1)+sigma
                # Compute the inverse of the shifted arrowhead (DPR1)
                Ainv,Qout1=inv(A,sigmav) # Ainv is Float64
                nu1=bisect(Ainv,side) 
                mu1 = 1.0/nu1
                D0=map(BigFloat,A.D)-sigmav
                D1=[D0[l] for l=1:length(D0)]
                if findfirst(D1-mu1,0.0)>0 
                    ind=find(D1-mu1==0.0)
                    v=zeros(n)
                    v[ind]=1.0
                else
                    v=[ A.z./(D1-mu1);-1.0]
                end
                # Shift the eigenvalue back in Double
                lam = BigFloat(1.0)/nu1+sigmav
                # Return this
                lambda,v = lam, v/norm(v)       
            else
                # Compute the inverse of the shifted arrowhead (DPR1)
                Ainv, Krho,Qout1=inv(A,sigma1,tols[4]) # Ainv is Float64
                # Compute the eigenvalue by bisect for DPR1
	            # Note: instead of bisection could use dlaed4 (need a wrapper) but
	            # it is not faster. There norm(u)==1
                nu1= bisect(Ainv,side)
                mu1=1.0/nu1 
                # standard v
                v=[ A.z./((A.D-sigma1)-mu1);-1.0]
                # Return this - shift the eigenvalue back and normalize the vector
                lambda, v = mu1+sigma1, v/norm(v)
            end
            Qout==1 && (Qout=Qout+2)
        end

        # Remedy according to Remark 1 - we recompute the the eigenvalue
        # near zero from the inverse of the original matrix (a DPR1 matrix).  
        if (abs(A.D[i])+abs(1.0/nu))/abs(lambda)>tols[5]

            if k==1 && A.D[1]<0.0 || k==n && A.D[n-1]>0.0 || i<n-1 && side=='L' && sign(A.D[i])+sign(A.D[i+1])==0 ||
                side=='R' && sign(A.D[i])+sign(A.D[i-1])==0
                # println("Remedy 1 ")
                # Compute the inverse of the original arrowhead (DPR1)
                Ainv,Krho,Qout1 = inv(A,0.0,tols[4]) # Ainv is Float64
                Qout==1 && (Qout=Qout+4)
                if isinf(Ainv.r) || isnan(Ainv.r)
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
    Float64(lambda), Qout
    
end # eig (k)

function rootsWDK{T,T1}(p::Poly{T},x0::Vector{T1},steps)
    # Polynomial roots with Weierstrass, Durand, Kerner method
    n=length(x0)
    x=x0
    for k=1:steps        
        for i=1:n
            x[i]=x[i]-polyval(p,x[i])/(p[end]*prod(x[i]-x[collect([1:i-1;i+1:n])]))
        end
    end
    x
end

