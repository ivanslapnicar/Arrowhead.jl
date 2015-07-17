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
    HalfArrow(rand(n-1)-0.5,rand(n-1)-0.5) :
    # This is one with the point
    HalfArrow(rand(n-1)-0.5,rand(n)-0.5)
end


#-------- Inverses  

function inv{T}(A::HalfArrow{T},i::Integer,tols::Vector{Float64})

    # COMPUTES: inverse of a SymArrow matrix B, where B=A'*A, and A =[A.D, A.z] 
    # is the HalfArrow. Here inv(B-B.D[i]*I) is again SymArrow.
    # Uses higher precision to compute top of the arrow element accurately, if
    # needed. 
    # tols=[tolb,tolz] are tolerances, usually [1e3, 10*n]
    # [0.0,0.0] forces DoubleDouble, [1e50,1e50] would never use it
    # RETURNS:  SymArrow(D,z,b,i), Kb, Kz, Qout
    # Kb - condition Kb, Kz - condition Kz, Qout = 1 / 0 - double was / was not used 
    
    n=length(A.D)
    D=Array(T,n)
    z=Array(T,n)
    z1=Array(T,n)
    z1=A.z[1:n].*A.D
    
    wz=1/(z1[i])
    shift=A.D[i]
    
    for k=1:i-1
        D[k]=one(T)/((A.D[k]-shift)*(A.D[k]+shift))
        z[k]=-z1[k]*D[k]*wz
    end
    for k=i+1:n
        D[k-1]=one(T)/((A.D[k]-shift)*(A.D[k]+shift))
        z[k-1]=-z1[k]*D[k-1]*wz
    end
    
    # D=[1./(A.D[1:i-1]-shift),1./(A.D[i+1:end]-shift)]
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
    a<0 ? P=P-a : Q=Q-a
    
    Kb=(P-Q)/abs(P+Q)
    Kz=maxabs(z1)*abs(wz)
    
    if Kb<tols[1] ||  Kz<tols[2]
        b=(P+Q)*wz*wz
    else  # recompute in Double
        Qout=1
        shiftd=map(Double,A.D[i])
        Dd=[Double{Float64}[(Double(A.D[k])+shiftd)*(Double(A.D[k])-shiftd) for k=1:i-1], 
            Double{Float64}[(Double(A.D[k])+shiftd)*(Double(A.D[k])-shiftd) for
                            k=i+1:length(A.D)]]
        wzd=Double(A.z[i])*shiftd
        
    #    ad=doubledot(A.z,A.z)-shiftd^2 # we already have this if we use doubledot above
        
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
        
        ad.hi<0 ?   Pd=Pd-ad : Qd=Qd-ad 
        
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



function inv{T}(A::HalfArrow{T}, shift::Float64, tolr::Float64)
    
    # COMPUTES: inverse of a SymArrow matrix B, where B=A'*A, and A =[A.D, A.z] 
    # is the HalfArrow. Here inv(B-shift^2*I) is a SymDPR1.
    # uses DoubleDouble to compute rho accurately, if needed. 
    # tolr is tolerance, usually 1e3,  0.0 forces Double, 1e50 would never use it
    # RETURNS: SymDPR1(D,u,rho), Krho, Qout
    # Krho - condition Krho, Qout = 1 / 0 - Double was / was not used 
    
    n=length(A.D)
    D=Array(T,n+1)
    u=Array(T,n+1)
    z1=Array(T,n)
    z1=A.z[1:n].*A.D
    
    # Ds=A.D-shift
    
    # D=[1./Ds[1:A.i-1];0.0;1./Ds[A.i:end]];
    # u=[A.z[1:A.i-1]./Ds[1:A.i-1];-1.0;A.z[A.i:end]./Ds[A.i:end]];
    
    for k=1:n
        D[k]=one(T)/((A.D[k]-shift)*(A.D[k]+shift))
        u[k]=z1[k]*D[k]
    end
    D[n+1]=zero(T)
    u[n+1]=-one(T)
    
    
    #=
    for k=A.i:n
        D[k+1]=one(T)/((A.D[k]-shift)*(A.D[k]-shift))
        u[k+1]=z1[k]*D[k+1]
    end
    =#
    
    # compute rho and Krho
    #--- compute the sum in a plain loop
    
    P=zero(T)
    Q=zero(T)
    Qout=0

    # a=dot(A.z,A.z)-shift^2 # try with more accuracy
    ad=doubledot([A.z;shift],[A.z;-shift])
    a=ad.hi+ad.lo   

    for k=1:n
        D[k]>0.0 ? P+=z1[k]^2*D[k] : Q+=z1[k]^2*D[k]
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
        Dd=Double{Float64}[(Double(A.D[k])+shiftd)*(Double(A.D[k])-shiftd) for k=1:n]
        # ad=doubledot(A.z,A.z)-shiftd^2 # we already have this if doubledot is used above
        #   Dd=[Double(A.D[k])-shiftd for k=1:length(A.D) ] 
        #   Dd=map(Double,A.D)-shiftd
        
        for k=1:n
            Dd[k].hi>0.0 ? Pd+=(Double(A.z[k])*Double(A.D[k]))^2/Dd[k] :
            Qd+=(Double(A.z[k])*Double(A.D[k]))^2/Dd[k]
        end
        
        ad.hi+ad.lo<0 ? Pd=Pd-ad : Qd=Qd-ad
        r=Double(1.0)/(Pd+Qd)
        rho=-(r.hi+r.lo)
    end
    
    # returns the following
    SymDPR1(D,u,rho), Krho, Qout
    
end # inv


function inv{T}(A::HalfArrow{T}, shift::Double)
    
    # COMPUTES: inverse of a SymArrow matrix B, where B=A'*A, and A =[A.D, A.z] 
    # is the HalfArrow. Here inv(B-shift^2*I) is SymDPR1.
    # Here shift is Double so it uses Double to compute everything 
    # RETURNS: SymDPR1(D1,u1,rho), Qout
    # Qout = 1 on exit meaning Double was used 
    
    n=length(A.D)
    D=Array(Double,n+1)
    u=Array(Double,n+1)
    
    oned=Double(1.0)
    zerod=Double(0.0)
    
    # z1=Array(T,n)
    # z1=A.z[1:n].*A.D
    
    # D[1:n]=map(Double,A.D)-shift
    for k=1:n
        D[k]=(Double(A.D[k])-shift)*(Double(A.D[k])+shift)
        u[k]=Double(A.z[k])*Double(A.D[k])/D[k] 
    end
    
    a=doubledot(A.z,A.z)-shift^2
    
    D[n+1]=zerod
    u[n+1]=Double(-1.0)
    
    # compute rho and Krho
    #--- compute the sum in a plain loop
    
    P,Q=zerod,zerod
    Qout=1
    
    for k=1:n
        D[k].hi>0.0 ? P+=(Double(A.z[k])*Double(A.D[k]))^2/D[k] :
        Q+=(Double(A.z[k])*Double(A.D[k]))^2/D[k]
        D[k]=oned/D[k]
    end
    
    a.hi+a.lo<0.0  ? P=P-a : Q=Q-a
    r=oned/(P+Q)
    
    # for k=1:length(A.D)
    #    D1[k].hi+D1[k].lo>0.0 ? P=P+z1[k]^2/D1[k] : Q=Q+z1[k]^2/D1[k]
    # end
    
    rho=-(r.hi+r.lo)
    
    # returns the following
    D1=Array(T,n+1)
    u1=Array(T,n+1)
    for k=1:n+1
        D1[k]=D[k].hi+D[k].lo
    end
    for k=1:n+1
        u1[k]=u[k].hi+u[k].lo
    end
    
    # SymDPR1(T[x.hi+x.lo for x=D],T[x.hi+x.lo for x=u],rho), Qout
    
    SymDPR1(D1,u1,rho), Qout
    
end # inv


function  svd( A::HalfArrow,k::Integer,tols::Vector{Float64})
    
    # COMPUTES: k-th singular value triple of an ordered irreducible HalfArrow
    # A = [diagm(A.D) A.z] with A.D > 0
    # tols=[tolb,tolz,tolnu,tolrho,tollambda] = [1e3,10.0*n,1e3,1e3,1e3]
    # RETURNS: lambda, u, v, Sind, Kb, Kz, Knu, Krho, Qout
    # lambda - k-th singular value in descending order
    # u - lambda's normalized left singular vector
    # v - lambda's normalized right singular vector
    # Kb, Kz, Knu, Krho - condition numbers
    # Qout = 1 / 0 - Double was / was not used 
    
    # If Qout>0, quad precision was used
    # i was the shift index
    
    # Set the dimension
    n = length(A.D) + 1
    
    # Set all conditions initially to zero
    Kb,Kz,Knu,Krho=0.0,0.0,0.0,0.0
    Qout=0
    u=zeros(length(A.z))
    v=zeros(n)
    z1=Array(Float64,n)
    z1=A.z[1:n-1].*A.D
    
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
        Dtemp = (A.D-A.D[k]).*(A.D+A.D[k])
        atemp = dot(A.z,A.z)-A.D[k]^2
        middle = Dtemp[k-1]/2.0
        Fmiddle = (atemp-middle)-sum(z1.^2./(Dtemp-middle))
        # middle=(A.D[k-1]-A.D[k])/2.0
        # Fmiddle=0.0
        # for l=1:n-1
        #     Fmiddle+=A.z[l]^2/((A.D[l]-A.D[k])-middle)
        # end
        # Fmiddle=((A.a-A.D[k])-middle)-Fmiddle
        # @show Fmiddle
        sigma,i,side = Fmiddle < 0.0 ? (A.D[k],k,'R') : (A.D[k-1],k-1,'L')
    end
    
    # Compute the inverse of the shifted matrix, (A'*A)_i^(-1), Kb and Kz
    
    Ainv,Kb,Kz,Qout = inv(A,i,tols[1:2])
    
    # Compute the eigenvalue of the inverse shifted matrix
    
    nu = bisect( Ainv,side )
    
    #  nu=fastzero([invD1; 0; invD2], [w1;wz;w2], b, side); # implement later
    #  [nu-nu1]/nu, [nu-nueva]/nu, pause  # and compare
    
    if abs(nu)==Inf
        # this is nonstandard
        # Deflation in aheig (nu=Inf)
        u[i]=1.0
        v[i]=1.0
        lambda=sigma
    else
        # standard case, full computation
        # nu1 is the F- or 1-norm of the inverse of the shifted matrix
        # nu10=maximum([sum(abs(Ainv.z))+abs(Ainv.a), maximum(abs(Ainv.D)+abs(Ainv.z))])
        nu1=0.0
        for k=1:n-1
            nu1=max(nu1,abs(Ainv.D[k])+abs(Ainv.z[k]))
        end
        nu1=max(sumabs(Ainv.z)+abs(Ainv.a), nu1)
        
        Knu=nu1/abs(nu)
        if Knu<tols[3]
            # Accuracy is fine, compute the eigenvector
            mu = 1.0/nu
            #         v=[ z1./((A.D-sigma).*(A.D+sigma)-mu);-1.0]
            for k=1:n-1
	        v[k] = z1[k]/((A.D[k]-sigma)*(A.D[k]+sigma)-mu)
            end 
            v[n]=-1.0
            v=v/norm(v)
            lambda=sqrt(mu+sigma^2) # this may have errors
            u[1:n-1]=lambda*v[1:n-1]./A.D
            if length(A.z)==n 
                u[n]=A.z[n]*v[n]/lambda
            end
        else
            # Remedies according to Remark 3 - we shift between original
            # eigenvalues and compute DPR1 matrix
            # 1/nu1+sigma, 1/nu+sigma
            println("Remedy 3 ")
            nu = side=='R' ? abs(nu) : -abs(nu)
            nu1=-sign(nu)*nu1
            sigma1=(nu1+nu)/(2.0*nu*nu1)+sigma^2
            
            if findfirst(abs(A.D)-sqrt(sigma1),0.0)>0 # we came back with a pole
                # recompute sigmav more accurately according with dekker
                sigmav=(Double(nu1)+Double(nu))/(Double(2.0)*Double(nu)*Double(nu1))+Double(sigma)^2
                # Compute the inverse of the shifted arrowhead (DPR1)
                Ainv,Qout1=inv(A,sqrt(sigmav)) # Ainv is Float64, here it was sqrt(sigmav)
                
                nu1=bisect(Ainv,side) 
                mu1 = 1.0/nu1
                D0=map(A.D,Double)
                D0=A.D.*A.D-sigmav
                D1=D0.hi+D0.lo
                
                # Shift the eigenvalue back in Double
                lam = Double(1.0)/Double(nu1)+sigmav
                lambda=sqrt(lam.hi+lam.lo)     
                
                if findfirst(D1-mu1,0.0)>0 
                    ind=find(D1-mu1==0.0)
                    u[ind]=1.0
                    v[ind]=1.0
                else
                    for k=1:n-1
	                v[k] = z1[k]/(D1[k]-mu1)
                    end 
                    v[n]=-1.0
                    v=v/norm(v)
                    u[1:n-1]=lambda*v[1:n-1]./A.D
                    if length(A.z)==n 
                        u[n]=A.z[n]*v[n]/lambda
                    end
                    # v=[ A.z./(D1-mu1);-1.0]
                end
            else
                # Compute the inverse of the shifted arrowhead (DPR1)
                ssigma=sqrt(sigma1)
                Ainv, Krho,Qout1=inv(A,ssigma,tols[4]) # Ainv is Float64
                # Compute the eigenvalue by bisect for DPR1
	        # Note: instead of bisection could use dlaed4 (need a wrapper) but
	        # it is not faster. There norm(u)==1
                nu1= bisect(Ainv,side)
                mu1=1.0/nu1 
                
                for k=1:n-1
	            v[k] = z1[k]/((A.D[k]-ssigma)*(A.D[k]+ssigma)-mu1)
                end 
                v[n]=-1.0
                v=v/norm(v)
                lambda=sqrt(mu1+sigma1)
                u[1:n-1]=lambda*v[1:n-1]./A.D
                if length(A.z)==n 
                    u[n]=A.z[n]*v[n]/lambda
                end
                
                # v=[ A.z./((A.D-sigma1)-mu1);-1.0]
                # Return this - shift the eigenvalue back and normalize the vector
            end
            Qout==1 && (Qout=Qout+2)
        end
        
        # Remedy according to Remark 1 - we recompute the the eigenvalue
        # near zero from the inverse of the original matrix (a DPR1 matrix).  
        # Here this can only happen for k==n

        if k==n && (abs(A.D[i])^2+abs(1.0/nu))/abs(lambda^2)>tols[5]
            
            # if k==1 && A.D[1]<0.0 || k==n && A.D[n-1]>0.0 || side=='L' && sign(A.D[i])+sign(A.D[i+1])==0 || side=='R' && sign(A.D[i])+sign(A.D[i-1])==0
            println("Remedy 1 ")
            # Compute the inverse of the original arrowhead (DPR1)
            Ainv,Krho,Qout1 = inv(A,0.0,tols[4]) # Ainv is Float64
            Qout==1 && (Qout=Qout+4)
            if abs(Ainv.r)==Inf
                lambda=0.0
            else                
                # Here we do not need bisection. We compute the Rayleigh
                # quotient by using already computed vectors which are
                # componentwise accurate
                nu1=sum(v.^2.*Ainv.D)+Ainv.r*sum(v.*Ainv.u)^2;
                lambda=sqrt(1.0/nu1)
            end
        end
    end
    
    # Return this
    lambda,u,v,i,Kb,Kz,Knu,Krho,Qout
    
end # svd (k)


function svd(A::HalfArrow, tols::Vector{Float64})

    # COMPUTES: all singular values and singular vectors of a real HalfArrow
    # A = [diagm(A.D) A.z]
    # tols=[tolb,tolz,tolnu,tolrho,tollambda] = [1e3,10.0*n,1e3,1e3,1e3]
    # RETURNS: U, S, V, Sind, Kb, Kz, Knu, Krho, Qout
    # U = left singular vectors
    # E = singular values in decreasing order
    # V = right singular vectors
    # Sind[k] - shift index i for the k-th eigenvalue
    # Kb, Kz, Knu, Krho [k] - respective conditions for the k-th eigenvalue
    # Qout[k] = 1 / 0 - Double was / was not used when computing k-th eigenvalue 
    
    nd=length(A.D)
    n=length(A.z)
    n0=n
    
    # Ordering the matrix
    signD=sign(A.D)
    D=abs(A.D)
    is=sortperm(D,rev=true)
    D=D[is]
    signD=signD[is]
    z=A.z[1:n] 
    z=z[is]
    if n==nd
        U=eye(n,n)
        V=eye(n+1,n)
    else
        U=eye(n,n)
        V=eye(n,n)
    end
    E=zeros(n)

    Kb=zeros(n); Kz=zeros(n); Knu=zeros(n); Krho=zeros(n)
    Qout=zeros(Int,n); Sind=zeros(Int,n) 
    
    # Quick return for 1x1 - this does not exist for HalfArrow
    # if n==1
    #     return U,A.a,V,Sind,Kb,Kz,Knu,Krho,Qout
    # end
    
    #  test for deflation in z
    z0=find(z.==0)
    zx=find(z.!=0)

    if isempty(zx)  # nothing to do
        if n==nd
            E=[A.D]
            isE=sortperm(E,rev=true)
            E=E[isE]
            V[:,1:n]=V[:,isE]
            U=U[invperm(isE),:]    
        else
            E=[A.D;z[n]]
            isE=sortperm(E,rev=true)
            E=E[isE]
            V=V[:,isE]
            U=U[invperm(isE),:]
        end
        return U,E,V,Sind,Kb,Kz,Knu,Krho,Qout
    end
    
    if !isempty(z0)
        E[z0]=D[z0]
        for k=1:length(z0)
            V[z0[k],z0[k]]=signD[z0[k]]
        end
        D=D[zx]
        # signD=signD[zx]
        z=z[zx]
        if isempty(z)
            if n>nd
                E[n]=abs(A.z[n])
                V[n,n]=sign(A.z[n])
            end
            # return U,E,V,Sind,Kb,Kz,Knu,Krho,Qout
        end
    else
        if n0==nd
            zxv=[zx,n0+1]
        else
            zx=zxv=[zx,n0]
            z=[z;A.z[n0]]
        end  
        n=length(z)
        
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
            R=Array(Givens{Float64},lg0)
            for l=lg0:-1:1
                R[l]=givens(z[g0[l]],z[g0[l]+1],zx[g0[l]],zx[g0[l]+1],n0)
                z[g0[l]]=R[l].r; z[g0[l]+1]=0.0
                # A_mul_Bc!(U,R) # multiply R'*U later
                E[zx[g0[l]+1]]=D[g0[l]+1]
            end    
            # remains
            gx=[0;gx]+1
            nn=length(gx)
            
            zxx=zx[[gx;n]]
     
            if n0==nd
                zxxv=[zxx,n0+1]
            else
                zxx=zxxv=[zxx,n0]
                # z=[z;A.z[n0]]
            end  

            for k=1:nn+1
                E[zxx[k]],U[zxx,zxx[k]],V[zxxv,zxx[k]],Sind[zxx[k]],Kb[zxx[k]],Kz[zxx[k]],Knu[zxx[k]],Krho[zxx[k]],
                Qout[zxx[k]]=svd(HalfArrow(D[gx],z[gx]),k,tols)
            end
            
            for l=1:lg0
                # manual transpose
                R1=Base.LinAlg.Givens(R[l].size, R[l].i1,R[l].i2,R[l].c,-R[l].s,R[l].r)
                U=R1*U
            end 
            
        else
            
            # No deflation in D
            # @show n, zx, zxv, U, V,D,z
            for k=1:n
                # V[zxv,zx[k]]
                E[zx[k]],U[zx,zx[k]],V[zxv,zx[k]],Sind[zx[k]],Kb[zx[k]],Kz[zx[k]],Knu[zx[k]],Krho[zx[k]],Qout[zx[k]]=
                svd(HalfArrow(D,z),k,tols)
                # if signD[zx[k]]<0.0
                #     V[zxv,zx[k]]=-V[zxv,zx[k]]
                # end
                # @show V[zxv,zx[k]],signD, zxv, zxv[1:end-1]
                V[zxv[1:end-1],zx[k]]=V[zxv[1:end-1],zx[k]].*signD[zxv[1:end-1]]
                # @show k,U,E,V
            end
        end
    end
    # back premutation of vectors
    isi=sortperm(is)
    
    # must sort E once more 
    es=sortperm(E,rev=true)
    E=E[es]
    if n0==nd  
        U=U[isi,es]
        V=V[[isi;n0+1],es] 
    else
        isiv=[isi;n0]
        U=U[isiv,es]
        V=V[isiv,es]
    end
    
    # Return this
    U,E,V,Sind[es],Kb[es],Kz[es],Knu[es],Krho[es],Qout[es]
    
end # svd (all)


# this is just for convenience
tols=[1e3,1000.0,1e3,1e3,1e3]
            
