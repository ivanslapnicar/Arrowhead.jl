#-------- Inverses of DPR1

function invA{T}(A::SymDPR1{T},i::Integer,tols::Vector{Float64})

# COMPUTES: inverse of a shifted SymDPR1 matrix A=diagm(A.D)+A.r*A.u*A.u',
# invA(A-A.D[i]*I) which is a SymArrow.
# Uses higher precision to compute top of the arrow element accurately, if
# needed. 
# tols=[tolb,tolz] are tolerances, usually [1e3, 10*n]
# [0.0,0.0] forces DoubleDouble, [1e50,1e50] would never use it
# RETURNS:  SymArrow(D,z,b,i), Kb, Kz, Qout
# Kb - condition Kb, Kz - condition Kz, Qout = 1 / 0 - double was / was not used 

n=length(A.D)
D=Array(T,n-1)
z=Array(T,n-1)
wz=one(T)/A.u[i]
shift=A.D[i]

for k=1:i-1
    D[k]=one(T)/(A.D[k]-shift)
    z[k]=-A.u[k]*D[k]*wz
end
for k=i+1:n
    D[k-1]=one(T)/(A.D[k]-shift)
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

Kz=(sumabs(A.u)-abs(A.u[i]))*abs(wz)


# Kz=maxabs(A.z)*abs(wz)

if Kb<tols[1] ||  Kz<tols[2]
   b=(P+Q)*wz*wz
else  # recompute in Double
   Qout=1
   shiftd=map(Double,A.D[i])
   Dd=[Double{Float64}[Double(A.D[k])-shiftd for k=1:i-1], 
      Double{Float64}[Double(A.D[k])-shiftd for
   k=i+1:length(A.D)]]
   wzd=Double(A.u[i])

   Pd,Qd=map(Double,(0.0,0.0))
   
   for k=1:i-1
      Dd[k].hi>0.0 ? Pd+=Double(A.u[k])^2/Dd[k] :
      Qd+=Double(A.u[k])^2/Dd[k]
   end
   
   for k=i+1:n
      Dd[k-1].hi>0.0 ? Pd+=Double(A.u[k])^2/Dd[k-1] :
      Qd+=Double(A.u[k])^2/Dd[k-1]
      # @show P,Q
   end 

   A.r > 0 ?   Pd+=Double(1.0)/Double(A.r)  : Qd+=Double(1.0)/Double(A.r)

   bd=(Pd+Qd)/(wzd*wzd)
   b=bd.hi+bd.lo
end

# return this
SymArrow(D,z,b,i),Kb,Kz,Qout

end # invA


#=

function invA{T}(A::SymArrow{T}, shift::Float64, tolr::Float64)

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
   r=Double(1.0)/(Pd+Qd)
   rho=-(r.hi+r.lo)
end

# returns the following
SymDPR1(D,u,rho), Krho, Qout

end # invA

function invA{T}(A::SymArrow{T}, shift::Double)

# COMPUTES: inverse of the shifted SymArrow A, inv(A-shift*I), which is a SymDPR1
# here shift is Double so it uses Double to compute everything 
# RETURNS: SymDPR1(D1,u1,rho), Qout
# Qout = 1 / 0 - Double was / was not used 

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

end # invA


function  aheig( A::SymArrow,k::Integer,tols::Vector{Float64})

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
v=zeros(n)
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
   Fmiddle = (atemp-middle)-sum(A.z.^2./(Dtemp-middle))
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

Ainv,Kb,Kz,Qout = invA(A,i,tols[1:2])

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
    for k=1:n-1
       nu1=max(nu1,abs(Ainv.D[k])+abs(Ainv.z[k]))
    end
    nu1=max(sumabs(Ainv.z)+abs(Ainv.a), nu1)

    Knu=nu1/abs(nu)
    if Knu<tols[3]
        # Accuracy is fine, compute the eigenvector
        mu = 1.0/nu
        v=[ A.z./((A.D-sigma)-mu);-1.0]
#        for k=1:n-1
#	   v[k] = A.z[k]/((A.D[k]-sigma)-mu)
#        end 
#        v[n]=-one
        lambda, v = mu+sigma, v/norm(v)
    else
        # Remedies according to Remark 3 - we shift between original
        # eigenvalues and compute DPR1 matrix
        # 1/nu1+sigma, 1/nu+sigma
        println("Remedy 3 ")
        nu = side=='R' ? abs(nu) : -abs(nu)
        nu1=-sign(nu)*nu1
        sigma1=(nu1+nu)/(2.0*nu*nu1)+sigma

        if findfirst(A.D-sigma1,0.0)>0 # we came back with a pole
            # recompute sigmav more accurately according with dekker
            sigmav=(Double(nu1)+Double(nu))/(Double(2.0)*Double(nu)*Double(nu1))+Double(sigma)
            # Compute the inverse of the shifted arrowhead (DPR1)
            Ainv,Qout1=invA(A,sigmav) # Ainv is Float64
            nu1=bisect(Ainv,side) 
            mu1 = 1.0/nu1
            D0=map(A.D,Double)-sigmav
            D1=D0.hi+D0.lo
            if findfirst(D1-mu1,0.0)>0 
                ind=find(D1-mu1==0.0);
                v=zeros(n)
                v[ind]=1.0;
            else
                v=[ A.z./(D1-mu1);-1.0]
            end
            # Shift the eigenvalue back in Double
            lam = Double(1.0)/Double(nu1)+sigmav
            # Return this
            lambda,v = lam.hi+lam.lo, v/norm(v)       
        else
            # Compute the inverse of the shifted arrowhead (DPR1)
            Ainv, Krho,Qout1=invA(A,sigma1,tols[4]) # Ainv is Float64
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

        if k==1 & A.D[1]<0.0 | k==n & A.D[n-1]>0.0 | side=='L' & sign(A.D[i])+sign(A.D[i+1])==0 | side=='R' & sign(A.D[i])+sign(A.D[i-1])==0
            println("Remedy 1 ")
            # Compute the inverse of the original arrowhead (DPR1)
            Ainv,Krho,Qout1 = invA(A,0.0,tols[4]) # Ainv is Float64
            Qout==1 && (Qout=Qout+4)
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
lambda,v,i,Kb,Kz,Knu,Krho,Qout

end # aheig


function aheigall(A::SymArrow, tols::Vector{Float64})

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
U=eye(n,n)
E=zeros(n)
Kb=zeros(n); Kz=zeros(n); Knu=zeros(n); Krho=zeros(n)
Qout=zeros(Int,n); Sind=zeros(Int,n) 

# Quick return for 1x1
if n==1
    return U,A.a,Sind,Kb,Kz,Knu,Krho,Qout
end

#  test for deflation in z
z0=find(z.==0)
zx=find(z.!=0)
if !isempty(z0)
   E[z0]=D[z0]
   D=D[zx]
   z=z[zx]
   if isempty(z)
      E[n]=A.a
      return U,E,Sind,Kb,Kz,Knu,Krho,Qout
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
   for k=1:nn+1
      E[zxx[k]],U[zxx,zxx[k]],Sind[zxx[k]],Kb[zxx[k]],Kz[zxx[k]],Knu[zxx[k]],Krho[zxx[k]],Qout[zxx[k]]=aheig(SymArrow(D[gx],z[gx],A.a,nn+1),k,tols)
   end

   for l=1:lg0
      # multiplication by eye is a temporary fix by Andreas
      U=(R[l]*eye(n0))'*U
   end 

else

   # No deflation in D
   for k=1:n
      E[zx[k]],U[zx,zx[k]],Sind[zx[k]],Kb[zx[k]],Kz[zx[k]],Knu[zx[k]],Krho[zx[k]],Qout[zx[k]]=aheig(SymArrow(D,z,A.a,n),k,tols)
   end
end

# back premutation of vectors
isi=sortperm(is)

# must sort E once more 
es=sortperm(E,rev=true)
E=E[es]

U=U[[isi[1:A.i-1],n0,isi[A.i:n0-1]],es]

# Return this
U,E,Sind,Kb,Kz,Knu,Krho,Qout

end # aheigall

=#

