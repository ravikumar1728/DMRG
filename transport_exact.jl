using LinearAlgebra

Sx=[0 0.5; 0.5 0]
Sy=[0 -0.5im; 0.5im 0]
Sz=[0.5 0; 0 -0.5]
Id=I(2)

SxSx=kron(Sx,Sx)
SySy=kron(Sy,Sy)
SzSz=kron(Sz,Sz)

function op_twosite(op,i,L)
    if i==1
        s=kron(op,I(2^(L-2)))
    elseif i==L-1
        s=kron(I(2^(L-2)),op)
    elseif i > 1 && i < L-1
        s=kron(I(2^(i)),kron(op,I(2^(L-i-2))))
    end
    return s
end

function op_onesite(op,i,L)
    if i==1
        s=kron(op,I(2^(L-1)))
    elseif i==L
        s=kron(I(2^(L-1)),op)
    elseif i > 1 && i < L
        s=kron(I(2^(i)),kron(op,I(2^(L-i-1))))
    end
    return s
end


L=10
gamma=pi/4
delta=[(i<=(L/2)) ? 0 : cos(gamma)  for i in 1:L]
h_i=zeros(2^L,2^L)
hlocal=[]
for i in 1:(L-1)
    h_i=h_i+op_twosite(SxSx,i,L)+op_twosite(SySy,i,L)+delta[i]*op_twosite(SzSz,i,L)
    h=SxSx+SySy+delta[i]*SzSz
    push!(hlocal,h)
end 


t=100

measurement1=[]
measurement3=[]

for j in 1:t
    U_t=exp(-j*0.1im*h_i)
    psi_t=U_t*psi[:,1]
    
    measurement2=[]
    measurement4=[]


    for i in 1:L
        rho=psi_t[:,1]*reshape(conj(psi_t[:,1]),(1,1024))
        mea=tr(rho*op_onesite(Sz,i,L))
        push!(measurement2,real(mea))
    end
    push!(measurement1,measurement2)
    
    for i in 1:L-1
        rho=psi_t[:,1]*reshape(conj(psi_t[:,1]),(1,1024))
        mea=tr(rho*op_twosite(hlocal[i],i,L))
        push!(measurement4,real(mea))
    end
    push!(measurement3,measurement4)
        
end
        