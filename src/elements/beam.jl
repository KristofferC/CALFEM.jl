using LinearAlgebra

"""
beam2e(ex,ey,elem_prop) -> Ke
beam2e(ex,ey,elem_prop,eq) -> Ke,fe

Compute the stiffness matrix for a two dimensional beam element.
The optional eq vector is an eventual uniformly distributed load over the element.
"""

function beam2e(ex::Union{LinearAlgebra.Transpose, Adjoint,Vector}, ey::Union{LinearAlgebra.Transpose, Adjoint,Vector}, elem_prop::Union{LinearAlgebra.Transpose, Adjoint,Vector},eq::Union{LinearAlgebra.Transpose, Adjoint,Vector,Nothing}=nothing)


    b=[[ex[2]-ex[1]],[ey[2]-ey[1]]]
    L = sqrt(b'*b)
    n = reshape(b'/L,2,)

    E=elem_prop[1]
    A=elem_prop[2]
    I=elem_prop[3]

    qx=0.0
    qy=0.0
    if  eq != nothing
        qx=eq[1]
        qy=eq[2]
    end
    Kle = [E*A/L       0.0           0.0     -E*A/L     0.0         0.0
           0.0     12*E*I/L^3.0  6*E*I/L^2.0     0.0  -12*E*I/L^3.0  6*E*I/L^2.0
           0.0     6*E*I/L^2.0   4*E*I/L       0.0  -6*E*I/L^2.0   2*E*I/L
           -E*A/L      0.0           0.0      E*A/L     0.0         0.0
           0.0    -12*E*I/L^3.0 -6*E*I/L^2.0     0.0   12*E*I/L^3.0 -6*E*I/L^2.0
           0.0     6*E*I/L^2.0   2*E*I/L       0.0   -6*E*I/L^2.0  4*E*I/L]'


    fle=L*[qx/2, qy/2, qy*L/12, qx/2, qy/2, -qy*L/12]

    G=[ n[1]  n[2]   0.0     0.0     0.0    0.0
        -n[2]  n[1]   0.0     0.0     0.0    0.0
        0.0     0.0     1.0     0.0     0.0    0.0
        0.0     0.0     0.0    n[1]   n[2]   0.0
        0.0     0.0     0.0   -n[2]   n[1]   0.0
        0.0     0.0     0.0     0.0     0.0    1.0]'


    Ke=G'*Kle*G
    fe=G'*fle

    if eq == nothing
        return Ke
    else
        return Ke,fe
    end

end






"""
beam2s(ex,ey,elem_prop,eq,nelem_prop) -> es,edi,eci

Compute section forces in two dimensional beam element (beam2e).

eq and nelem_prop are optional arguments, specifying an element UDL
and number of evaluation points of the element respectively


Returned values:

es = [ N V M ]    section forces, local directions

edi = [ u1 v1 ]     element displacements, local directions
eci = [ x1 ]    local x-coordinates of the evaluation points on element

"""

function beam2s(ex::Union{LinearAlgebra.Transpose, Adjoint,Vector}, ey::Union{LinearAlgebra.Transpose, Adjoint,Vector}, elem_prop::Union{LinearAlgebra.Transpose, Adjoint,Vector},el_disp::Union{LinearAlgebra.Transpose, Adjoint,Vector},eq::Union{LinearAlgebra.Transpose, Adjoint,Vector,Nothing}=nothing,nelem_prop::Union{Number,Nothing}=nothing)

    EA=elem_prop[1]*elem_prop[2]
    EI=elem_prop[1]*elem_prop[3]

    dx = ex[2]-ex[1]
    dy = ey[2]-ey[1]
    L = sqrt(dx^2 + dy^2)

    b=[dx dy]
    n = b/L#Kolla reshape

    qx=0.0
    qy=0.0

    if eq !=nothing
        qx=eq[1]
        qy=eq[2]
    end

    ne=2

    if nelem_prop!=nothing
        ne = nelem_prop
    end


    C=[0   0   0    1   0   0
     0   0   0    0   0   1
     0   0   0    0   1   0
     L   0   0    1   0   0
     0   L^3  L^2 0   L   1
     0 3*L^2 2*L  0   1   0]

    #n=b/L

    G=[n[1] n[2]  0    0    0   0
      -n[2] n[1]  0    0    0   0
        0    0    1    0    0   0
        0    0    0   n[1] n[2] 0
        0    0    0  -n[2] n[1] 0
        0    0    0    0    0   1]'

    M=inv(C)*(G*el_disp-[0 0 0 -qx*L^2/(2*EA) qy*L^4/(24*EI) qy*L^3/(6*EI)]' )

    A=[M[1] M[4]]';  B=[M[2] M[3] M[5] M[6]]';

    x=collect(0:L/(ne-1):L);   zero=zeros(length(x));    one=ones(length(x));
    u=[x one]*A-(x.^2)*qx/(2*EA);
    du=[one zero]*A-x*qx/EA;
    v=[x.^3 x.^2 x one]*B+(x.^4)*qy/(24*EI);
    d2v=[6*x 2*one zero zero]*B+(x.^2)*qy/(2*EI);
    d3v=[6*one zero zero zero]*B+x*qy/EI;

    N=EA*du
    M=EI*d2v
    V=-EI*d3v
    edi=hcat(u,v)
    eci=x
    es=hcat(N,V,M)

    return (es,edi,eci)

end
