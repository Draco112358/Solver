using MKL

function ComputeMatrixVectorNew(x, w, incidence_selection,
    P_data, Lp_data, Z_self, Yle, invZ, invP, F, resProd)

    m = size(incidence_selection["A"], 1)
    ns = size(incidence_selection["Gamma"], 2)

    # Extract vectors from x
    I = x[1:m]
    Q = x[m+1:m+ns]
    Phi = x[m+ns+1:end]

    mx = incidence_selection["mx"]
    my = incidence_selection["my"]
    mz = incidence_selection["mz"]

    # Initialize Y1
    Y1 = zeros(ComplexF64, m)
    Y1[1:mx] .= Lp_data[:Lp_x] * I[1:mx]
    Y1[mx+1:mx+my] .= Lp_data[:Lp_y] * I[mx+1:mx+my]
    Y1[mx+my+1:mx+my+mz] .= Lp_data[:Lp_z] * I[mx+my+1:mx+my+mz]

    # Update Y1 with additional terms
    Y1 .= 1im * w * Y1 .+ Z_self .* I .+ incidence_selection["A"] * Phi

    # Compute Y2
    Y2 = P_data[:P] * Q .- transpose(incidence_selection["Gamma"]) * Phi

    # Compute Y3
    Y3 = -transpose(incidence_selection["A"]) * I .+ Yle * Phi .+ 1im * w * (incidence_selection["Gamma"] * Q)

    # Combine results
    MatrixVector = precond_3_3_vector(
        F, invZ, invP, incidence_selection["A"], incidence_selection["Gamma"], w, Y1, Y2, Y3
    )

    return MatrixVector
end


function precond_3_3_vector(lu,invZ,invP,A,Gamma,w,X1,X2,X3)

    n1=length(X1)
    n2=length(X2)
    n3=length(X3)

    i1=range(1, stop=n1)
    i2=range(n1+1,stop=n1+n2)
    i3=range(n1+n2+1,stop=n1+n2+n3)

    Y=zeros(ComplexF64 , n1+n2+n3)

    M1 = prod_real_complex(invZ, X1)
    M2 = lu\(prod_real_complex(transpose(A), M1))
    M3 = prod_real_complex(invP, X2)
    M4 = lu\prod_real_complex(Gamma, M3)
    M5 = lu\X3

    Y[i1] .= Y[i1] .+ M1-1.0*(prod_real_complex((invZ),prod_real_complex((A), M2)))
    Y[i1] .= Y[i1] .+ 1im*w*(prod_real_complex((invZ),prod_real_complex((A), M4)))
    Y[i1] .= Y[i1] .- 1.0*(prod_real_complex((invZ),prod_real_complex((A), M5)))

    Y[i2] .= Y[i2] .+ (prod_real_complex(invP,prod_real_complex((transpose(Gamma)), M2)))
    Y[i2] .= Y[i2] .+ M3 - 1im*w*(prod_real_complex(invP,prod_real_complex(transpose(Gamma), M4)))
    Y[i2] .= Y[i2] .+ (prod_real_complex(invP,prod_real_complex(transpose(Gamma), M5)))

    Y[i3] .= Y[i3] .+ M2
    Y[i3] .= Y[i3] .- 1im*w*M4
    Y[i3] .= Y[i3] .+ M5

    Y=convert(Array{Complex{Float64}}, Y)
    return Y
end

function prod_real_complex(A,x)
    # A is a N x N real matrix and x is a complex matrix

    N=size(A,1);
    y=zeros(ComplexF64 , N, 1)
    y=*(A,real.(x))+1im * *(A,imag.(x))
    return y
end

function prod_complex_real(A,x)
    # A is a N x N complex matrix and x is a real matrix

    N=size(A,1);
    y=zeros(ComplexF64 , N, 1)
    y=*(real.(A),x)+1im * *(imag.(A),x)
    return y
end


# function precond_3_3_vector(F,invZ,invP,A,Gamma,w,X1,X2,X3, resProd)

#     n1=length(X1)
#     n2=length(X2)
#     n3=length(X3)

#     i1=range(1, stop=n1)
#     i2=range(n1+1,stop=n1+n2)
#     i3=range(n1+n2+1,stop=n1+n2+n3)

#     Y=zeros(ComplexF64 , n1+n2+n3)
    
#     invZ_view = @view resProd[1:size(invZ,1)]
#     mul!(invZ_view, invZ, X1)
#     Yi1 = @view Y[i1]
#     Y[i1] .= Yi1 .+ invZ_view
   
#     A_view = @view resProd[size(resProd,1)-size(A,2)+1:end]
#     mul!(A_view, transpose(A), invZ_view)
#     M2 = F\A_view
    
#     invP_view = @view resProd[1:size(invP,1)]
#     mul!(invP_view, invP, X2)
#     Yi2 = @view Y[i2]
#     Y[i2] .= Yi2 .+ invP_view
    
#     Gamma_view = @view resProd[size(resProd,1)-size(Gamma,1)+1:end]
#     mul!(Gamma_view, Gamma, invP_view)
#     M4 = F\Gamma_view
    
#     M5 = F\X3
    

#     # Yi1 = @view Y[i1]  
#     A_view = @view resProd[1:size(A,1)]
#     mul!(A_view, A, M2)
#     invZ_view = @view resProd[size(resProd,1)-size(invZ,1)+1:end]
#     mul!(invZ_view, invZ, A_view)
#     Y[i1] .= Yi1 .-lmul!(1.0, invZ_view)
#     mul!(A_view, A, M4)
#     mul!(invZ_view, invZ, A_view)
#     Y[i1] .= Yi1 .+ lmul!(1im*w, invZ_view) 
#     mul!(A_view, A, M5)
#     mul!(invZ_view, invZ, A_view)
#     Y[i1] .= Yi1 .- lmul!(1.0,invZ_view) 

#     # Yi2 = @view Y[i2]
#     Gamma_view = @view resProd[size(resProd,1)-size(Gamma,2)+1:end]
#     mul!(Gamma_view, transpose(Gamma), M2)
#     mul!(invP_view, invP, Gamma_view)
#     Y[i2] .= Yi2 .+ invP_view 
#     mul!(Gamma_view, transpose(Gamma), M4)
#     mul!(invP_view, invP, Gamma_view)
#     Y[i2] .= Yi2 .- lmul!(1im*w,invP_view)
#     mul!(Gamma_view, transpose(Gamma), M5)
#     mul!(invP_view, invP, Gamma_view) 
#     Y[i2] .= Yi2 .+ invP_view
    
#     Yi3 = @view Y[i3]
#     Y[i3] .= Yi3 .+ M2 .- lmul!(1im*w,M4) .+ M5

#     return Y
# end
