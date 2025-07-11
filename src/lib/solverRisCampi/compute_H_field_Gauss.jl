function compute_H_field_Gauss(Lambda_x, Lambda_y, Lambda_z, I)

    mu0 = 4 * pi * 1e-7
    eps0 = 8.854187816997944e-12

    Hx = (1.0 / mu0) * (Lambda_x * I)
    Hy = (1.0 / mu0) * (Lambda_y * I)
    Hz = (1.0 / mu0) * (Lambda_z * I)

    return Hx, Hy, Hz

end