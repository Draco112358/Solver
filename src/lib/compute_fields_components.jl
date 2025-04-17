using LinearAlgebra

function compute_fields_components(phi, tetha, E_theta, E_phi)
    T = [sin(tetha)*cos(phi) cos(tetha)*cos(phi) -sin(phi);
         sin(tetha)*sin(phi) cos(tetha)*sin(phi)  cos(phi);
         cos(tetha)          -sin(tetha)          0]

         #println(T)

    K = T * [1; 0; 0]
    E = T * [0; E_theta; E_phi]
    E_theta_v = T * [0; E_theta; 0]
    E_phi_v = T * [0; 0; E_phi]

    Hm = norm(E) / (120 * pi)
    E_hat = E / norm(E)
    H_hat = cross(E_hat, K)
    H = Hm * H_hat

    return E, K, H, E_theta_v, E_phi_v
end