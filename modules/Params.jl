# -*- coding: utf-8 -*-
module Params

"Evaluate the 'Redfield' ratios for particulate organic matter, normalised to C."
function redfield(dtPO4_w::Float64, rho_sw::Float64)
    RC = @. 1.0 / (6.9e-3dtPO4_w / 1e-6rho_sw + 6e-3)
    # ^P:C computed as a function of SRP from Galbraith and Martiny PNAS 2015
    RN = 11.0 # value at 60 degS from Martiny et al. Nat G 2013
    RP = 1.0 # Redfield ratio for P in the deep sea
    return 1.0, RN / RC, RP / RC
end  # function redfield

"Evaluate the 'Redfield' ratios for particulate organic matter, normalised to C."
function redfield()
    RC = 106.0
    RN = 16.0
    RP = 1.0
    return 1.0, RN / RC, RP / RC
end  # function redfield

"Calculate the relative molar mass of POM in g/mol."
function rmm_pom(RC::Float64, RN::Float64, RP::Float64)
    rmm_CH2O = 30.031 # g/mol
    rmm_NH3 = 17.031 # g/mol
    rmm_H3PO4 = 97.994 # g/mol
    return RC*rmm_CH2O + RN*rmm_NH3 + RP*rmm_H3PO4
end  # function rmm_pom

"Calculate the depth-dependent porosity."
function phi(phi0::Float64, phiInf::Float64, beta::Float64,
        depths::Array{Float64})
    return @. (phi0 - phiInf)*exp(-beta*depths) + phiInf
end  # function phi

"Calculate the depth-dependent porosity - solid volume fraction."
phiS(phi::Array{Float64}) = 1.0 .- phi

"Calculate the tortuousity squared following Boudreau (1996, GCA)."
tort2(phi::Array{Float64}) = @. 1.0 - 2.0log(phi)

"Calculate the 1st derivative of the depth-dependent porosity w.r.t. depth."
function delta_phi(phi0::Float64, phiInf::Float64, beta::Float64,
        depths::Array{Float64})
    return @. -beta*(phi0 - phiInf)*exp(-beta*depths)
end  # function delta_phi

"""Calculate the 1st derivative of the depth-dependent porosity - solid volume
fraction - w.r.t. depth.
"""
delta_phiS(delta_phi::Array{Float64}) = -delta_phi

"""Calculate the 1st derivative of the inverse of the tortuousity squared
w.r.t. depth.
"""
function delta_tort2i(delta_phi::Array{Float64}, phi::Array{Float64},
        tort2::Array{Float64})
    return @. 2.0delta_phi/(phi*tort2^2)
end  # function delta_tort2i

"depth-dependent tortuosity gain"
function delta_tort2(delta_phi::Array{Float64}, phi::Array{Float64})
    return @. -2.0delta_phi./phi
end

"prepare sediment depth vector"
function prepdepth(depthSed::Float64, z_res::Float64)
    z_range = 0.0:z_res:depthSed
    z = array_of_floats = collect(z_range)
    return z
end

# "Centers z with spacings varying from dz_top to dz_bot across Nz layers."
# function prepdepth_stretched(depthSed; dz_top::Float64=1e-3, dz_bot::Float64=1e-2, Nz::Int=41)
#     @assert dz_top > 0 && dz_bot > 0 && depthSed > 0 && Nz ≥ 2
#     # raw spacings (one per layer)
#     r = range(dz_top, dz_bot; length=Nz) |> collect
#     # scale to hit total depth exactly
#     s = depthSed / sum(r)
#     dx = s .* r
#     # centers; match your legacy style (start at 0.0)
#     z = cumsum(vcat(0.0, dx))[1:end-1]           # length Nz, last < depthSed
#     z[end] = depthSed                            # force exact bottom
#     return z
# end

"""
Geometric stretch from dz_top to ~dz_bot over Nz layers with gentle ratio r.
Returns: (zc, dx, ze, r)
- r is the actual growth factor used (clamped to r_max for smoothness)
"""
function prepdepth_geometric(depthSed; dz_top::Float64=1e-4, dz_bot::Float64=5e-2,
                             Nz::Int=84, r_max::Float64=1.10)
    @assert dz_top>0 && dz_bot>0 && depthSed>0 && Nz≥2
    # target ratio from top/bottom; clamp to keep growth gentle
    r = (dz_bot/dz_top)^(1/(Nz-1))
    r = max(1.0, min(r, r_max))
    # raw geometric spacings and scale to exact depth
    dx_raw = dz_top .* (r .^ (0:Nz-1))
    s      = depthSed / sum(dx_raw)
    dx     = s .* dx_raw
    ze     = cumsum(vcat(0.0, dx))
    zc     = 0.5 .* (ze[1:end-1] .+ ze[2:end])
    return zc, dx, ze, r
end

# # given non-uniform centers z (length Nz)
# function zres_from_z(z::AbstractVector{<:Real})
#     Nz = length(z)
#     zres = similar(z, Float64)
#     zres[1]  = z[2] - z[1]                      # one-sided at the top
#     zres[Nz] = z[Nz] - z[Nz-1]                  # one-sided at the bottom
#     @inbounds for k in 2:Nz-1
#         zres[k] = 0.5 * (z[k+1] - z[k-1])       # center-to-center spacing
#     end
#     return zres
# end


"Assemble depth-dependent porosity parameters."
function porosity(
    phi0::Float64, phiInf::Float64, beta::Float64, depths::Array{Float64},
)
    phi = Params.phi(phi0, phiInf, beta, depths)
    phiS = Params.phiS(phi)  # solid volume fraction
    phiS_phi = phiS./phi
    tort2 = Params.tort2(phi)  # tortuosity squared from Boudreau (1996, GCA)
    delta_phi = Params.delta_phi(phi0, phiInf, beta, depths)
    delta_phiS = Params.delta_phiS(delta_phi)
    delta_tort2i = Params.delta_tort2i(delta_phi, phi, tort2)
    delta_tort2i_tort2 = delta_tort2i.*tort2
    return phi, phiS, phiS_phi, tort2, delta_phi, delta_phiS, delta_tort2i_tort2
end  # function porosity

"Calculate the diffusion-by-bioturbation constant coefficient."
D_bio_0(Fpoc::Float64) = @. 0.0232e-4*(1e2Fpoc)^0.85

"Calculate diffusion by bioturbation vs depth."
function D_bio(depths::Array{Float64}, D_bio_0::Float64, lambda_b::Float64,
        dO2_w::Float64)
    return @. D_bio_0*exp(-(depths/lambda_b)^2)*dO2_w/(dO2_w + 0.02)
end  # function D_bio

"Calculate 1st derivative of diffusion by bioturbation w.r.t. depth."
function delta_D_bio(depths::Array{Float64}, D_bio::Array{Float64},
        lambda_b::Float64)
    return @. -2.0depths*D_bio/lambda_b^2
end  # function delta_D_bio

"Calculate POC degradation parameter."
function krefractory(depths::Array{Float64,1}, D_bio_0::Float64)
    @. 80.25D_bio_0*exp(-depths)
end  # function krefractory

"Calculate fast-degrading POC degradation parameter." 
function kfast(Fpoc::Float64, depths, Q10_primary::Float64, T::Float64, Tref::Float64)
    kfast_0 = 1.5e-1*(1e2*Fpoc)^0.85 * Q10_primary^((T - Tref)/10)
    return fill(kfast_0, depths)
end  # function kfast

"Calculate slow-degrading POC degradation parameter." 
function kslow(Fpoc::Float64, depths, Q10_primary::Float64, T::Float64, Tref::Float64)
    kslow_0 = 1.3e-4*(1e2*Fpoc)^0.85 * Q10_primary^((T - Tref)/10)
    return fill(kslow_0, depths)
end  # function kslow

"Calculate bulk burial velocity at the sediment-water interface in m/a."
x0(Fp::Float64, rho_p::Float64, phiS_2::Float64) = Fp/(rho_p*phiS_2)

"Calculate bulk burial velocity at infinite depth in the sediment in m/a."
xinf(x0::Float64, phiS_2::Float64, phiS_e2::Float64) = x0*phiS_2/phiS_e2

"Calculate porewater burial velocity in m/a."
u(xinf::Float64, phi::Array{Float64}) = xinf*phi[end]./phi

"Calculate solid burial velocity in m/a."
w(xinf::Float64, phiS::Array{Float64}) = xinf*phiS[end]./phiS

# ------------------------
# Half-cell Peclet number
# ------------------------

@inline function Peh(w::AbstractVector{<:Real},
                     z_res::AbstractVector{<:Real},
                     D_bio::AbstractVector{<:Real})
    @assert length(w) == length(z_res) == length(D_bio)
    @. (w * z_res) / (2 * D_bio)
end

@inline function Peh!(out::AbstractVector,
                      w::AbstractVector{<:Real},
                      z_res::AbstractVector{<:Real},
                      D_bio::AbstractVector{<:Real})
    @assert length(out) == length(w) == length(z_res) == length(D_bio)
    @. out = (w * z_res) / (2 * D_bio)
    return out
end

# ------------------------
# sigma(Peh) from Boudreau eq. 96
# sigma = coth(Peh) - 1/Peh
# Use a numerically stable form for small Peh: (x/tanh(x) - 1)/x ~ x/3
# ------------------------

# scalar, numerically-stable
@inline function sigma_scalar(x::Real)
    ax = abs(x)
    if ax < 1e-6
        #  x/3 + x^3/45 + O(x^5)
        return x*(1/3) + x^3*(1/45)
    else
        # (x/tanh(x) - 1)/x  == 1/tanh(x) - 1/x
        return (x / tanh(x) - 1) / x
    end
end

@inline sigma(Peh::AbstractVector{<:Real}) = map(sigma_scalar, Peh)

@inline function sigma!(out::AbstractVector, Peh::AbstractVector{<:Real})
    @assert length(out) == length(Peh)
    @inbounds for i in eachindex(out, Peh)
        out[i] = sigma_scalar(Peh[i])
    end
    return out
end

# convenience wrapper to go directly from (w, z_res, D_bio) → sigma
@inline function sigma(w::AbstractVector{<:Real},
                       z_res::AbstractVector{<:Real},
                       D_bio::AbstractVector{<:Real})
    sigma(Peh(w, z_res, D_bio))
end


"Calculate T-dependent 'free solution' diffusion coeff. for O2 in m^2/a."
D_dO2(T::Float64) = 0.031558 + 0.001428T

"""Calculate T-dependent 'free solution' diffusion coeff. for tCO2 in m^2/a, as
approximated by the bicarbonate diffusion coefficient of Hulse et al (2018).
"""
D_dtCO2(T::Float64) = 0.015179 + 0.000795T

"Nitrate diffusion coefficient in m^2/a."
D_dtNO3(T::Float64) = 0.030863 + 0.001153T

"Sulfate diffusion coefficient  in m^2/a."
D_dtSO4(T::Float64) = 0.015779 + 0.000712T

"Phosphate diffusion coefficient  in m^2/a."
D_dtPO4(T::Float64) = 0.009783 + 0.000513T

"Ammonium diffusion coefficient  in m^2/a."
D_dtNH4(T::Float64) = 0.030926 + 0.001225T

"""Hydrogen sulfide diffusion coefficient in m^2/a."""
D_dtH2S(T::Float64) = 0.028938 + 0.001314T

"Manganese diffusion coefficient in m^2/a."
D_dMn(T::Float64) = 0.009625 + 0.000481T

"Iron diffusion coefficient in m^2/a."
D_dFe(T::Float64) = 0.010761 + 0.000466T

"Bicarbonate diffusion coefficient in m^2/a."
D_dHCO3(T::Float64) = 0.015179 + 0.000795T

"Calcium diffusion coefficient in m^2/a."
D_dCa(T::Float64) = 0.011771 + 0.000529T

"Methane diffusion coefficient  in m^2/a."
D_dCH4(T::Float64)=0.041628576888000035+0.00076369T


"Function is already included in RADIv2 and can be used, but not used in published results"
function disp_coeff(wave_height::Float64, wavelength::Float64, period::Float64, permeability::Float64, depth::Float64, depths::Array{Float64})
    Wnumb = 2π/wavelength
    ang_freq = 2π/period
    return @. 3600*24*365.25*(π*permeability*wave_height/(wavelength*ang_freq^0.5*cosh(Wnumb*depth))).^2*exp(-2*Wnumb.*depths)  
    #[m2/a] dipersion coefficient Eq. 4.146 in Boudreau (1997)
end

"Calculate alpha_0 parameter for irrigation (Archer et al. 2002)."
function alpha_0(Fpoc::Float64, dO2_w::Float64)
    return @. 11.0*(atan((1e2Fpoc*5.0 - 400.0)/400.0)/pi + 0.5) - 0.9 +
        20.0*(dO2_w/(dO2_w + 0.01))*exp(-dO2_w/0.01)*1e2Fpoc/(1e2Fpoc + 30.0)
end  # function alpha_0

"Calculate alpha parameter for irrigation (Archer et al. 2002)."
function alpha(alpha_0::Float64, depths::Array{Float64}, lambda_i::Float64)
    return @. alpha_0*exp(-(depths/lambda_i)^2)
end  # function alpha

"Calculate alpha parameter for irrigation (Archer et al. 2002)."
function alpha(Fpoc::Float64, dO2_w::Float64, depths::Array{Float64},
        lambda_i::Float64)
    return @. alpha_0(Fpoc, dO2_w)*exp(-(depths/lambda_i)^2)
end  # function alpha

"Calculate APPW convenience term."
function APPW(w::Array{Float64}, delta_D_bio::Array{Float64},
        delta_phiS::Array{Float64}, D_bio::Array{Float64}, phiS::Array{Float64})
    return @. w - delta_D_bio - delta_phiS*D_bio/phiS
end  # function APPW

function DFF(tort2, delta_phi, phi, delta_tort2)
    return (tort2 .* delta_phi ./ phi - delta_tort2)./ (tort2.^2)
end

"Calculate TR convenience term."
TR(z_res::Float64, tort2_2::Float64, dbl::Float64) = 2.0 * z_res * tort2_2 / dbl

"Calculate solute flux"
function solute_flux(phi0,diff_coef, v0, vw,dbl)
    Jv = phi0*diff_coef*((v0-vw)/dbl)
    return (Jv)
end

# total sulfate from salinity (mol/kg); same slope you used earlier
@inline TSO4_ref_kg_from_S(S::Float64) = 29264.2e-6 * (S/35.0)  # ≈ 29.264 mmol/kg at S=35

end  # module Params
