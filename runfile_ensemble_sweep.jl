# runfile_ensemble_sweep.jl
# HPC-friendly sweep version of runfile_ensemble.jl.
# Usage:
#   julia --project=. runfile_ensemble_sweep.jl [IC_FILE] [ABSTOL] [RELTOL]
# Defaults:
#   IC_FILE = setup/IC_CS2_2_shallow.jl
#   ABSTOL  = 1e-7
#   RELTOL  = 1e-5

# Parse command-line arguments (all optional with defaults)
const _ic_file  = length(ARGS) >= 1 ? ARGS[1]                  : "setup/IC_CS2_2_shallow.jl"
const abstol_v  = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1e-7
const reltol_v  = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 1e-5

println("[sweep] IC file : ", _ic_file)
println("[sweep] abstol  : ", abstol_v)
println("[sweep] reltol  : ", reltol_v)
println("[sweep] threads : ", Threads.nthreads())
flush(stdout)

Threads.nthreads()

# %% [markdown]
# ### Load packages for for plotting and the Differential Equation solver

# %%
using Distributed
@everywhere using MAT
@everywhere using DifferentialEquations, ProgressLogging
@everywhere using Plots; gr()
@everywhere using JLD2
using StaticArrays

# %% [markdown]
# ### RADI modules

# %%
include("modules/gsw_rho.jl")
include("modules/CO2System.jl")
include("modules/React.jl")
include("modules/Equilibrate.jl")
include("modules/Params.jl");

# %% [markdown]
# ### Initial conditions

# %%
# @everywhere include("setup/IC_IberianMargin_shallow.jl");
# @everywhere include("setup/IC_HF2_shallow.jl");
@everywhere include(_ic_file);

# %% [markdown]
# ### Functions running inside the ODE Solver

# %%
# ---------------- SOLUTES ----------------

"Diffusion at SWI (water–sediment interface) for a solute."
@inline function diffuse_SWI(
    then_z1p::Float64,  # u[k+1]
    then_z::Float64,    # u[k]
    then_w::Float64,    # water-side value
    D_var::Float64,
    TR::Float64,
    z_res::Float64,     # vertical resolution
)::Float64
    # (2*(u[k+1]-u[k]) + TR*(uw - u[k])) / dz^2
    return D_var * (2.0*(then_z1p - then_z) + TR*(then_w - then_z)) / (z_res^2)
end

@inline function G_series(Dw::Float64, Deff::Float64, phi0::Float64,
                          dbl::Float64, dz::Float64)::Float64
    δ  = dbl
    Δ  = 0.5 * dz
    ϕ0 = phi0
    return 1.0 / (δ / Dw + Δ / (ϕ0 * Deff))   # [m/yr]
end

"""
Top-cell diffusion tendency with clean 2-resistor SWI exchange.
then_z1p = u[k+1] (cell 2), then_z = u[k] (cell 1), then_w = Cw
Deff0 = porewater effective diffusivity in top cell (e.g. mp.D_dalk_tort2[1])
Dw    = free-water molecular diffusivity in DBL (e.g. mp.D_dalk)
dbl   = DBL thickness for that species (e.g. mp.dbl_dalk)
phi0  = porosity at top cell (mp.phi[1])
dz    = mp.z_resvec[1]
Returns mol m^-3 yr^-1 contribution (same units as your RHS).
"""
@inline function diffuse_SWI_2R(
    then_z1p::Float64, then_z::Float64, then_w::Float64,
    Deff0::Float64, Dw::Float64, dbl::Float64,
    phi0::Float64,
    dz_cell::Float64,   # cell 1 thickness (control volume scaling)
    dz12::Float64,      # center-to-center distance between cell 1 and 2
)::Float64
    G = G_series(Dw, Deff0, phi0, dbl, dz_cell)                 # still uses Δ = dz_cell/2
    interior = 2.0 * Deff0 * (then_z1p - then_z) / (dz12^2)     # use center-to-center spacing
    swi      = 2.0 * G * (then_w - then_z) / (phi0 * dz_cell)
    return interior + swi
end


"Centered advection–diffusion drift term for a solute (interior)."
@inline function advectsolute(
    then_z1p::Float64,  # u[k+1]
    then_z1m::Float64,  # u[k-1]
    u_z::Float64,
    D_var::Float64,
    DFF::Float64,
    z_res::Float64,
)::Float64
    # -(u_bur - D*DFF) * (u[k+1]-u[k-1]) / (2*dz)
    return -(u_z - D_var*DFF) * (then_z1p - then_z1m) / (2.0*z_res)
end

"Advection term at SWI for a solute."
@inline function advectsolute_SWI(
    u_z::Float64,
    then_z::Float64,
    then_w::Float64,
    D_var::Float64,
    DFF::Float64,
    TR::Float64,
    z_res::Float64,
)::Float64
    # -(u_bur - D*DFF) * (-TR*(uw - u[k])) / (2*dz)
    return  (u_z - D_var*DFF) * TR * (then_w - then_z) / (2.0*z_res)
end

"""
SWI transport boundary contribution (the term that used to use TR/(2dz)*(Cw-C1)),
but using the *2-resistor* SWI conductance to get a consistent pore-side SWI gradient.

Returns mol m^-3 yr^-1 (same units as RHS tendencies).

u_z  : burial/advective velocity at top cell [m/yr]
Dw   : free-water molecular diffusivity in DBL [m^2/yr] (e.g. mp.D_dO2)
Deff0: porewater effective diffusivity in top cell [m^2/yr] (e.g. mp.D_dO2_tort2[1])
phi0 : porosity at top cell
dbl  : DBL thickness [m]
dz   : top cell thickness [m]
DFF  : your dispersion/throughflow factor (as used before)
"""
@inline function advectsolute_SWI_2R(
    u_z::Float64,
    C1::Float64,
    Cw::Float64,
    Dw::Float64,
    Deff0::Float64,
    DFF::Float64,
    phi0::Float64,
    dbl::Float64,
    dz::Float64,
)::Float64
    G = G_series(Dw, Deff0, phi0, dbl, dz)               # [m/yr]
    dCdz_swi = (G / (phi0 * Deff0)) * (Cw - C1)          # [mol/m^4]
    return (u_z - Dw * DFF) * dCdz_swi                   # [mol/m^3/yr]
end


"Irrigation exchange for a solute."
@inline function irrigate(then_z::Float64, above::Float64, alpha_z::Float64)::Float64
    return alpha_z * (above - then_z)
end

"Centered diffusion (interior) for a solute/solid."
@inline function diffuse(
    then_z1m::Float64,  # u[k-1]
    then_z::Float64,    # u[k]
    then_z1p::Float64,  # u[k+1]
    D_var::Float64,
    z_res::Float64,     # vertical resolution
)::Float64
    # D * (u[k+1] - 2u[k] + u[k-1]) / dz^2
    return D_var * (then_z1p - 2.0*then_z + then_z1m) / (z_res^2)
end

"Neumann (zero-gradient) diffusion at bottom boundary for a solute/solid."
@inline function diffuse_BBC(
    then_z::Float64,     # u[end]
    then_z1m::Float64,   # u[end-1]
    D_var::Float64,
    z_res::Float64,     # vertical resolution
)::Float64
    # 2D * (u[end-1] - u[end]) / dz^2
    return 2.0 * D_var * ((then_z1m - then_z) / (z_res^2))
end

# ---------------- SOLIDS ----------------

"Diffusion at SWI for a solid (includes throughflow term)."
@inline function diffuseSolid_SWI(
    then_z1p::Float64,
    then_z::Float64,
    F::Float64,          # flux
    D_bio::Float64,
    phiS::Float64,
    w::Float64,
    z_res::Float64,     # vertical resolution
)::Float64
    # D_bio * ( 2(u[k+1]-u[k]) + 2*dz * (F - phiS*w*u[k]) / (D_bio*phiS) ) / dz^2
    # -> D_bio cancels in last term’s numerator/denominator; keep form for clarity
    return D_bio * ( 2.0*(then_z1p - then_z) +
                     (2.0*z_res) * (F - phiS*w*then_z) / (D_bio*phiS) ) / (z_res^2)
end

"Advection at SWI for a solid."
@inline function advectsolid_SWI(
    then_z::Float64,
    then_z1p::Float64,
    F::Float64,
    D_bio::Float64,
    APPW::Float64,
    sigma1m::Float64,
    sigma::Float64,
    sigma1p::Float64,
    phiS::Float64,
    w::Float64,
    z_res::Float64,     # vertical resolution
)::Float64
    # -APPW * [ σ1m*u[k+1] + 2σ*u[k] - σ1p*(u[k+1] + 2*dz/D_bio*(F/φS - w*u[k])) ] / (2dz)
    return -APPW * ( sigma1m*then_z1p + 2.0*sigma*then_z -
                     sigma1p*( then_z1p + (2.0*z_res/D_bio)*(F/phiS - w*then_z) )
                   ) / (2.0*z_res)
end

"Centered advection for a solid (interior)."
@inline function advectsolid(
    then_z::Float64,     # u[k]
    then_z1p::Float64,   # u[k+1]
    then_z1m::Float64,   # u[k-1]
    APPW_z::Float64,
    sigma_z::Float64,
    sigma1p_z::Float64,
    sigma1m_z::Float64,
    z_res::Float64,      # vertical resolution
)::Float64
    # -APPW * ( σ1m*u[k+1] + 2σ*u[k] - σ1p*u[k-1] ) / (2dz)
    return -APPW_z * ( sigma1m_z*then_z1p + 2.0*sigma_z*then_z - sigma1p_z*then_z1m ) / (2.0*z_res)
end

"Bottom boundary advection for a solid."
@inline function advectsolid_BBC(
    then_z::Float64,    # u[end]
    then_z1m::Float64,  # u[end-1]
    APPW_z::Float64,
    sigma_z::Float64,
    z_res::Float64,  # vertical resolution
)::Float64
    # -APPW * σ * (u[end] - u[end-1]) / dz
    return -APPW_z * sigma_z * (then_z - then_z1m) / z_res
end


function input_P(omegaCa, T)
    P = (0.9012882388755719 + 0.01814331*T - log(omegaCa)) / 0.00016455
    return P
end

function input_omegaCa(P, T)
    omegaCa = exp(0.9012882388755719 + 0.01814331*T - 0.00016455*P)
    return omegaCa
end


"Function to calculate the DBL thickness based on the current speed, following: Sulpis, O., Boudreau, B. P., Mucci, A., Jenkins, C., Trossman, D. S., Arbic, B. K., & Key, R. M. (2018). Current CaCO3 dissolution at the seafloor caused by anthropogenic CO2. Proceedings of the National Academy of Sciences, 115(46), 11700-11705."
function DBL(U, T; constant_DBL=false, user_DBL=nothing)
    # temperature dependent friction velocity
    u_star = 0.00136 - 2.19598542e-5*T + 2.35862843e-2 * U
    # Kinematic viscosity 
    nu = 1.747567451381780806e-6 - 3.23886387e-8*T

    # Helper function for diffusion and Schmidt number calculations
    function calc_diffusion_and_schmidt(param_func)
        D = param_func(T) / (60 * 60 * 24 * 365.25)
        Sc = nu / D
        beta = 0.0417 * u_star * Sc^(-2/3)
        return D / beta
    end

    # List of parameter functions
    param_funcs = [Params.D_dO2, Params.D_dtCO2, Params.D_dtNO3, Params.D_dtSO4, Params.D_dtPO4, Params.D_dtNH4, Params.D_dtH2S, Params.D_dMn, Params.D_dFe, Params.D_dCH4, Params.D_dHCO3, Params.D_dCa]

    if constant_DBL
        # Apply the user-provided DBL value to all species
        return fill(user_DBL, length(param_funcs))
    else
        # Apply the helper function to each parameter function
        return [calc_diffusion_and_schmidt(func) for func in param_funcs]
    end
end


"apparent bulk sediment diffusivity, from: McGinnis, D. F., Sommer, S., Lorke, A., Glud, R. N., & Linke, P. (2014). Quantifying tidally driven benthic oxygen exchange across permeable sediments: An aquatic eddy correlation study. Journal of Geophysical Research: Oceans, 119(10), 6918-6932."
function Ksed(U, T, permeability)
    u_star=0.00136 - 2.19598542e-5*T + 2.35862843e-2 * U #mm/s, from: Sulpis, O., Boudreau, B. P., Mucci, A., Jenkins, C., Trossman, D. S., Arbic, B. K., & Key, R. M. (2018). Current CaCO3 dissolution at the seafloor caused by anthropogenic CO2. Proceedings of the National Academy of Sciences, 115(46), 11700-11705.
    if permeability >= 1e-12 # Threshold value identified in [Huettel et al., 2014]
        alpha = 1.0
    else
        alpha = 0.0
    end
    von_k = 0.4
    # Helper to calculate Temperature-Dependent diffusion, accounting for tortuosity 
    function diff_tort(param_func)
            D = param_func(T) / (60 * 60 * 24 * 365.25)
        return D*exp(alpha*von_k*(u_star*1000))*60*60*24*365.25
    end
    
    param_funcs = [Params.D_dO2, Params.D_dtCO2, Params.D_dtNO3, Params.D_dtSO4, Params.D_dtPO4, Params.D_dtNH4, Params.D_dtH2S, Params.D_dMn, Params.D_dFe, Params.D_dCH4, Params.D_dHCO3, Params.D_dCa]
    
    return [diff_tort(func) for func in param_funcs]
end

# --- tuning knobs ---
const HMIN = 1e-12     # clamp for H+ concentration (≈ pH 12 upper bound for [H+])
const HMAX = 1.0       # clamp for H+ (≈ pH 0 lower bound)
const ETA_ALK  = 1e-6  # relative tol to skip re-solving H
const DALK_EPS = 1e-18 # guard for tiny derivatives

# convert helpers
@inline to_kg(x, ρ) = x/ρ
@inline to_m3(x, ρ) = x*ρ

# CO2-system constants on *kg* basis for Equilibrate.*
@inline function chem_constants_kg(mp)
    ρ = mp.rho_sw
    return (
        K1  = mp.K1/ρ,   K2  = mp.K2/ρ,   Kw  = mp.Kw/(ρ*ρ),
        KB  = mp.KB/ρ,   KF  = mp.KF/ρ,   KSO4= mp.KSO4/ρ,
        KP1 = mp.KP1/ρ,  KP2 = mp.KP2/ρ,  KP3 = mp.KP3/ρ,
        KSi = mp.KSi/ρ,  KNH3= mp.KNH3/ρ, KH2S= mp.KH2S/ρ,
        TB  = mp.TB/ρ,   TF  = mp.TF/ρ,
        # (KCa, KAr are only used in your precipitation kernel; keep them on m³ basis there)
    )
end

# H_from_alk_m3: all inputs/outputs in mol/m^3; internally switches to kg for Equilibrate.*
@inline function H_from_alk_m3(
    TA_m3::Float64, tCO2_m3::Float64, tNH4_m3::Float64, tPO4_m3::Float64,
    tSO4_m3::Float64, tH2S_m3::Float64, mp;
    H0::Float64 = 1e-8 * mp.rho_sw,   # ~pH 8 guess on m^3 basis
    iters::Int = 1,                   # like your original: 1 (or 2) steps
    damp::Float64 = 1.0
)::Float64
    ρ  = mp.rho_sw
    cs = chem_constants_kg(mp)

    # convert state/species to kg-basis for Equilibrate.*
    TA   = to_kg(TA_m3, ρ)
    tCO2 = to_kg(tCO2_m3, ρ)
    tNH4 = to_kg(tNH4_m3, ρ)
    tPO4 = to_kg(tPO4_m3, ρ)
    tSO4 = to_kg(tSO4_m3, ρ)
    tH2S = to_kg(tH2S_m3, ρ)
    dSi  = to_kg(mp.dSi_inp, ρ)

    # kg-basis H guess + clamp in reasonable [1e-12, 1] kg-range
    H = clamp(to_kg(H0, ρ), 1e-12, 1.0)

    @inbounds for _ = 1:iters
        # total alkalinity pieces on kg-basis
        alk_bor = Equilibrate.alk_borate(H, cs.TB, cs.KB)
        alk_nc  = Equilibrate.alk_ammonia(H, tNH4, cs.KNH3) + alk_bor +
                  Equilibrate.alk_fluoride(H, cs.TF, cs.KF) +
                  Equilibrate.alk_phosphate(H, tPO4, cs.KP1, cs.KP2, cs.KP3) +
                  Equilibrate.alk_silicate(H, dSi,  cs.KSi) +
                  Equilibrate.alk_sulfate(H, tSO4, cs.KSO4) +
                  Equilibrate.alk_sulfide(H, tH2S, cs.KH2S) +
                  Equilibrate.alk_water(H, cs.Kw)
        alk_c   = Equilibrate.alk_carbonate(H, tCO2, cs.K1, cs.K2)

        res = TA - (alk_c + alk_nc)  # kg-basis residual

        # early exit if close (same relative criterion but on kg-basis)
        if abs(res) <= ETA_ALK * max(1.0, abs(TA))
            break
        end

        # dAlk_dH = Equilibrate.dalk_dh(H, tCO2, alk_bor, cs.K1, cs.K2, cs.KB, cs.Kw)
        dAlk_dH = Equilibrate.dalk_dh(
            H, tCO2, alk_bor, tNH4, cs.KNH3, cs.TF, cs.KF, tPO4, cs.KP1, cs.KP2, cs.KP3, dSi, cs.KSi, tSO4, cs.KSO4, tH2S, 
            cs.KH2S, cs.K1, cs.K2, cs.KB, cs.Kw
        )
        if !isfinite(dAlk_dH) || abs(dAlk_dH) < DALK_EPS
            break
        end
        Δ = res / dAlk_dH
        if !isfinite(Δ)
            break
        end
        H = clamp(H + damp*Δ, 1e-12, 1.0)
    end

    return to_m3(H, ρ)  # back to mol/m^3 for your RHS
end

@inline function rates_wD!(
    r::AbstractVector{Float64}, mp::MP, k::Int,
    H::Float64,
    dO2::Float64, dtCO2::Float64, dtNO3::Float64, dtSO4::Float64, dtPO4::Float64,
    dtNH4::Float64, dtH2S::Float64, dFeII::Float64, dMnII::Float64, dCH4::Float64,
    dalk::Float64, dCa::Float64, pfoc::Float64, psoc::Float64, pFeOH3::Float64,
    pMnO2::Float64, pcalcite::Float64, paragonite::Float64, pFeS::Float64, pFeS2::Float64,
    pS0::Float64, pFeOH3_PO4::Float64
)::Nothing where {MP}
    @inbounds begin
        # ---- pull *typed* params from mp (index vector fields at k) ----
        kfast_k    = (mp.kfast[k])::Float64
        kslow_k    = (mp.kslow[k])::Float64
        phiS_phi_k = (mp.phiS_phi[k])::Float64

        K1   = (mp.K1)::Float64
        K2   = (mp.K2)::Float64
        KCa  = (mp.KCa)::Float64
        KAr  = (mp.KAr)::Float64
        RC   = (mp.RC)::Float64
        RN   = (mp.RN)::Float64
        RP   = (mp.RP)::Float64

        # ---- cheap carbonate helper you had inline ----
        denom::Float64 = K1*K2 + K1*H + H*H
        dCO3::Float64  = (dtCO2 * K1 * K2) / denom

        # ---- call your reaction kernel (must be scalar & type-stable) ----
        tup = React.rates(
            dO2, dtNO3, pMnO2, pFeOH3, dtSO4, dtNH4, dtH2S, dFeII, dMnII, dCH4, dtPO4, pFeOH3_PO4, pFeS, pS0, pFeS2,
            pfoc*kfast_k, psoc*kslow_k, pcalcite, paragonite, dCa,
            dCO3, KCa, KAr, phiS_phi_k, RC, RN, RP, T, calcite_diss_scheme, aragonite_diss_scheme, Q10_secondary, Tref
        )  # should be NTuple{27,Float64}

        # destructure without allocations
        rate_dO2, rate_dtCO2, rate_dtNO3, rate_dtSO4, rate_dtPO4, rate_dtNH4, rate_dtH2S,
        rate_dFeII, rate_dMnII, rate_dCH4, rate_dalk, rate_dCa, rate_pfoc, rate_psoc,
        rate_pFeOH3, rate_pMnO2, rate_pFeS, rate_pFeS2, rate_pFeOH3_PO4, rate_S0, rate_pcalcite, rate_paragonite,
        Rdeg_dO2, Rdeg_dtNO3, Rdeg_dtSO4, Rdeg_pFeOH3, Rdeg_pMnO2, Rdeg_dCH4, Rdeg_total = tup

        # write directly to r (no broadcasts)
        r[1]  = rate_dO2
        r[2]  = rate_dtCO2
        r[3]  = rate_dtNO3
        r[4]  = rate_dtSO4
        r[5]  = rate_dtPO4
        r[6]  = rate_dtNH4
        r[7]  = rate_dtH2S
        r[8]  = rate_dFeII
        r[9]  = rate_dMnII
        r[10] = rate_dCH4
        r[11] = rate_dalk
        r[12] = rate_dCa
        r[13] = rate_pfoc
        r[14] = rate_psoc
        r[15] = rate_pFeOH3
        r[16] = rate_pMnO2
        r[17] = rate_pcalcite
        r[18] = rate_paragonite
        r[19] = rate_pFeS
        r[20] = rate_pFeS2
        r[21] = rate_S0
        r[22] = rate_pFeOH3_PO4
        # r[19] presumably unused in your layout
        # r[23] = Rdeg_dO2
        # r[24] = Rdeg_dtNO3
        # r[25] = Rdeg_dtSO4
        # r[26] = Rdeg_pFeOH3
        # r[27] = Rdeg_pMnO2
        # r[28] = Rdeg_dCH4
        # r[29] = Rdeg_total
    end
    return nothing
end

@inline function make_H_cache(dH_i, Nz::Int)
    if dH_i isa AbstractVector
        if length(dH_i) == Nz
            return Float64.(dH_i)            # copy, ensure Float64
        else
            # fallback: use first value everywhere
            return fill(Float64(dH_i[1]), Nz)
        end
    else
        return fill(Float64(dH_i), Nz)
    end
end

@inline nz(x::T) where {T<:Real} = ifelse(x < zero(T), zero(T), x);

# %% [markdown]
#

# %% [markdown]
# ### Parameters needed to run the model

# %%
@everywhere struct ModelParams
    zc
    phi
    phiS
    phiS_phi
    tort2
    delta_phi
    delta_phiS
    delta_tort2i_tort2
    rho_sw
    RC
    RN
    RP
    Mpom
    Fpom_mol
    Fpoc
    Ffoc
    Fsoc
    Froc
    Fcalcite_out
    M_MnO2
    M_FeOH3 
    M_CaCO3
    M_clay
    Fp
    D_bio_0
    D_bio
    delta_D_bio
    krefractory
    kfast
    kslow
    x0
    xinf
    u_bur
    w
    Peh
    sigma
    sigma1m
    sigma1p
    D_dO2
    D_dtCO2
    D_dtNO3
    D_dtSO4
    D_dtPO4
    D_dtNH4
    D_dtH2S
    D_dMnII
    D_dFeII
    D_dCH4
    D_dalk
    D_dCa
    D_dO2_tort2
    D_dtCO2_tort2
    D_dtNO3_tort2
    D_dtSO4_tort2
    D_dtPO4_tort2
    D_dtNH4_tort2
    D_dtH2S_tort2
    D_dMnII_tort2
    D_dFeII_tort2
    D_dCH4_tort2
    D_dalk_tort2
    D_dCa_tort2
    alpha_0
    alpha
    dbl_dO2
    dbl_dtCO2
    dbl_dtNO3
    dbl_dtSO4
    dbl_dtPO4
    dbl_dtNH4
    dbl_dtH2S
    dbl_dMnII
    dbl_dFeII
    dbl_dCH4
    dbl_dalk
    dbl_dCa
    APPW
    DFF
    TR_dO2
    TR_dtCO2
    TR_dtNO3
    TR_dtSO4
    TR_dtPO4
    TR_dtNH4
    TR_dtH2S
    TR_dMnII
    TR_dFeII
    TR_dCH4
    TR_dalk
    TR_dCa
    zr_Db_0
    K1
    K2
    Kw
    KB
    KF
    KSO4
    KP1
    KP2
    KP3
    KSi
    KNH3
    KH2S
    TB
    TF
    KCa
    KAr
    dH_i
    dSi_inp
    omegaCa
    rates_scratch
    H_cache
    H_diag
    z_resvec
    z_rescc
end

# %%
@everywhere begin   
    function calculate_constants(T, U, P, Fpom, Fcalcite, Fpom_s, Fpom_f, kfast, kslow, dO2_w_local)
        # sediment depth vector
        zc, dx, ze, r = Params.prepdepth_geometric(depthSed; dz_top=dz_top, dz_bot=dz_bot, Nz=Nz, r_max=1.10)
        z_resvec = dx

        # Calculate depth-dependent porosity
        phi, phiS, phiS_phi, tort2, delta_phi, delta_phiS, delta_tort2i_tort2 =
            Params.porosity(phi0, phiInf, beta, zc)
        # Define 'Redfield' ratios and OM stoichiometry
        rho_sw = gsw_rho(S, T, P)  # seawater density [kg/m^3]
        dSi_inp = dSi_w*rho_sw
        # RC, RN, RP = Params.redfield(dtPO4_w, rho_sw)  # for P-variable ratios
        RC, RN, RP = Params.redfield()  # for constant, canonical Redfield values
        Mpom = Params.rmm_pom(RC, RN, RP)  # g/mol
        Fpom_mol = Fpom / Mpom  # mol/m^2/a
        Fpoc = Fpom_mol * RC  # mol/m^2/a
        # Split total flux into fast-slow-refractory portions
        Ffoc = Fpoc * Fpom_f #mol
        Fsoc = Fpoc * Fpom_s #mol
        Froc = Fpoc * Fpom_r #mol
        if !(Fpom_f + Fpom_s + Fpom_r ≈ 1.0)
            println("\nRadi WARNING: the fractions of POM do not add up to 1!\n")
        end
        Fcalcite_out = Fcalcite  # mol/m^2/a
        # `Fp` = total sediment flux to bottom in g/m^2/a
        M_MnO2 = 86.9368  # g/mol
        M_FeOH3 = 106.867  # g/mol
        M_CaCO3 = 100.0869  # g/mol
        M_clay = 360.31  # g/mol (montmorillonite)
        Fp = Fpom + FMnO2 * M_MnO2 + FFeOH3 * M_FeOH3 + (Fcalcite_out + Faragonite) * M_CaCO3 + Fclay * M_clay

        # Bioturbation (for solids)
        D_bio_0 = Params.D_bio_0(Fpoc)
        # ^[m2/a] surf bioturb coeff, Archer et al. (2002)
        D_bio = Params.D_bio(zc, D_bio_0, lambda_b, dO2_w_local*rho_sw)
        # ^[m2/a] bioturb coeff, Archer et al (2002)
        delta_D_bio = Params.delta_D_bio(zc, D_bio, lambda_b)

        # Organic matter degradation parameters
        krefractory = 0.0
        kfast = fill(kfast, Nz)
        kslow = fill(kslow, Nz)
        # ^[/a] from Archer et al (2002)

        # Solid fluxes and solid initial conditions
        x0 = Params.x0(Fp, rho_p, phiS[1])
        # ^[m/a] bulk burial velocity at sediment-water interface
        xinf = Params.xinf(x0, phiS[1], phiS[end])
        u_bur = Params.u(xinf, phi)  # [m/a] porewater burial velocity
        w = Params.w(xinf, phiS)  # [m/a] solid burial velocity

        # Biodiffusion depth-attenuation: see Boudreau (1996); Fiadeiro & Veronis (1977)
        Peh = Params.Peh(w, z_resvec, D_bio)
        # ^one half the cell Peclet number (Eq. 97 in Boudreau 1996)
        # When Peh<<1, biodiffusion dominates, when Peh>>1, advection dominates
        sigma = Params.sigma(Peh)
        sigma1m = 1.0 .- sigma
        sigma1p = 1.0 .+ sigma

        D_dO2 =  Params.D_dO2(T)
        D_dtCO2 = Params.D_dtCO2(T)
        D_dtNO3 =  Params.D_dtNO3(T)
        D_dtSO4 =  Params.D_dtSO4(T)
        D_dtPO4 =  Params.D_dtPO4(T)
        D_dtNH4 =  Params.D_dtNH4(T)
        D_dtH2S =  Params.D_dtH2S(T)
        D_dMnII =  Params.D_dMn(T)
        D_dFeII =  Params.D_dFe(T)
        D_dCH4 = Params.D_dCH4(T)
        D_dalk = Params.D_dHCO3(T)
        D_dCa =  Params.D_dCa(T)

        D_dO2_tort2 = Ksed(U, T, permeability)[1] ./ tort2
        D_dtCO2_tort2 = Ksed(U, T, permeability)[2] ./ tort2
        D_dtNO3_tort2 = Ksed(U, T, permeability)[3] ./ tort2
        D_dtSO4_tort2 =  Ksed(U, T, permeability)[4] ./ tort2
        D_dtPO4_tort2 = Ksed(U, T, permeability)[5] ./ tort2
        D_dtNH4_tort2 = Ksed(U, T, permeability)[6] ./ tort2
        D_dtH2S_tort2 = Ksed(U, T, permeability)[7] ./ tort2
        D_dMnII_tort2 = Ksed(U, T, permeability)[8] ./ tort2
        D_dFeII_tort2 = Ksed(U, T, permeability)[9] ./ tort2
        D_dCH4_tort2 = Ksed(U, T, permeability)[10] ./ tort2
        D_dalk_tort2 = Ksed(U, T, permeability)[11] ./ tort2
        D_dCa_tort2 = Ksed(U, T, permeability)[12] ./ tort2

        # Irrigation (for solutes)
        alpha_0 = Params.alpha_0(Fpoc, dO2_w_local*rho_sw)  # [/a] from Archer et al (2002)
        alpha = Params.alpha(alpha_0, zc, lambda_i) # [/a] Archer et al (2002)

        # Species and temperature dependent DBL thickness
        dbl_dO2, dbl_dtCO2, dbl_dtNO3, dbl_dtSO4, dbl_dtPO4, dbl_dtNH4, dbl_dtH2S, dbl_dMnII, dbl_dFeII, dbl_dCH4, dbl_dalk, dbl_dCa = DBL(U, T) 

        # Miscellaneous convenience variables
        APPW = Params.APPW(w, delta_D_bio, delta_phiS, D_bio, phiS)
        delta_tort2 = Params.delta_tort2(delta_phi, phi)
        DFF = Params.DFF(tort2, delta_phi, phi, delta_tort2)
        TR_dO2 = Params.TR(z_resvec[1], tort2[1], dbl_dO2)
        TR_dtCO2 = Params.TR(z_resvec[1], tort2[1], dbl_dtCO2)
        TR_dtNO3 = Params.TR(z_resvec[1], tort2[1], dbl_dtNO3)
        TR_dtSO4 = Params.TR(z_resvec[1], tort2[1], dbl_dtSO4)
        TR_dtPO4 = Params.TR(z_resvec[1], tort2[1], dbl_dtPO4)
        TR_dtNH4 = Params.TR(z_resvec[1], tort2[1], dbl_dtNH4)
        TR_dtH2S = Params.TR(z_resvec[1], tort2[1], dbl_dtH2S)
        TR_dMnII = Params.TR(z_resvec[1], tort2[1], dbl_dMnII)
        TR_dFeII = Params.TR(z_resvec[1], tort2[1], dbl_dFeII)
        TR_dCH4 = Params.TR(z_resvec[1], tort2[1], dbl_dCH4)
        TR_dalk = Params.TR(z_resvec[1], tort2[1], dbl_dalk)
        TR_dCa = Params.TR(z_resvec[1], tort2[1], dbl_dCa)
        zr_Db_0 = 2.0 * z_resvec[1] / D_bio[1]

        # Get it all from CO2System.jl instead, with pH all on Free scale
        co2s = CO2System.CO2SYS(
            1e6dalk_i,
            1e6dtCO2_i,
            1,
            2,
            S,
            T,
            T,
            P,
            P,
            1e6dSi_w,
            1e6dtPO4_i,
            1e6dtNH4_i,
            1e6dtH2S_i,
            3,
            10,
            1,)[1]
        K1 = co2s[1, 54][1] * rho_sw
        K2 = co2s[1, 55][1] * rho_sw
        Kw = co2s[1, 58][1] * rho_sw ^ 2
        KB = co2s[1, 59][1] * rho_sw
        KF = co2s[1, 60][1] * rho_sw
        KSO4 = co2s[1, 61][1] * rho_sw
        KP1 = co2s[1, 62][1] * rho_sw
        KP2 = co2s[1, 63][1] * rho_sw
        KP3 = co2s[1, 64][1] * rho_sw
        KSi = co2s[1, 65][1] * rho_sw
        KNH3 = co2s[1, 66][1] * rho_sw
        KH2S = co2s[1, 67][1] * rho_sw
        TB = co2s[1, 83][1] * 1e-6rho_sw
        TF = co2s[1, 84][1] * 1e-6rho_sw
        KCa = co2s[1, 86][1] * rho_sw ^ 2
        KAr = co2s[1, 87][1] * rho_sw ^ 2
        dH_i = @. (10.0 ^ -co2s[:, 35]) * rho_sw
        dH_i = length(dH_i) == 1 ? dH_i[1] : dH_i

        H_cache = make_H_cache(dH_i, Nz)
        H_diag = make_H_cache(dH_i, Nz)

        omegaCa=omegaCa=co2s[1, 30][1]

        rates_scratch = zeros(Float64, 26)

        z_rescc = similar(z_resvec)
        @inbounds for k = 2:Nz-1
            z_rescc[k] = 0.5*(z_resvec[k-1] + z_resvec[k])   # center-to-center distance
        end

        return ModelParams(zc,
        phi,
        phiS,
        phiS_phi,
        tort2,
        delta_phi,
        delta_phiS,
        delta_tort2i_tort2,
        rho_sw,
        RC,
        RN,
        RP,
        Mpom,
        Fpom_mol,
        Fpoc,
        Ffoc,
        Fsoc,
        Froc,
        Fcalcite_out,
        M_MnO2,
        M_FeOH3, 
        M_CaCO3,
        M_clay,
        Fp,
        D_bio_0,
        D_bio,
        delta_D_bio,
        krefractory,
        kfast,
        kslow,
        x0,
        xinf,
        u_bur,
        w,
        Peh,
        sigma,
        sigma1m,
        sigma1p,
        D_dO2,
        D_dtCO2,
        D_dtNO3,
        D_dtSO4,
        D_dtPO4,
        D_dtNH4,
        D_dtH2S,
        D_dMnII,
        D_dFeII,
        D_dCH4,
        D_dalk,
        D_dCa,
        D_dO2_tort2,
        D_dtCO2_tort2,
        D_dtNO3_tort2,
        D_dtSO4_tort2,
        D_dtPO4_tort2,
        D_dtNH4_tort2,
        D_dtH2S_tort2,
        D_dMnII_tort2,
        D_dFeII_tort2,
        D_dCH4_tort2,
        D_dalk_tort2,
        D_dCa_tort2,
        alpha_0,
        alpha,
        dbl_dO2,
        dbl_dtCO2,
        dbl_dtNO3,
        dbl_dtSO4,
        dbl_dtPO4,
        dbl_dtNH4,
        dbl_dtH2S,
        dbl_dMnII,
        dbl_dFeII,
        dbl_dCH4,
        dbl_dalk,
        dbl_dCa,
        APPW,
        DFF,
        TR_dO2,
        TR_dtCO2,
        TR_dtNO3,
        TR_dtSO4,
        TR_dtPO4,
        TR_dtNH4,
        TR_dtH2S,
        TR_dMnII,
        TR_dFeII,
        TR_dCH4,
        TR_dalk,
        TR_dCa,
        zr_Db_0,
        K1,
        K2,
        Kw,
        KB,
        KF,
        KSO4,
        KP1,
        KP2,
        KP3,
        KSi,
        KNH3,
        KH2S,
        TB,
        TF,
        KCa,
        KAr,
        dH_i,
        dSi_inp,
        omegaCa,
        rates_scratch,
        H_cache,
        H_diag,
        z_resvec,
        z_rescc,)
    end
end

# %%
kfast = Params.kfast(Fpom, Nz, Q10_primary, T, Tref)
kslow = Params.kslow(Fpom, Nz, Q10_primary, T, Tref)

model_params=calculate_constants(T, U, P, Fpom, Fcalcite, Fpom_s, Fpom_f, kfast, kslow, dO2_w)

# %%
const ParamsT = typeof(model_params) # or typeof(proto) if you prefer a separate proto

# %% [markdown]
# ### Prepare initial conditions and parameters for the DifferentialEquation.jl solver

# %%
# Create variables to model
#make a depth vector filled with the initial conditions from the setup file
dO2 = fill(dO2_i*model_params.rho_sw, length(model_params.zc))
dtCO2 = fill(dtCO2_i*model_params.rho_sw, length(model_params.zc))
dtNO3 = fill(dtNO3_i*model_params.rho_sw, length(model_params.zc))
dtSO4 = fill(dtSO4_i*model_params.rho_sw, length(model_params.zc))
dtPO4 = fill(dtPO4_i*model_params.rho_sw, length(model_params.zc))
dtNH4 = fill(dtNH4_i*model_params.rho_sw, length(model_params.zc))
dtH2S = fill(dtH2S_i*model_params.rho_sw, length(model_params.zc))
dFeII = fill(dFeII_i*model_params.rho_sw, length(model_params.zc))
dMnII = fill(dMnII_i*model_params.rho_sw, length(model_params.zc))
dCH4 = fill(dCH4_i*model_params.rho_sw, length(model_params.zc))
dalk = fill(dalk_i*model_params.rho_sw, length(model_params.zc))
dCa_w = 0.02128 / 40.087 * S / 1.80655 #mol/kg
dCa = fill(dCa_w*model_params.rho_sw, length(model_params.zc))
pfoc = fill(pfoc_i, length(model_params.zc))
psoc = fill(psoc_i, length(model_params.zc))
proc = fill(proc_i, length(model_params.zc))
pFeOH3 = fill(pFeOH3_i, length(model_params.zc))
pMnO2 = fill(pMnO2_i, length(model_params.zc))
pcalcite = fill(pcalcite_i, length(model_params.zc))
paragonite = fill(paragonite_i, length(model_params.zc))
pFeS = fill(pFeS_i, length(model_params.zc))
pFeS2 = fill(pFeS2_i, length(model_params.zc))
pS0 = fill(pS0_i, length(model_params.zc))
pFeOH3_PO4 = fill(pFeOH3_PO4_i, length(model_params.zc))
dH = fill(model_params.dH_i, length(model_params.zc))

#create an input matrix for the solver
u0 = zeros(22, length(model_params.zc))
u0[1, :] = dO2
u0[2, :] = dtCO2
u0[3, :] = dtNO3
u0[4, :] = dtSO4
u0[5, :] = dtPO4
u0[6, :] = dtNH4
u0[7, :]= dtH2S
u0[8, :]= dFeII
u0[9, :]= dMnII
u0[10, :]= dCH4
u0[11, :]=dalk
u0[12, :]=dCa
u0[13, :]=pfoc
u0[14, :]=psoc 
u0[15, :]=pFeOH3
u0[16, :]=pMnO2
u0[17, :]=pcalcite
u0[18, :]=paragonite
u0[19, :]=pFeS
u0[20, :]=pFeS2
u0[21, :]=pS0
u0[22, :]=pFeOH3_PO4

# %%
const RHS_CALLS = Threads.Atomic{Int}(0)
const H_CALLS   = Threads.Atomic{Int}(0)

# %% [markdown]
# ### Differential equation solver

# %%
function physics_ensamble!(du, u , p ,t)

    mp = (p.model_params)::ParamsT     # <- one cast: everything below is inferred
    dO2_w_local = p.dO2_w::Float64

    # --- aliases (typed) ---
    rho_sw  = mp.rho_sw::Float64
    alpha   = mp.alpha::Vector{Float64}
    DFF     = mp.DFF::Vector{Float64}
    u_bur   = mp.u_bur::Vector{Float64}
    phiS    = mp.phiS::Vector{Float64}
    sigma   = mp.sigma::Vector{Float64}
    sigma1m = mp.sigma1m::Vector{Float64}
    sigma1p = mp.sigma1p::Vector{Float64}
    D_bio   = mp.D_bio::Vector{Float64}
    w       = mp.w::Vector{Float64}

    D_dO2_tort2   = mp.D_dO2_tort2::Vector{Float64}
    D_dtCO2_tort2 = mp.D_dtCO2_tort2::Vector{Float64}
    D_dtNO3_tort2 = mp.D_dtNO3_tort2::Vector{Float64}
    D_dtSO4_tort2 = mp.D_dtSO4_tort2::Vector{Float64}
    D_dtPO4_tort2 = mp.D_dtPO4_tort2::Vector{Float64}
    D_dtNH4_tort2 = mp.D_dtNH4_tort2::Vector{Float64}
    D_dtH2S_tort2 = mp.D_dtH2S_tort2::Vector{Float64}
    D_dFeII_tort2 = mp.D_dFeII_tort2::Vector{Float64}
    D_dMnII_tort2 = mp.D_dMnII_tort2::Vector{Float64}
    D_dCH4_tort2  = mp.D_dCH4_tort2::Vector{Float64}
    D_dalk_tort2  = mp.D_dalk_tort2::Vector{Float64}
    D_dCa_tort2   = mp.D_dCa_tort2::Vector{Float64}

    D_dO2   = mp.D_dO2::Float64
    D_dtCO2 = mp.D_dtCO2::Float64
    D_dtNO3 = mp.D_dtNO3::Float64
    D_dtSO4 = mp.D_dtSO4::Float64
    D_dtPO4 = mp.D_dtPO4::Float64
    D_dtNH4 = mp.D_dtNH4::Float64
    D_dtH2S = mp.D_dtH2S::Float64
    D_dFeII = mp.D_dFeII::Float64
    D_dMnII = mp.D_dMnII::Float64
    D_dCH4  = mp.D_dCH4::Float64
    D_dalk  = mp.D_dalk::Float64
    D_dCa   = mp.D_dCa::Float64

    TR_dO2   = mp.TR_dO2::Float64
    TR_dtCO2 = mp.TR_dtCO2::Float64
    TR_dtNO3 = mp.TR_dtNO3::Float64
    TR_dtSO4 = mp.TR_dtSO4::Float64
    TR_dtPO4 = mp.TR_dtPO4::Float64
    TR_dtNH4 = mp.TR_dtNH4::Float64
    TR_dtH2S = mp.TR_dtH2S::Float64
    TR_dFeII = mp.TR_dFeII::Float64
    TR_dMnII = mp.TR_dMnII::Float64
    TR_dCH4  = mp.TR_dCH4::Float64
    TR_dalk  = mp.TR_dalk::Float64
    TR_dCa   = mp.TR_dCa::Float64
    APPW     = mp.APPW::Vector{Float64}

    Ffoc = mp.Ffoc::Float64
    Fsoc = mp.Fsoc::Float64

    z_res = (mp.z_resvec)::Vector{Float64}
    z_rescc = (mp.z_rescc)::Vector{Float64}

    φ0  = mp.phi[1]
    dz0 = z_res[1]   # mp.z_resvec[1]
    dz12 = z_rescc[2] # center-to-center distance for top two cells

    Fcalcite::Float64   = mp.Fcalcite_out
    Faragonite::Float64 = Main.Faragonite
    FMnO2::Float64      = Main.FMnO2
    FFeOH3::Float64     = Main.FFeOH3

    # precompute water-side values as plain Float64
    O2w   = (dO2_w_local * rho_sw)::Float64
    tCO2w = (dtCO2_w * rho_sw)::Float64
    tNO3w = (dtNO3_w * rho_sw)::Float64
    tSO4w = (dtSO4_w * rho_sw)::Float64
    tPO4w = (dtPO4_w * rho_sw)::Float64
    tNH4w = (dtNH4_w * rho_sw)::Float64
    tH2Sw = (dtH2S_w * rho_sw)::Float64
    FeIIw = (dFeII_w * rho_sw)::Float64
    MnIIw = (dMnII_w * rho_sw)::Float64
    CH4w  = (dCH4_w  * rho_sw)::Float64
    alkalw= (dalk_w  * rho_sw)::Float64
    Caw   = (dCa_w   * rho_sw)::Float64

    # in-place reaction scratch (keep it concrete!)
    r = (mp.rates_scratch)::Vector{Float64}
    Nz = size(u, 2)


   # top layer of the sediment
    @inbounds begin
                k = 1

                # clamp just for the reaction inputs
                dO2c   = nz(u[1,k]);  dtCO2c = nz(u[2,k]);  dtNO3c = nz(u[3,k])
                dtSO4c = nz(u[4,k]);  dtPO4c = nz(u[5,k]);  dtNH4c = nz(u[6,k])
                dtH2Sc = nz(u[7,k]);  dFeIIc = nz(u[8,k]);  dMnIIc = nz(u[9,k])
                dCH4c  = nz(u[10,k]); dalkc  = nz(u[11,k]); dCac   = nz(u[12,k])
                pfocc  = nz(u[13,k]); psocc  = nz(u[14,k]); pFeOH3c = nz(u[15,k])
                pMnO2c = nz(u[16,k]); pcalcitec = nz(u[17,k]); paragonitec = nz(u[18,k])
                pFeSc  = nz(u[19,k]); pFeS2c = nz(u[20,k]); pS0c   = nz(u[21,k])
                pFeOH3_PO4c = nz(u[22,k]);

                Hk = clamp(mp.H_cache[k], HMIN, HMAX)   # lagged H

                # compute reaction rates using Hk (no H dynamics)
                rates_wD!(r, mp, k, Hk,
                        dO2c, dtCO2c, dtNO3c, dtSO4c, dtPO4c, dtNH4c, dtH2Sc, dFeIIc, dMnIIc, dCH4c,
                        dalkc, dCac, pfocc, psocc, pFeOH3c, pMnO2c, pcalcitec, paragonitec, pFeSc, pFeS2c, pS0c, pFeOH3_PO4c)

                mp.H_diag[k] = Hk

                # 3) Top boundary residuals (SWI) + irrigation + reactions
                du[1,k]  = diffuse_SWI_2R(u[1,k+1],  u[1,k],  O2w,
                                        D_dO2_tort2[1],  D_dO2,  mp.dbl_dO2,  φ0, dz0, dz12) +
                        advectsolute_SWI_2R(u_bur[k], u[1,k],  O2w,
                                                D_dO2, D_dO2_tort2[1], DFF[k], φ0, mp.dbl_dO2, dz0) +
                        irrigate(u[1,k],  O2w,  alpha[k]) + r[1]

                du[2,k]  = diffuse_SWI_2R(u[2,k+1],  u[2,k],  tCO2w,
                                        D_dtCO2_tort2[1], D_dtCO2, mp.dbl_dtCO2, φ0, dz0, dz12) +
                        advectsolute_SWI_2R(u_bur[k], u[2,k],  tCO2w,
                                                D_dtCO2, D_dtCO2_tort2[1], DFF[k], φ0, mp.dbl_dtCO2, dz0) +
                        irrigate(u[2,k],  tCO2w, alpha[k]) + r[2]

                du[3,k]  = diffuse_SWI_2R(u[3,k+1],  u[3,k],  tNO3w,
                                        D_dtNO3_tort2[1], D_dtNO3, mp.dbl_dtNO3, φ0, dz0, dz12) +
                        advectsolute_SWI_2R(u_bur[k], u[3,k],  tNO3w,
                                                D_dtNO3, D_dtNO3_tort2[1], DFF[k], φ0, mp.dbl_dtNO3, dz0) +
                        irrigate(u[3,k],  tNO3w, alpha[k]) + r[3]

                du[4,k]  = diffuse_SWI_2R(u[4,k+1],  u[4,k],  tSO4w,
                                        D_dtSO4_tort2[1], D_dtSO4, mp.dbl_dtSO4, φ0, dz0, dz12) +
                        advectsolute_SWI_2R(u_bur[k], u[4,k],  tSO4w,
                                                D_dtSO4, D_dtSO4_tort2[1], DFF[k], φ0, mp.dbl_dtSO4, dz0) +
                        irrigate(u[4,k],  tSO4w, alpha[k]) + r[4]

                du[5,k]  = diffuse_SWI_2R(u[5,k+1],  u[5,k],  tPO4w,
                                        D_dtPO4_tort2[1], D_dtPO4, mp.dbl_dtPO4, φ0, dz0, dz12) +
                        advectsolute_SWI_2R(u_bur[k], u[5,k],  tPO4w,
                                                D_dtPO4, D_dtPO4_tort2[1], DFF[k], φ0, mp.dbl_dtPO4, dz0) +
                        irrigate(u[5,k],  tPO4w, alpha[k]) + r[5]

                du[6,k]  = diffuse_SWI_2R(u[6,k+1],  u[6,k],  tNH4w,
                                        D_dtNH4_tort2[1], D_dtNH4, mp.dbl_dtNH4, φ0, dz0, dz12) +
                        advectsolute_SWI_2R(u_bur[k], u[6,k],  tNH4w,
                                                D_dtNH4, D_dtNH4_tort2[1], DFF[k], φ0, mp.dbl_dtNH4, dz0) +
                        irrigate(u[6,k],  tNH4w, alpha[k]) + r[6]

                du[7,k]  = diffuse_SWI_2R(u[7,k+1],  u[7,k],  tH2Sw,
                                        D_dtH2S_tort2[1], D_dtH2S, mp.dbl_dtH2S, φ0, dz0, dz12) +
                        advectsolute_SWI_2R(u_bur[k], u[7,k],  tH2Sw,
                                                D_dtH2S, D_dtH2S_tort2[1], DFF[k], φ0, mp.dbl_dtH2S, dz0) +
                        irrigate(u[7,k],  tH2Sw, alpha[k]) + r[7]

                du[8,k]  = diffuse_SWI_2R(u[8,k+1],  u[8,k],  FeIIw,
                                        D_dFeII_tort2[1], D_dFeII, mp.dbl_dFeII, φ0, dz0, dz12) +
                        advectsolute_SWI_2R(u_bur[k], u[8,k],  FeIIw,
                                                D_dFeII, D_dFeII_tort2[1], DFF[k], φ0, mp.dbl_dFeII, dz0) +
                        irrigate(u[8,k],  FeIIw, alpha[k]) + r[8]

                du[9,k]  = diffuse_SWI_2R(u[9,k+1],  u[9,k],  MnIIw,
                                        D_dMnII_tort2[1], D_dMnII, mp.dbl_dMnII, φ0, dz0, dz12) +
                        advectsolute_SWI_2R(u_bur[k], u[9,k],  MnIIw,
                                                D_dMnII, D_dMnII_tort2[1], DFF[k], φ0, mp.dbl_dMnII, dz0) +
                        irrigate(u[9,k],  MnIIw, alpha[k]) + r[9]

                du[10,k] = diffuse_SWI_2R(u[10,k+1], u[10,k], CH4w,
                                        D_dCH4_tort2[1],  D_dCH4,  mp.dbl_dCH4,  φ0, dz0, dz12) +
                        advectsolute_SWI_2R(u_bur[k], u[10,k], CH4w,
                                                D_dCH4, D_dCH4_tort2[1], DFF[k], φ0, mp.dbl_dCH4, dz0) +
                        irrigate(u[10,k], CH4w, alpha[k]) + r[10]

                du[11,k] = diffuse_SWI_2R(u[11,k+1], u[11,k], alkalw,
                                        D_dalk_tort2[1],  D_dalk,  mp.dbl_dalk,  φ0, dz0, dz12) +
                        advectsolute_SWI_2R(u_bur[k], u[11,k], alkalw,
                                                D_dalk, D_dalk_tort2[1], DFF[k], φ0, mp.dbl_dalk, dz0) +
                        irrigate(u[11,k], alkalw, alpha[k]) + r[11]

                du[12,k] = diffuse_SWI_2R(u[12,k+1], u[12,k], Caw,
                                        D_dCa_tort2[1],   D_dCa,   mp.dbl_dCa,   φ0, dz0, dz12) +
                        advectsolute_SWI_2R(u_bur[k], u[12,k], Caw,
                                                D_dCa, D_dCa_tort2[1], DFF[k], φ0, mp.dbl_dCa, dz0) +
                        irrigate(u[12,k], Caw, alpha[k]) + r[12]

                # Solids at the SWI
                du[13,k] = diffuseSolid_SWI(u[13,k+1], u[13,k], Ffoc, D_bio[k], phiS[k], w[k], z_res[k]) +
                        advectsolid_SWI(u[13,k], u[13,k+1], Ffoc, D_bio[k], APPW[k], sigma1m[k], sigma[k], sigma1p[k], phiS[k], w[k], z_res[k]) +
                        r[13]

                du[14,k] = diffuseSolid_SWI(u[14,k+1], u[14,k], Fsoc, D_bio[k], phiS[k], w[k], z_res[k]) +
                        advectsolid_SWI(u[14,k], u[14,k+1], Fsoc, D_bio[k], APPW[k], sigma1m[k], sigma[k], sigma1p[k], phiS[k], w[k], z_res[k]) +
                        r[14]

                du[15,k] = diffuseSolid_SWI(u[15,k+1], u[15,k], FFeOH3, D_bio[k], phiS[k], w[k], z_res[k]) +
                        advectsolid_SWI(u[15,k], u[15,k+1], FFeOH3, D_bio[k], APPW[k], sigma1m[k], sigma[k], sigma1p[k], phiS[k], w[k], z_res[k]) +
                        r[15]

                du[16,k] = diffuseSolid_SWI(u[16,k+1], u[16,k], FMnO2, D_bio[k], phiS[k], w[k], z_res[k]) +
                        advectsolid_SWI(u[16,k], u[16,k+1], FMnO2, D_bio[k], APPW[k], sigma1m[k], sigma[k], sigma1p[k], phiS[k], w[k], z_res[k]) +
                        r[16]

                du[17,k] = diffuseSolid_SWI(u[17,k+1], u[17,k], Fcalcite, D_bio[k], phiS[k], w[k], z_res[k]) +
                        advectsolid_SWI(u[17,k], u[17,k+1], Fcalcite, D_bio[k], APPW[k], sigma1m[k], sigma[k], sigma1p[k], phiS[k], w[k], z_res[k]) +
                        r[17]

                du[18,k] = diffuseSolid_SWI(u[18,k+1], u[18,k], Faragonite, D_bio[k], phiS[k], w[k], z_res[k]) +
                        advectsolid_SWI(u[18,k], u[18,k+1], Faragonite, D_bio[k], APPW[k], sigma1m[k], sigma[k], sigma1p[k], phiS[k], w[k], z_res[k]) +
                        r[18]

                du[19,k] = diffuseSolid_SWI(u[19,k+1], u[19,k], FFeS, D_bio[k], phiS[k], w[k], z_res[k]) +
                        advectsolid_SWI(u[19,k], u[19,k+1], FFeS, D_bio[k], APPW[k], sigma1m[k], sigma[k], sigma1p[k], phiS[k], w[k], z_res[k]) +
                        r[19]

                du[20,k] = diffuseSolid_SWI(u[20,k+1], u[20,k], FFeS2, D_bio[k], phiS[k], w[k], z_res[k]) +
                        advectsolid_SWI(u[20,k], u[20,k+1], FFeS2, D_bio[k], APPW[k], sigma1m[k], sigma[k], sigma1p[k], phiS[k], w[k], z_res[k]) +
                        r[20]

                du[21,k] = diffuseSolid_SWI(u[21,k+1], u[21,k], FS0, D_bio[k], phiS[k], w[k], z_res[k]) +
                        advectsolid_SWI(u[21,k], u[21,k+1], FS0, D_bio[k], APPW[k], sigma1m[k], sigma[k], sigma1p[k], phiS[k], w[k], z_res[k]) +
                        r[21]

                du[22,k] = diffuseSolid_SWI(u[22,k+1], u[22,k], FFeOH3_PO4, D_bio[k], phiS[k], w[k], z_res[k]) +
                        advectsolid_SWI(u[22,k], u[22,k+1], FFeOH3_PO4, D_bio[k], APPW[k], sigma1m[k], sigma[k], sigma1p[k], phiS[k], w[k], z_res[k]) +
                        r[22]

    end

    @inbounds for k in 2:Nz-1
        
                # clamp just for the reaction inputs
                dO2c   = nz(u[1,k]);  dtCO2c = nz(u[2,k]);  dtNO3c = nz(u[3,k])
                dtSO4c = nz(u[4,k]);  dtPO4c = nz(u[5,k]);  dtNH4c = nz(u[6,k])
                dtH2Sc = nz(u[7,k]);  dFeIIc = nz(u[8,k]);  dMnIIc = nz(u[9,k])
                dCH4c  = nz(u[10,k]); dalkc  = nz(u[11,k]); dCac   = nz(u[12,k])
                pfocc  = nz(u[13,k]); psocc  = nz(u[14,k]); pFeOH3c = nz(u[15,k])
                pMnO2c = nz(u[16,k]); pcalcitec = nz(u[17,k]); paragonitec = nz(u[18,k]) 
                pFeSc  = nz(u[19,k]); pFeS2c = nz(u[20,k]); pS0c   = nz(u[21,k]); pFeOH3_PO4c = nz(u[22,k]);

                Hk = clamp(mp.H_cache[k], HMIN, HMAX)   # lagged H

                # compute reaction rates using Hk (no H dynamics)
                rates_wD!(r, mp, k, Hk,
                        dO2c, dtCO2c, dtNO3c, dtSO4c, dtPO4c, dtNH4c, dtH2Sc, dFeIIc, dMnIIc, dCH4c,
                        dalkc, dCac, pfocc, psocc, pFeOH3c, pMnO2c, pcalcitec, paragonitec, pFeSc, pFeS2c, pS0c, pFeOH3_PO4c)

                mp.H_diag[k] = Hk

                # ---------- 3) assemble du[:,k] ----------
                # SOLUTES (centered diffusion + centered advection + irrigation + reactions)
                du[1, k]  = diffuse(u[1,k-1],  u[1,k],  u[1,k+1],  D_dO2_tort2[k], z_rescc[k])  +
                        advectsolute(u[1,k+1], u[1,k-1], u_bur[k], D_dO2,  DFF[k], z_rescc[k]) +
                        irrigate(u[1,k],  O2w,   alpha[k]) + r[1]

                du[2, k]  = diffuse(u[2,k-1],  u[2,k],  u[2,k+1],  D_dtCO2_tort2[k], z_rescc[k]) +
                        advectsolute(u[2,k+1], u[2,k-1], u_bur[k], D_dtCO2, DFF[k], z_rescc[k]) +
                        irrigate(u[2,k],  tCO2w, alpha[k]) + r[2]

                du[3, k]  = diffuse(u[3,k-1],  u[3,k],  u[3,k+1],  D_dtNO3_tort2[k], z_rescc[k]) +
                        advectsolute(u[3,k+1], u[3,k-1], u_bur[k], D_dtNO3, DFF[k], z_rescc[k]) +
                        irrigate(u[3,k],  tNO3w, alpha[k]) + r[3]

                du[4, k]  = diffuse(u[4,k-1],  u[4,k],  u[4,k+1],  D_dtSO4_tort2[k], z_rescc[k]) +
                        advectsolute(u[4,k+1], u[4,k-1], u_bur[k], D_dtSO4, DFF[k], z_rescc[k]) +
                        irrigate(u[4,k],  tSO4w, alpha[k]) + r[4]

                du[5, k]  = diffuse(u[5,k-1],  u[5,k],  u[5,k+1],  D_dtPO4_tort2[k], z_rescc[k]) +
                        advectsolute(u[5,k+1], u[5,k-1], u_bur[k], D_dtPO4, DFF[k], z_rescc[k]) +
                        irrigate(u[5,k],  tPO4w, alpha[k]) + r[5]

                du[6, k]  = diffuse(u[6,k-1],  u[6,k],  u[6,k+1],  D_dtNH4_tort2[k], z_rescc[k]) +
                        advectsolute(u[6,k+1], u[6,k-1], u_bur[k], D_dtNH4, DFF[k], z_rescc[k]) +
                        irrigate(u[6,k],  tNH4w, alpha[k]) + r[6]

                du[7, k]  = diffuse(u[7,k-1],  u[7,k],  u[7,k+1],  D_dtH2S_tort2[k], z_rescc[k]) +
                        advectsolute(u[7,k+1], u[7,k-1], u_bur[k], D_dtH2S, DFF[k], z_rescc[k]) +
                        irrigate(u[7,k],  tH2Sw, alpha[k]) + r[7]

                du[8, k]  = diffuse(u[8,k-1],  u[8,k],  u[8,k+1],  D_dFeII_tort2[k], z_rescc[k]) +
                        advectsolute(u[8,k+1], u[8,k-1], u_bur[k], D_dFeII, DFF[k], z_rescc[k]) +
                        irrigate(u[8,k],  FeIIw, alpha[k]) + r[8]

                du[9, k]  = diffuse(u[9,k-1],  u[9,k],  u[9,k+1],  D_dMnII_tort2[k], z_rescc[k]) +
                        advectsolute(u[9,k+1], u[9,k-1], u_bur[k], D_dMnII, DFF[k], z_rescc[k]) +
                        irrigate(u[9,k],  MnIIw, alpha[k]) + r[9]

                du[10, k] = diffuse(u[10,k-1], u[10,k], u[10,k+1], D_dCH4_tort2[k], z_rescc[k]) +
                        advectsolute(u[10,k+1], u[10,k-1], u_bur[k], D_dCH4, DFF[k], z_rescc[k]) +
                        irrigate(u[10,k], CH4w,  alpha[k]) + r[10]

                du[11, k] = diffuse(u[11,k-1], u[11,k], u[11,k+1], D_dalk_tort2[k], z_rescc[k]) +
                        advectsolute(u[11,k+1], u[11,k-1], u_bur[k], D_dalk, DFF[k], z_rescc[k]) +
                        irrigate(u[11,k], alkalw, alpha[k]) + r[11]

                du[12, k] = diffuse(u[12,k-1], u[12,k], u[12,k+1], D_dCa_tort2[k], z_rescc[k]) +
                        advectsolute(u[12,k+1], u[12,k-1], u_bur[k], D_dCa,  DFF[k], z_rescc[k]) +
                        irrigate(u[12,k], Caw,    alpha[k]) + r[12]

                # SOLIDS
                du[13, k] = diffuse(u[13,k-1], u[13,k], u[13,k+1], D_bio[k], z_rescc[k]) +
                        advectsolid(u[13,k], u[13,k+1], u[13,k-1],
                                        APPW[k], sigma[k], sigma1p[k], sigma1m[k], z_rescc[k]) +
                        r[13]

                du[14, k] = diffuse(u[14,k-1], u[14,k], u[14,k+1], D_bio[k], z_rescc[k]) +
                        advectsolid(u[14,k], u[14,k+1], u[14,k-1],
                                        APPW[k], sigma[k], sigma1p[k], sigma1m[k], z_rescc[k]) +
                        r[14]

                du[15, k] = diffuse(u[15,k-1], u[15,k], u[15,k+1], D_bio[k], z_rescc[k]) +
                        advectsolid(u[15,k], u[15,k+1], u[15,k-1],
                                        APPW[k], sigma[k], sigma1p[k], sigma1m[k], z_rescc[k]) +
                        r[15]

                du[16, k] = diffuse(u[16,k-1], u[16,k], u[16,k+1], D_bio[k], z_rescc[k]) +
                        advectsolid(u[16,k], u[16,k+1], u[16,k-1],
                                        APPW[k], sigma[k], sigma1p[k], sigma1m[k], z_rescc[k]) +
                        r[16]

                du[17, k] = diffuse(u[17,k-1], u[17,k], u[17,k+1], D_bio[k], z_rescc[k]) +
                        advectsolid(u[17,k], u[17,k+1], u[17,k-1],
                                        APPW[k], sigma[k], sigma1p[k], sigma1m[k], z_rescc[k]) +
                        r[17]

                du[18, k] = diffuse(u[18,k-1], u[18,k], u[18,k+1], D_bio[k], z_rescc[k]) +
                        advectsolid(u[18,k], u[18,k+1], u[18,k-1],
                                        APPW[k], sigma[k], sigma1p[k], sigma1m[k], z_rescc[k]) +
                        r[18]

                du[19, k] = diffuse(u[19,k-1], u[19,k], u[19,k+1], D_bio[k], z_rescc[k]) +
                        advectsolid(u[19,k], u[19,k+1], u[19,k-1],
                                        APPW[k], sigma[k], sigma1p[k], sigma1m[k], z_rescc[k]) +
                        r[19]

                du[20, k] = diffuse(u[20,k-1], u[20,k], u[20,k+1], D_bio[k], z_rescc[k]) +
                        advectsolid(u[20,k], u[20,k+1], u[20,k-1],
                                        APPW[k], sigma[k], sigma1p[k], sigma1m[k], z_rescc[k]) +
                        r[20]

                du[21, k] = diffuse(u[21,k-1], u[21,k], u[21,k+1], D_bio[k], z_rescc[k]) +
                        advectsolid(u[21,k], u[21,k+1], u[21,k-1],
                                        APPW[k], sigma[k], sigma1p[k], sigma1m[k], z_rescc[k]) +
                        r[21]

                du[22, k] = diffuse(u[22,k-1], u[22,k], u[22,k+1], D_bio[k], z_rescc[k]) +
                        advectsolid(u[22,k], u[22,k+1], u[22,k-1],
                                        APPW[k], sigma[k], sigma1p[k], sigma1m[k], z_rescc[k]) +
                        r[22]

                # H (see earlier note re algebraic vs relaxation)
    end

        # --- BOTTOM layer (k = Nz) ---

@inbounds begin
                k = Nz

                # clamp just for the reaction inputs
                dO2c   = nz(u[1,k]);  dtCO2c = nz(u[2,k]);  dtNO3c = nz(u[3,k])
                dtSO4c = nz(u[4,k]);  dtPO4c = nz(u[5,k]);  dtNH4c = nz(u[6,k])
                dtH2Sc = nz(u[7,k]);  dFeIIc = nz(u[8,k]);  dMnIIc = nz(u[9,k])
                dCH4c  = nz(u[10,k]); dalkc  = nz(u[11,k]); dCac   = nz(u[12,k])
                pfocc  = nz(u[13,k]); psocc  = nz(u[14,k]); pFeOH3c = nz(u[15,k])
                pMnO2c = nz(u[16,k]); pcalcitec = nz(u[17,k]); paragonitec = nz(u[18,k])
                pFeSc  = nz(u[19,k]); pFeS2c = nz(u[20,k]); pS0c   = nz(u[21,k]); pFeOH3_PO4c = nz(u[22,k])

                Hk = clamp(mp.H_cache[k], HMIN, HMAX)   # lagged H

                # compute reaction rates using Hk (no H dynamics)
                rates_wD!(r, mp, k, Hk,
                        dO2c, dtCO2c, dtNO3c, dtSO4c, dtPO4c, dtNH4c, dtH2Sc, dFeIIc, dMnIIc, dCH4c,
                        dalkc, dCac, pfocc, psocc, pFeOH3c, pMnO2c, pcalcitec, paragonitec, pFeSc, pFeS2c, pS0c, pFeOH3_PO4c)

                mp.H_diag[k] = Hk


                # ---- SOLUTES (Neumann diffusion + irrigation + reactions) ----
                du[1, k]  = diffuse_BBC(u[1,k],  u[1,k-1],  D_dO2_tort2[k], z_res[k]) + irrigate(u[1,k],  O2w,   alpha[k]) + r[1]
                du[2, k]  = diffuse_BBC(u[2,k],  u[2,k-1],  D_dtCO2_tort2[k], z_res[k]) + irrigate(u[2,k],  tCO2w, alpha[k]) + r[2]
                du[3, k]  = diffuse_BBC(u[3,k],  u[3,k-1],  D_dtNO3_tort2[k], z_res[k]) + irrigate(u[3,k],  tNO3w, alpha[k]) + r[3]
                du[4, k]  = diffuse_BBC(u[4,k],  u[4,k-1],  D_dtSO4_tort2[k], z_res[k]) + irrigate(u[4,k],  tSO4w, alpha[k]) + r[4]
                du[5, k]  = diffuse_BBC(u[5,k],  u[5,k-1],  D_dtPO4_tort2[k], z_res[k]) + irrigate(u[5,k],  tPO4w, alpha[k]) + r[5]
                du[6, k]  = diffuse_BBC(u[6,k],  u[6,k-1],  D_dtNH4_tort2[k], z_res[k]) + irrigate(u[6,k],  tNH4w, alpha[k]) + r[6]
                du[7, k]  = diffuse_BBC(u[7,k],  u[7,k-1],  D_dtH2S_tort2[k], z_res[k]) + irrigate(u[7,k],  tH2Sw, alpha[k]) + r[7]
                du[8, k]  = diffuse_BBC(u[8,k],  u[8,k-1],  D_dFeII_tort2[k], z_res[k]) + irrigate(u[8,k],  FeIIw, alpha[k]) + r[8]
                du[9, k]  = diffuse_BBC(u[9,k],  u[9,k-1],  D_dMnII_tort2[k], z_res[k]) + irrigate(u[9,k],  MnIIw, alpha[k]) + r[9]
                du[10,k]  = diffuse_BBC(u[10,k], u[10,k-1], D_dCH4_tort2[k], z_res[k]) + irrigate(u[10,k], CH4w,  alpha[k]) + r[10]
                du[11,k]  = diffuse_BBC(u[11,k], u[11,k-1], D_dalk_tort2[k], z_res[k]) + irrigate(u[11,k], alkalw, alpha[k]) + r[11]
                du[12,k]  = diffuse_BBC(u[12,k], u[12,k-1], D_dCa_tort2[k], z_res[k]) + irrigate(u[12,k], Caw,    alpha[k]) + r[12]

                # ---- SOLIDS (Neumann diffusion + throughflow advection + reactions) ----
                du[13,k]  = diffuse_BBC(u[13,k], u[13,k-1], D_bio[k], z_res[k]) +
                                advectsolid_BBC(u[13,k], u[13,k-1], APPW[k], sigma[k], z_res[k]) + r[13]
                du[14,k]  = diffuse_BBC(u[14,k], u[14,k-1], D_bio[k], z_res[k]) +
                                advectsolid_BBC(u[14,k], u[14,k-1], APPW[k], sigma[k], z_res[k]) + r[14]
                du[15,k]  = diffuse_BBC(u[15,k], u[15,k-1], D_bio[k], z_res[k]) +
                                advectsolid_BBC(u[15,k], u[15,k-1], APPW[k], sigma[k], z_res[k]) + r[15]
                du[16,k]  = diffuse_BBC(u[16,k], u[16,k-1], D_bio[k], z_res[k]) +
                                advectsolid_BBC(u[16,k], u[16,k-1], APPW[k], sigma[k], z_res[k]) + r[16]
                du[17,k]  = diffuse_BBC(u[17,k], u[17,k-1], D_bio[k], z_res[k]) +
                                advectsolid_BBC(u[17,k], u[17,k-1], APPW[k], sigma[k], z_res[k]) + r[17]
                du[18,k]  = diffuse_BBC(u[18,k], u[18,k-1], D_bio[k], z_res[k]) +
                                advectsolid_BBC(u[18,k], u[18,k-1], APPW[k], sigma[k], z_res[k]) + r[18]
                du[19,k]  = diffuse_BBC(u[19,k], u[19,k-1], D_bio[k], z_res[k]) +
                                advectsolid_BBC(u[19,k], u[19,k-1], APPW[k], sigma[k], z_res[k]) + r[19]
                du[20,k]  = diffuse_BBC(u[20,k], u[20,k-1], D_bio[k], z_res[k]) +
                                advectsolid_BBC(u[20,k], u[20,k-1], APPW[k], sigma[k], z_res[k]) + r[20]
                du[21,k]  = diffuse_BBC(u[21,k], u[21,k-1], D_bio[k], z_res[k]) +
                                advectsolid_BBC(u[21,k], u[21,k-1], APPW[k], sigma[k], z_res[k]) + r[21]
                du[22,k]  = diffuse_BBC(u[22,k], u[22,k-1], D_bio[k], z_res[k]) +
                                advectsolid_BBC(u[22,k], u[22,k-1], APPW[k], sigma[k], z_res[k]) + r[22]


        end
    
end #end function

# %%
prob = ODEProblem(physics_ensamble!, u0, tspan, (
    T = T,
    U = U,
    P = P,
    Fpom = Fpom,
    Fcalcite = Fcalcite,
    model_params = model_params,
    dO2_w = dO2_w
))

# %%
using SparseArrays, LinearSolve, OrdinaryDiffEq

# Block-tridiagonal sparsity pattern (nvar per layer, Nz layers)
function jac_prototype(mp; nvar=22)
    Nz = length(mp.zc)
    N  = nvar*Nz
    rows = Int[]; cols = Int[]
    for k in 1:Nz
        base = (k-1)*nvar
        # dense-ish block on the diagonal (reactions couple vars within a cell)
        for i in 1:nvar, j in 1:nvar
            push!(rows, base+i); push!(cols, base+j)
        end
        # off-diagonal: transport couples same species to k±1 (cheap approximation)
        if k < Nz
            for i in 1:nvar
                push!(rows, base+i);       push!(cols, base+nvar+i)   # k -> k+1
                push!(rows, base+nvar+i);  push!(cols, base+i)         # k+1 -> k
            end
        end
    end
    return sparse(rows, cols, ones(length(rows)), N, N)
end

Jp = jac_prototype(model_params)

f  = ODEFunction(physics_ensamble!;
                 jac_prototype=Jp)   # let the solver keep/factor a matrix of this shape

prob = ODEProblem(f, u0, tspan, (model_params=model_params, dO2_w=dO2_w,))
# sol  = solve(prob,
#             Rosenbrock23(autodiff=false);
#              abstol=1e-6, reltol=1e-4, save_everystep=false, save_on=false, dense=false)

# %% [markdown]
# ### Parallel ensamble run

# %%
using DiffEqCallbacks

# -----------------------------
# SWI flux diagnostics (drop-in)
# -----------------------------

"""
Return SWI diffusive exchange flux for a solute (positive upward, out of sediment),
consistent with your TR_* series conductance formulation.

Units: if TR_* is in m/yr and C in mol/m^3, then J is mol/m^2/yr (same time basis as RHS).
"""
@inline function k_swi_series(Dw::Float64, Deff::Float64, phi0::Float64,
                             dbl::Float64, dz::Float64)
    δ = max(dbl, eps())
    Δ = max(0.5*dz, eps())
    ϕ0 = max(phi0, eps())
    return 1.0 / (δ/Dw + Δ/(ϕ0*Deff))  # m/yr
end

@inline function J_swi_diff(C1::Float64, Cw::Float64,
                              Dw::Float64, Deff::Float64,
                              phi0::Float64, dbl::Float64, dz::Float64)
    k = k_swi_series(Dw, Deff, phi0, dbl, dz)
    return k * (C1 - Cw)  # mol m^-2 yr^-1, positive upward
end

"""
Net irrigation exchange with the overlying water (positive upward, out of sediment),
consistent with `irrigate(C, Cw, alpha) = alpha*(Cw - C)` used in the RHS.

RHS irrigation adds +alpha*(Cw - C) to du, so if Cw > C it is an input to sediments.
The corresponding *export* (upward) is negative. Hence J_irr_out = sum(phi*dz*alpha*(C - Cw)).
"""
function J_irr_net(u_var::AbstractVector{<:Real}, Cw::Float64,
                   phi::AbstractVector{<:Real}, alpha::AbstractVector{<:Real},
                   dz::AbstractVector{<:Real})
    @inbounds begin
        s = 0.0
        for k in eachindex(u_var)
            s += phi[k] * dz[k] * alpha[k] * (u_var[k] - Cw)  # positive upward (out of sediment)
        end
        return s
    end
end

function compute_swi_fluxes(u, mp, dO2_w_local)
    ρ = mp.rho_sw
    phi0 = mp.phi[1]
    dz0 = mp.z_resvec[1]

    # water-side values (mol/m^3)
    O2w   = dO2_w_local * ρ
    tCO2w = Main.dtCO2_w * ρ
    tNO3w = Main.dtNO3_w * ρ
    tSO4w = Main.dtSO4_w * ρ
    tPO4w = Main.dtPO4_w * ρ
    tNH4w = Main.dtNH4_w * ρ
    tH2Sw = Main.dtH2S_w * ρ
    FeIIw = Main.dFeII_w * ρ
    MnIIw = Main.dMnII_w * ρ
    CH4w  = Main.dCH4_w  * ρ
    ALKw  = Main.dalk_w  * ρ
    Caw   = Main.dCa_w   * ρ

    # top-cell values
    O2_1   = u[1,1]
    tCO2_1 = u[2,1]
    tNO3_1 = u[3,1]
    tSO4_1 = u[4,1]
    tPO4_1 = u[5,1]
    tNH4_1 = u[6,1]
    tH2S_1 = u[7,1]
    FeII_1 = u[8,1]
    MnII_1 = u[9,1]
    CH4_1  = u[10,1]
    ALK_1  = u[11,1]
    Ca_1   = u[12,1]

    # SWI diffusive exchange (positive upward / out of sediment)
    Jdiff = (
        O2   = J_swi_diff(O2_1,   O2w,   mp.D_dO2,   mp.D_dO2_tort2[1],   phi0, mp.dbl_dO2,   dz0),
        tCO2 = J_swi_diff(tCO2_1, tCO2w, mp.D_dtCO2, mp.D_dtCO2_tort2[1], phi0, mp.dbl_dtCO2, dz0),
        tNO3 = J_swi_diff(tNO3_1, tNO3w, mp.D_dtNO3, mp.D_dtNO3_tort2[1], phi0, mp.dbl_dtNO3, dz0),
        tSO4 = J_swi_diff(tSO4_1, tSO4w, mp.D_dtSO4, mp.D_dtSO4_tort2[1], phi0, mp.dbl_dtSO4, dz0),
        tPO4 = J_swi_diff(tPO4_1, tPO4w, mp.D_dtPO4, mp.D_dtPO4_tort2[1], phi0, mp.dbl_dtPO4, dz0),
        tNH4 = J_swi_diff(tNH4_1, tNH4w, mp.D_dtNH4, mp.D_dtNH4_tort2[1], phi0, mp.dbl_dtNH4, dz0),
        tH2S = J_swi_diff(tH2S_1, tH2Sw, mp.D_dtH2S, mp.D_dtH2S_tort2[1], phi0, mp.dbl_dtH2S, dz0),
        FeII = J_swi_diff(FeII_1, FeIIw, mp.D_dFeII, mp.D_dFeII_tort2[1], phi0, mp.dbl_dFeII, dz0),
        MnII = J_swi_diff(MnII_1, MnIIw, mp.D_dMnII, mp.D_dMnII_tort2[1], phi0, mp.dbl_dMnII, dz0),
        CH4  = J_swi_diff(CH4_1,  CH4w,  mp.D_dCH4,  mp.D_dCH4_tort2[1],  phi0, mp.dbl_dCH4,  dz0),
        Alk  = J_swi_diff(ALK_1,  ALKw,  mp.D_dalk,  mp.D_dalk_tort2[1],  phi0, mp.dbl_dalk,  dz0),
        Ca   = J_swi_diff(Ca_1,   Caw,   mp.D_dCa,   mp.D_dCa_tort2[1],   phi0, mp.dbl_dCa,   dz0),
    )

    # net irrigation exchange (depth-integrated), same sign convention
    ϕ  = mp.phi
    α  = mp.alpha
    dz = mp.z_resvec

    Jirr = (
        O2   = J_irr_net(view(u,1,:),   O2w,   ϕ, α, dz),
        tCO2 = J_irr_net(view(u,2,:),   tCO2w, ϕ, α, dz),
        tNO3 = J_irr_net(view(u,3,:),   tNO3w, ϕ, α, dz),
        tSO4 = J_irr_net(view(u,4,:),   tSO4w, ϕ, α, dz),
        tPO4 = J_irr_net(view(u,5,:),   tPO4w, ϕ, α, dz),
        tNH4 = J_irr_net(view(u,6,:),   tNH4w, ϕ, α, dz),
        tH2S = J_irr_net(view(u,7,:),   tH2Sw, ϕ, α, dz),
        FeII = J_irr_net(view(u,8,:),   FeIIw, ϕ, α, dz),
        MnII = J_irr_net(view(u,9,:),   MnIIw, ϕ, α, dz),
        CH4  = J_irr_net(view(u,10,:),  CH4w,  ϕ, α, dz),
        Alk  = J_irr_net(view(u,11,:),  ALKw,  ϕ, α, dz),
        Ca   = J_irr_net(view(u,12,:),  Caw,   ϕ, α, dz),
    )

    # total net fluxes as local scalars first
    J_O2   = Jdiff.O2   + Jirr.O2
    J_tCO2 = Jdiff.tCO2 + Jirr.tCO2
    J_tNO3 = Jdiff.tNO3 + Jirr.tNO3
    J_tSO4 = Jdiff.tSO4 + Jirr.tSO4
    J_tPO4 = Jdiff.tPO4 + Jirr.tPO4
    J_tNH4 = Jdiff.tNH4 + Jirr.tNH4
    J_tH2S = Jdiff.tH2S + Jirr.tH2S
    J_FeII = Jdiff.FeII + Jirr.FeII
    J_MnII = Jdiff.MnII + Jirr.MnII
    J_CH4  = Jdiff.CH4  + Jirr.CH4
    J_Alk  = Jdiff.Alk  + Jirr.Alk
    J_Ca   = Jdiff.Ca   + Jirr.Ca

    Jnet = (
        O2   = J_O2,
        tCO2 = J_tCO2,
        tNO3 = J_tNO3,
        tSO4 = J_tSO4,
        tPO4 = J_tPO4,
        tNH4 = J_tNH4,
        tH2S = J_tH2S,
        FeII = J_FeII,
        MnII = J_MnII,
        CH4  = J_CH4,
        Alk  = J_Alk,
        Ca   = J_Ca,
        Alk_eff_H2Scorr   = J_Alk - 2.0*J_tH2S,
        Alk_eff_totalcorr = J_Alk - 2.0*J_tH2S - 2.0*J_FeII - 2.0*J_MnII
    )

    Nz_ = size(u, 2)


    OmegaCa_profile = Vector{Float64}(undef, Nz_)
    @inbounds for k in 1:Nz_
        dCa_k   = max(u[12, k], 0.0)
        dtCO2_k = max(u[2,  k], 0.0)
        Hk      = clamp(mp.H_cache[k], HMIN, HMAX)
        denom_k = mp.K1*mp.K2 + mp.K1*Hk + Hk*Hk
        dCO3_k  = (dtCO2_k * mp.K1 * mp.K2) / denom_k
        OmegaCa_profile[k] = dCa_k * dCO3_k / mp.KCa
    end

    H_profile = copy(mp.H_cache)

    return (Jdiff=Jdiff, Jirr=Jirr, Jnet=Jnet, OmegaCa=OmegaCa_profile, H = H_profile)
end

"""
Save-ready reaction snapshot for one model state.

Returns a nested NamedTuple with:
- `reactions`: individual reaction rates by depth
- `rates`: net source/sink terms for each model variable by depth
- `aggregates`: aggregated degradation rates by depth
- `carbonate`: `H`, `dCO3`, and `OmegaCa` on the same save grid
"""
function compute_reaction_rate_snapshot(u, mp)
    Nz_ = size(u, 2)

    reaction_names = (
        :Rfast_dO2, :Rslow_dO2, :Rfast_dtNO3, :Rslow_dtNO3,
        :Rfast_pMnO2, :Rslow_pMnO2, :Rfast_pFeOH3, :Rslow_pFeOH3,
        :Rfast_dtSO4, :Rslow_dtSO4, :Rfast_dCH4, :Rslow_dCH4,
        :Rfast_total, :Rslow_total,
        :R_dMnII, :R_dFeII, :R_dNH3, :R_dH2S, :R_CH4_O2redox, :R_CH4_SO4redox,
        :R_FEOH3_PO4_adsorp, :R_Fe_MnO2_red, :R_H2S_FeOOH_PO4_red,
        :R_H2S_FeOOH_red, :R_H2S_MnO2_red, :R_FeS_H2S_Fe,
        :R_FeS2_FeS_S0, :R_FeS2_SO4_H2S_FeS, :R_FeS_ox, :R_FeS2_O2,
        :R_S0_H20, :R_MnO2a_MnO2b, :R_FeOOHa_FeOOHb,
        :Rdiss_calcite, :Rdiss_aragonite, :Rprec_calcite, :Rprec_aragonite,
    )

    rate_names = (
        :dO2, :dtCO2, :dtNO3, :dtSO4, :dtPO4, :dtNH4, :dtH2S,
        :dFeII, :dMnII, :dCH4, :dalk, :dCa,
        :pfoc, :psoc, :pFeOH3, :pMnO2, :pFeS, :pFeS2,
        :pFeOH3_PO4, :pS0, :pcalcite, :paragonite,
    )

    aggregate_names = (
        :Rdeg_dO2, :Rdeg_dtNO3, :Rdeg_dtSO4,
        :Rdeg_pFeOH3, :Rdeg_pMnO2, :Rdeg_dCH4, :Rdeg_total,
    )

    reaction_profiles = NamedTuple{reaction_names}(Tuple(Vector{Float64}(undef, Nz_) for _ in reaction_names))
    rate_profiles = NamedTuple{rate_names}(Tuple(Vector{Float64}(undef, Nz_) for _ in rate_names))
    aggregate_profiles = NamedTuple{aggregate_names}(Tuple(Vector{Float64}(undef, Nz_) for _ in aggregate_names))

    H_profile = copy(mp.H_cache)
    dCO3_profile = Vector{Float64}(undef, Nz_)
    OmegaCa_profile = Vector{Float64}(undef, Nz_)

    @inbounds for k in 1:Nz_
        dO2   = nz(u[1,  k])
        dtCO2 = nz(u[2,  k])
        dtNO3 = nz(u[3,  k])
        dtSO4 = nz(u[4,  k])
        dtPO4 = nz(u[5,  k])
        dtNH4 = nz(u[6,  k])
        dtH2S = nz(u[7,  k])
        dFeII = nz(u[8,  k])
        dMnII = nz(u[9,  k])
        dCH4  = nz(u[10, k])
        dCa   = nz(u[12, k])
        pfoc  = nz(u[13, k])
        psoc  = nz(u[14, k])
        pFeOH3 = nz(u[15, k])
        pMnO2  = nz(u[16, k])
        pcalcite = nz(u[17, k])
        paragonite = nz(u[18, k])
        pFeS = nz(u[19, k])
        pFeS2 = nz(u[20, k])
        pS0 = nz(u[21, k])
        pFeOH3_PO4 = nz(u[22, k])

        Hk = clamp(H_profile[k], HMIN, HMAX)
        denom_k = mp.K1 * mp.K2 + mp.K1 * Hk + Hk * Hk
        dCO3_k = (dtCO2 * mp.K1 * mp.K2) / denom_k
        dCO3_profile[k] = dCO3_k
        OmegaCa_profile[k] = dCa * dCO3_k / mp.KCa

        kfast_k = (mp.kfast[k])::Float64
        kslow_k = (mp.kslow[k])::Float64
        phiS_phi_k = (mp.phiS_phi[k])::Float64

        rxn_vals = React.getreactions(
            dO2, dtNO3, pMnO2, pFeOH3, dtSO4, dtNH4, dtH2S, dFeII, dMnII, dCH4,
            dtPO4, pFeOH3_PO4, pFeS, pS0, pFeS2,
            pfoc * kfast_k, psoc * kslow_k,
            pcalcite, paragonite, dCa, dCO3_k, mp.KCa, mp.KAr,
            T, calcite_diss_scheme, aragonite_diss_scheme, Q10_secondary, Tref,
        )

        for (name, val) in zip(reaction_names, rxn_vals)
            getproperty(reaction_profiles, name)[k] = val
        end

        net_vals = React.reactions2rates(rxn_vals..., phiS_phi_k, mp.RC, mp.RN, mp.RP)
        for (name, val) in zip(rate_names, net_vals[1:length(rate_names)])
            getproperty(rate_profiles, name)[k] = val
        end
        for (name, val) in zip(aggregate_names, net_vals[(length(rate_names) + 1):end])
            getproperty(aggregate_profiles, name)[k] = val
        end
    end

    return (
        reactions = reaction_profiles,
        rates = rate_profiles,
        aggregates = aggregate_profiles,
        carbonate = (H = H_profile, dCO3 = dCO3_profile, OmegaCa = OmegaCa_profile),
    )
end

# -----------------------------
# SavingCallback setup
# -----------------------------

# Choose when to save diagnostics (in model time units; your notebook uses years)
# Example: save every 0.1 yr. Adjust as needed.
flux_saveat = 0.0:0.001:(Main.tspan[2])

# Optional detailed reaction-rate saving.
# Set `save_reaction_rates = true` to store a full depth profile of all reactions.
# Warning: this can generate very large outputs when `reaction_rates_save_everystep = true`.
save_reaction_rates = true
reaction_rates_save_everystep = true
reaction_rate_saveat = flux_saveat

flux_saved = SavedValues(Float64, Any)

flux_cb = SavingCallback(
    (u, t, integrator) -> compute_swi_fluxes(u, integrator.p.model_params, dO2_w),
    flux_saved;
    saveat = flux_saveat
)

reaction_saved_single = save_reaction_rates ? SavedValues(Float64, Any) : nothing
reaction_cb = if save_reaction_rates
    if reaction_rates_save_everystep
        SavingCallback(
            (u, t, integrator) -> compute_reaction_rate_snapshot(u, integrator.p.model_params),
            reaction_saved_single;
            save_everystep = true,
            save_start = true,
            save_end = true,
        )
    else
        SavingCallback(
            (u, t, integrator) -> compute_reaction_rate_snapshot(u, integrator.p.model_params),
            reaction_saved_single;
            saveat = reaction_rate_saveat,
            save_start = true,
            save_end = true,
        )
    end
else
    nothing
end

function refresh_H_cache_each_step!(u, t, integrator)
    mp = integrator.p.model_params
    Nz = size(u, 2)

    @inbounds for k in 1:Nz
        dtCO2c = nz(u[2,k])
        dtSO4c = nz(u[4,k])
        dtPO4c = nz(u[5,k])
        dtNH4c = nz(u[6,k])
        dtH2Sc = nz(u[7,k])
        dalkc  = nz(u[11,k])

        H0k = clamp(mp.H_cache[k], HMIN, HMAX)

        Hk = H_from_alk_m3(
            dalkc, dtCO2c, dtNH4c, dtPO4c, dtSO4c, dtH2Sc, mp;
            H0 = H0k, iters = 1, damp = 0.4
        )

        Hk = clamp(Hk, HMIN, HMAX)
        mp.H_cache[k] = Hk
        mp.H_diag[k]  = Hk
    end
    return nothing
end

cb_refreshH_i = FunctionCallingCallback(refresh_H_cache_each_step!; func_everystep=true)


function commit_H_cache!(u, t, integrator)
    mp = integrator.p.model_params
    mp.H_cache .= mp.H_diag
    return nothing
end

const cb_commitH = FunctionCallingCallback(commit_H_cache!;
    func_everystep = true,   # runs after each accepted step
)

# %%
using LinearAlgebra, SparseArrays, SparseDiffTools
using LinearSolve, OrdinaryDiffEq, SciMLBase

# -----------------------------
# 0) Threading / BLAS
# -----------------------------
BLAS.set_num_threads(1)  # avoid over-subscription with EnsembleThreads()

# -----------------------------
# 1) Jacobian sparsity helper
# -----------------------------
function jac_prototype(mp; nvar::Int=22)
    Nz = length(mp.zc)
    N = nvar * Nz
    rows = Int[]
    cols = Int[]
    @inbounds for k in 1:Nz
        base = (k - 1) * nvar
        for i in 1:nvar, j in 1:nvar
            push!(rows, base + i)
            push!(cols, base + j)
        end
        if k < Nz
            for i in 1:nvar
                push!(rows, base + i)
                push!(cols, base + nvar + i)
                push!(rows, base + nvar + i)
                push!(cols, base + i)
            end
        end
    end
    return sparse(rows, cols, ones(length(rows)), N, N)
end

# Cache: (nvar, Nz) -> (Jp, colors)
if !isdefined(@__MODULE__, :J_CACHE)
    const J_CACHE = Dict{Tuple{Int,Int},Tuple{SparseMatrixCSC{Float64,Int},Vector{Int}}}()
else
    empty!(J_CACHE)  # optional: wipe previous entries when you rerun the cell
end


function cached_jac(mp; nvar::Int=22)
    key = (nvar, length(mp.zc))
    get!(J_CACHE, key) do
        Jp = jac_prototype(mp; nvar=nvar)
        colors = SparseDiffTools.matrix_colors(Jp)
        (Jp, colors)
    end
end

# -----------------------------
# 2) Build ICs exactly like your single run (mol/kg * rho_sw)
# -----------------------------
function make_u0_from_IC(mp, dO2_w_local)
    Nz = length(mp.zc)
    ρ = mp.rho_sw
    u0 = zeros(22, Nz)
    u0[1, :] .= dO2_w_local * ρ
    u0[2, :] .= dtCO2_w * ρ
    u0[3, :] .= dtNO3_w * ρ
    u0[4, :] .= dtSO4_w * ρ
    u0[5, :] .= dtPO4_w * ρ
    u0[6, :] .= dtNH4_w * ρ
    u0[7, :] .= dtH2S_w * ρ
    u0[8, :] .= dFeII_w * ρ
    u0[9, :] .= dMnII_w * ρ
    u0[10, :] .= dCH4_w * ρ
    u0[11, :] .= dalk_w * ρ
    dCa_w = 0.02128 / 40.087 * S / 1.80655
    u0[12, :] .= dCa_w * ρ
    u0[13, :] .= pfoc_i
    u0[14, :] .= psoc_i
    u0[15, :] .= pFeOH3_i
    u0[16, :] .= pMnO2_i
    u0[17, :] .= pcalcite_i
    u0[18, :] .= paragonite_i
    u0[19, :] .= pFeS_i
    u0[20, :] .= pFeS2_i
    u0[21, :] .= pS0_i
    u0[22, :] .= pFeOH3_PO4_i
    return u0
end

# -----------------------------
# 3) Build parameter sets
# -----------------------------
# Example: one trajectory; extend to N as needed

# trajectories = 1  # or set explicitly

# # Option A: hard-code levels here
# T_levels        = [10.0, 20.0, 30.0]
# Fpom_levels     = [0.01, 0.1, 1.0]
# U_levels        = [0.01, 0.02, 0.05]          # keep fixed unless you want to vary
# P_levels        = [5, 10, 15]          # keep fixed unless you want to vary
# Fcalcite_levels = [1]   # keep fixed unless you want to vary

# new_T = [11.7 + rand() * (12.8 - 11.7) for i in 1:trajectories]
# new_U = [0.02 + rand() * (0.08 - 0.02) for i in 1:trajectories]
# new_omegaCa = [2.5 + rand() * (5.0 - 2.5) for i in 1:trajectories]
# new_P = [10.0 * 10.325 + rand() * (100.0 * 10.325 - 10.0 * 10.325) for i in 1:trajectories]
# new_Fpom = [26.0 + rand() * (55.0 - 26.0) for i in 1:trajectories]
# new_Fcalcite = [0.2 + rand() * (0.8 - 0.2) for i in 1:trajectories]

T_levels        = @isdefined(factorial_T_levels)        ? factorial_T_levels        : [T]
Fpom_levels     = @isdefined(factorial_Fpom_levels)     ? factorial_Fpom_levels     : [Fpom]
U_levels        = @isdefined(factorial_U_levels)        ? factorial_U_levels        : [U]
P_levels        = @isdefined(factorial_P_levels)        ? factorial_P_levels        : [P]
Fcalcite_levels = @isdefined(factorial_Fcalcite_levels) ? factorial_Fcalcite_levels : [Fcalcite]
O_levels        = @isdefined(factorial_O_levels) ? factorial_O_levels : [dO2_w]

grid = collect(Iterators.product(T_levels, U_levels, P_levels, Fpom_levels, Fcalcite_levels, O_levels))
trajectories = length(grid)

# Keep your existing downstream structure
new_T        = [g[1] for g in grid]
new_U        = [g[2] for g in grid]
new_P        = [g[3] for g in grid]
new_Fpom     = [g[4] for g in grid]
new_Fcalcite = [g[5] for g in grid]
new_O = [g[6] for g in grid]

# new_T = 11.7 
# new_U = 0.02 
# new_omegaCa = 2.5
# new_P = 10.0
# new_Fpom = 26.0 
# new_Fcalcite = 0.2 

# Fix a concrete params type
proto = calculate_constants(new_T[1], new_U[1], new_P[1], new_Fpom[1], 
    new_Fcalcite[1], Fpom_s, Fpom_f, kfast, kslow, dO2_w
)
const ParamsT = typeof(proto)

# --- keep simple per-trajectory metadata (like before) ---
const ParamMetaT = NamedTuple{(:T, :U, :P, :Fpom, :Fcalcite, :O),NTuple{6,Float64}}
param_meta_list = Vector{ParamMetaT}(undef, trajectories)
for i in 1:trajectories
    param_meta_list[i] = (
        T=new_T[i],
        U=new_U[i],
        P=new_P[i],
        Fpom=new_Fpom[i],
        Fcalcite=new_Fcalcite[i],
        O=new_O[i],
    )
end

model_params_list = Vector{ParamsT}(undef, trajectories)
for i in 1:trajectories
    model_params_list[i] = calculate_constants(
        new_T[i], new_U[i], new_P[i],
        new_Fpom[i], new_Fcalcite[i],
        0.2, 0.8, 18.0, 0.05, new_O[i]
    )
end

# Prebuild u0 and p for all trajectories (no allocs in prob_func)
u0_list = [make_u0_from_IC(model_params_list[i], new_O[i]) for i in 1:trajectories]
p_list = [(model_params=model_params_list[i], dO2_w=new_O[i]) for i in 1:trajectories]

# Pre-warm/calc J-cache for all unique Nz (optional but nice)
for mp in model_params_list
    cached_jac(mp)  # fills J_CACHE
end

# -----------------------------
# 4) Base problem (reused when Nz matches)
# -----------------------------
mp0 = model_params_list[1]
Jp0, colors0 = cached_jac(mp0)
f0 = ODEFunction(physics_ensamble!; jac_prototype=Jp0, colorvec=colors0)

prob_base = ODEProblem(f0, u0_list[1], tspan, p_list[1])

# -----------------------------
# 5) prob_func
# -----------------------------
# Save times (choose what you want; years in your setup)
flux_saveat = 0.0:0.001:tspan[2]
# flux_saveat = vcat(
#     0.0:0.01:1.0, 
#     1.1:0.1:8.0, 
#     9.0:1.0:Main.tspan[2]
#     )

# One SavedValues per trajectory
flux_saved = [SavedValues(Float64, Any) for _ in 1:trajectories]
reaction_saved = save_reaction_rates ? [SavedValues(Float64, Any) for _ in 1:trajectories] : nothing

prob_func = function (prob, i, repeat)
    mp   = model_params_list[i]
    u0_i = u0_list[i]
    p_i  = p_list[i]

    fcb = SavingCallback(
        (u, t, integrator) -> compute_swi_fluxes(u, integrator.p.model_params, integrator.p.dO2_w),
        flux_saved[i];
        save_everystep = true,
        save_start = true,
        save_end = true,
        # saveat = flux_saveat
    )

    cb_refreshH_i = FunctionCallingCallback(refresh_H_cache_each_step!; func_everystep=true)
    callbacks = Any[cb_refreshH_i, fcb]

    if save_reaction_rates
        rcb = if reaction_rates_save_everystep
            SavingCallback(
                (u, t, integrator) -> compute_reaction_rate_snapshot(u, integrator.p.model_params),
                reaction_saved[i];
                save_everystep = true,
                save_start = true,
                save_end = true,
            )
        else
            SavingCallback(
                (u, t, integrator) -> compute_reaction_rate_snapshot(u, integrator.p.model_params),
                reaction_saved[i];
                saveat = reaction_rate_saveat,
                save_start = true,
                save_end = true,
            )
        end
        push!(callbacks, rcb)
    end

    cb_i = CallbackSet(callbacks...)

    if size(prob.u0, 2) == length(mp.zc)
        return remake(prob; u0=u0_i, p=p_i, callback=cb_i)
    else
        Jp_i, colors_i = cached_jac(mp)
        f_i = ODEFunction(physics_ensamble!; jac_prototype=Jp_i, colorvec=colors_i)
        return ODEProblem(f_i, u0_i, prob.tspan, p_i; callback=cb_i)
    end
end

# -----------------------------
# 5b) Super-light per-trajectory progress (thread-safe)
# -----------------------------
using Base.Threads: Atomic, atomic_add!

ENSEMBLE_BAR_LOCK = ReentrantLock()
ENSEMBLE_DONE = Atomic{Int}(0)

# Format seconds as HH:MM:SS
fmt_hms(sec) = begin
    s = Int(round(max(sec, 0.0)))
    h = s ÷ 3600
    m = (s % 3600) ÷ 60
    ss = s % 60
    lpad(h, 2, '0') * ":" * lpad(m, 2, '0') * ":" * lpad(ss, 2, '0')
end

function make_output_with_progress(total::Int)
    ENSEMBLE_DONE[] = 0
    t0 = time()
    step = max(1, Int(ceil(total / 20)))  # ~5% steps

    return function (sol, i)
        done = atomic_add!(ENSEMBLE_DONE, 1) + 1  # atomic increment

        if done == 1 || done % step == 0 || done == total
            width = 30
            filled = Int(clamp(round(width * done / total), 0, width))
            bar = "[" * repeat("█", filled) * repeat(" ", width - filled) * "]"
            pct = 100 * done / total
            elapsed = time() - t0
            rate = done / max(elapsed, 1e-9)
            eta = (total - done) / max(rate, 1e-9)

            lock(ENSEMBLE_BAR_LOCK) do
                print("\r", bar, " ",
                    lpad(done, ndigits(total)), "/", total,
                    "  (", round(pct; digits=1), "%)  ",
                    "Elapsed: ", fmt_hms(elapsed), "  ",
                    "ETA: ", fmt_hms(eta))
                flush(stdout)
                if done == total
                    println()
                end
            end
        end

        return (sol, false)   # keep only final state
    end
end

outf = make_output_with_progress(trajectories)

ens = EnsembleProblem(
    prob_base;
    prob_func=prob_func,
    output_func=outf,
    safetycopy=true
)

# -----------------------------
# 6) Algorithm + warm-up (unchanged)
# -----------------------------
linsolve_alg = try
    KLUFactorization()
catch
    UMFPACKFactorization()
end

alg = FBDF(autodiff=false, linsolve=linsolve_alg)

_ = solve(remake(prob_base); alg,
    abstol=abstol_v, reltol=reltol_v,
    save_everystep=false, saveat=flux_saveat, save_on=false, dense=false)


# -----------------------------
# 7) Run ensemble (progress prints here)
# -----------------------------

sols = solve(
    ens, alg, EnsembleThreads();
    trajectories=trajectories,
    abstol=abstol_v / 5, reltol=reltol_v / 5,
    save_on=true, save_everystep=false, saveat=flux_saveat, save_start=true, save_end=true, dense=false, maxiters=5_000_000, dtmin=1e-12
)

# %%
size(sols[1])
sols[1].t

# %%
using DataFrames

const FLUX_CONV = 1000.0 / 365.25   # mol m^-2 yr^-1 -> mmol m^-2 d^-1

"""
Extract final flux diagnostics from one saved flux object.

Assumes:
- flux_saved[i].saveval is a vector of saved diagnostic snapshots
- each snapshot has fields Jnet, Jdiff, Jirr
- each of those is a NamedTuple

Returns a NamedTuple with flat column names ready for DataFrame use.
"""
function extract_final_fluxes(saved_run; conv=FLUX_CONV)
    # final saved diagnostic snapshot for this run
    v = saved_run.saveval[end]

    out = NamedTuple()

    # helper to append fields from one NamedTuple with a prefix
    function add_prefixed(nt::NamedTuple, prefix::String)
        pairs_vec = Pair{Symbol,Float64}[]
        for name in keys(nt)
            col = Symbol(prefix * "_" * String(name))
            val = getproperty(nt, name) * conv
            push!(pairs_vec, col => val)
        end
        return (; pairs_vec...)
    end

    nt_net  = add_prefixed(v.Jnet,  "Jnet")
    nt_diff = add_prefixed(v.Jdiff, "Jdiff")
    nt_irr  = add_prefixed(v.Jirr,  "Jirr")

    return merge(nt_net, nt_diff, nt_irr)
end

# %%
flux_rows = [extract_final_fluxes(flux_saved[i]) for i in eachindex(flux_saved)]
df_fluxes = DataFrame(flux_rows)

# %%
# extract parameter metadata
T_vals        = [param_meta_list[i].T        for i in 1:length(param_meta_list)]
U_vals        = [param_meta_list[i].U        for i in 1:length(param_meta_list)]
P_vals        = [param_meta_list[i].P        for i in 1:length(param_meta_list)]
Fpom_vals     = [param_meta_list[i].Fpom     for i in 1:length(param_meta_list)]
Fcalcite_vals = [param_meta_list[i].Fcalcite for i in 1:length(param_meta_list)]
O_vals        = [param_meta_list[i].O        for i in 1:length(param_meta_list)]
omegaCa_vals  = [model_params_list[i].omegaCa for i in 1:length(model_params_list)]

df_fluxes[!, "T"] = T_vals
df_fluxes[!, "U"] = U_vals
df_fluxes[!, "P"] = P_vals
df_fluxes[!, "Fpom"] = Fpom_vals
df_fluxes[!, "Fcalcite"] = Fcalcite_vals
df_fluxes[!, "O"] = O_vals
df_fluxes[!, "omegaCa"] = omegaCa_vals;


# %%
flux_saved[1].t          # times
flux_saved[1].saveval    # saved values (NamedTuples)
length(flux_saved[1].t)

# %%
times = flux_saved[1].t
JAlk_diff = [v.Jdiff.Alk for v in flux_saved[1].saveval]
JAlk_irr = [v.Jirr.Alk for v in flux_saved[1].saveval]
JAlk_net = [v.Jnet.Alk for v in flux_saved[1].saveval]

# %%
size(sols[1].u[1])

# Variable names in the same order used when populating u0[1:22, :].
var_names = [
    "dO2", "dtCO2", "dtNO3", "dtSO4", "dtPO4", "dtNH4", "dtH2S", "dFeII", "dMnII", "dCH4", "dalk", "dCa",
    "pfoc", "psoc", "pFeOH3", "pMnO2", "pcalcite", "paragonite", "pFeS", "pFeS2", "pS0", "pFeOH3_PO4"
 ]

# # Build one heatmap per variable index x = 1:22.
# heatmaps = Vector{Any}(undef, 22)
# for x in 1:22
#     hmap_data = [u[x, :] for u in sols[1].u]
#     hmap_matrix = vcat([d' for d in hmap_data]...)  # rows = time, columns = depth

#     hm = heatmap(
#         sols[1].t, model_params.zc, hmap_matrix',
#         xlabel="Time (y)",
#         ylabel="Depth (m)",
#         title="$(var_names[x]) Profile Heatmap",
#         yflip=true
#     )

#     heatmaps[x] = hm
#     display(hm)
# end

# heatmaps

# %%
using Dates

stamp = Dates.format(now(), "yyyymmdd_HHMMSS")   # e.g. 20260326_142530
fname = "sols_all_$(stamp).mat"

# -- pack fluxes: each species -> (ntimes_flux x ntrajectories) matrix ---
flux_t = flux_saved[1].t   # same time grid for all trajectories
flux_nmax = maximum(flux_lens)
for i in 1:trajectories
    if flux_lens[i] < flux_nmax
        @warn "Trajectory $i saved only $(flux_lens[i])/$flux_nmax flux time points (aborted early?)"
    end
end
flux_t = flux_saved[argmax(flux_lens)].t

# _pack_scalar(field, subfield) = hcat([
#     [getproperty(v[field], subfield) for v in flux_saved[i].saveval]
#     for i in 1:trajectories
# ]...)   # -> (ntimes_flux, ntrajectories)

function _pack_scalar(field, subfield)
    cols = map(1:trajectories) do i
        vals = [getproperty(getproperty(v, field), subfield) for v in flux_saved[i].saveval]
        n = length(vals)
        n < flux_nmax && append!(vals, fill(NaN, flux_nmax - n))
        vals
    end
    return hcat(cols...)
end

jnet_species  = keys(flux_saved[1].saveval[1].Jnet)
jdiff_species = keys(flux_saved[1].saveval[1].Jdiff)
jirr_species  = keys(flux_saved[1].saveval[1].Jirr)

Jnet_mat  = Dict("Jnet_$(String(sp))"  => _pack_scalar(:Jnet,  sp) for sp in jnet_species)
Jdiff_mat = Dict("Jdiff_$(String(sp))" => _pack_scalar(:Jdiff, sp) for sp in jdiff_species)
Jirr_mat  = Dict("Jirr_$(String(sp))"  => _pack_scalar(:Jirr,  sp) for sp in jirr_species)

# --- OmegaCa:  (Nz x ntimes_flux x ntrajectories) ---
OmegaCa_arr = cat([
    hcat([v.OmegaCa for v in flux_saved[i].saveval]...)  # (Nz, ntimes_flux)
    for i in 1:trajectories
]..., dims=3)

H_arr = cat([
    hcat([v.H for v in flux_saved[i].saveval]...)
    for i in 1:trajectories
]..., dims=3)

reaction_export = Dict{String,Any}()
if save_reaction_rates && reaction_saved !== nothing && !isempty(reaction_saved[1].saveval)
    reaction_t = reaction_saved[1].t

    _pack_reaction(group, field) = cat([
        hcat([getproperty(getproperty(v, group), field) for v in reaction_saved[i].saveval]...)
        for i in 1:trajectories
    ]..., dims=3)  # (Nz, ntimes_reaction, ntrajectories)

    _pack_carbonate(field) = cat([
        hcat([getproperty(v.carbonate, field) for v in reaction_saved[i].saveval]...)
        for i in 1:trajectories
    ]..., dims=3)

    reaction_names = keys(reaction_saved[1].saveval[1].reactions)
    rate_names = keys(reaction_saved[1].saveval[1].rates)
    aggregate_names = keys(reaction_saved[1].saveval[1].aggregates)

    reaction_export = merge(
        Dict(
            "reaction_t" => reaction_t,
            "reaction_H" => _pack_carbonate(:H),
            "reaction_dCO3" => _pack_carbonate(:dCO3),
            "reaction_OmegaCa" => _pack_carbonate(:OmegaCa),
        ),
        Dict("rxn_$(String(name))" => _pack_reaction(:reactions, name) for name in reaction_names),
        Dict("rate_$(String(name))" => _pack_reaction(:rates, name) for name in rate_names),
        Dict("aggregate_$(String(name))" => _pack_reaction(:aggregates, name) for name in aggregate_names),
    )
end

matwrite(fname, merge(
    Dict(
        "t"             => sols[1].t,
        "u"             => cat([cat(sols[i].u..., dims=3) for i in 1:trajectories]..., dims=4),
        "flux_t"        => flux_t,
        "OmegaCa"       => OmegaCa_arr,
        "H"             => H_arr,
        "T_vals"        => new_T,
        "U_vals"        => new_U,
        "Fpom_vals"     => new_Fpom,
        "Fcalcite_vals" => new_Fcalcite,
        "P_vals"        => new_P,
        "O_vals"        => new_O,
        "zc"            => model_params_list[1].zc,
        "S"             => S
    ),
    Jnet_mat,
    Jdiff_mat,
    Jirr_mat,
    reaction_export,
))


println("[sweep] Done. Output written to ", fname)
flush(stdout)
