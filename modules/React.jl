module React

# =========================
# Organic-matter pathways (unchanged)
# =========================

# Monod scheme constants in mol/m^3
# For oxygen and nitrate from Soetaert et al. 1996 (GCA)
const KM_dO2    = 0.003
const KMi_dO2   = 0.01
const KM_dtNO3  = 0.03
const KMi_dtNO3 = 0.005
# For others from van Cappellen and Wang (1996)
const KM_pMnO2   = 42.4
const KMi_pMnO2  = 42.4
const KM_pFeOH3  = 265.0
const KMi_pFeOH3 = 265.0
const KM_dtSO4   = 1.6
const KMi_dtSO4  = 1.6

# Calcite dissolution constants
const ETCH_PIT_RATE_CONSTANT = 7.891e-11
const DEFECT_ASSISTED_RATE_CONSTANT = 6.504e-13
const TRANSITION_OMEGA = 0.9
const LOW_TEMP_THRESHOLD = 5.0
const TIME_SA_CONVERSION = 4e6 *60^2 *24 *365 /4.709197671429600e+02

# Monod scheme inhibition factors
inhibition_dO2(dO2::Float64)      = KMi_dO2   / (KMi_dO2   + dO2)
inhibition_dtNO3(dtNO3::Float64)  = KMi_dtNO3 / (KMi_dtNO3 + dtNO3)
inhibition_pMnO2(pMnO2::Float64)  = KMi_pMnO2 / (KMi_pMnO2 + pMnO2)
inhibition_pFeOH3(pFeOH3::Float64)= KMi_pFeOH3/ (KMi_pFeOH3+ pFeOH3)
inhibition_dtSO4(dtSO4::Float64)  = KMi_dtSO4 / (KMi_dtSO4 + dtSO4)

"""
Organic matter degradation pathway factors.
From the code of Couture et al. (EST 2010), following Boudreau (1996).
"""
function degradationfactors(
    dO2::Float64, dtNO3::Float64, pMnO2::Float64, pFeOH3::Float64, dtSO4::Float64,
)
    # Inhibition
    Mi_dO2   = inhibition_dO2(dO2)
    Mi_dtNO3 = inhibition_dtNO3(dtNO3)
    Mi_pMnO2 = inhibition_pMnO2(pMnO2)
    Mi_pFeOH3= inhibition_pFeOH3(pFeOH3)
    Mi_dtSO4 = inhibition_dtSO4(dtSO4)

    # Pathway Monod factors (sequential inhibition)
    fdO2   = dO2  / (KM_dO2  + dO2)
    fdtNO3 = Mi_dO2       * dtNO3 / (KM_dtNO3 + dtNO3)
    fpMnO2 = Mi_dO2*Mi_dtNO3      * pMnO2 / (KM_pMnO2 + pMnO2)
    fpFeOH3= Mi_dO2*Mi_dtNO3*Mi_pMnO2 * pFeOH3 / (KM_pFeOH3 + pFeOH3)
    fdtSO4 = Mi_dO2*Mi_dtNO3*Mi_pMnO2*Mi_pFeOH3 * dtSO4 / (KM_dtSO4 + dtSO4)
    fdCH4  = Mi_dO2*Mi_dtNO3*Mi_pMnO2*Mi_pFeOH3*Mi_dtSO4

    fox = fdO2 + fdtNO3 + fpMnO2 + fpFeOH3 + fdtSO4 + fdCH4

    # Normalize
    fdO2N    = fdO2   /= fox
    fdtNO3N  = fdtNO3 /= fox
    fpMnO2N  = fpMnO2 /= fox
    fpFeOH3N = fpFeOH3/= fox
    fdtSO4N  = fdtSO4 /= fox
    fdCH4N   = fdCH4  /= fox
    foxN     = fdO2N + fdtNO3N + fpMnO2N + fpFeOH3N + fdtSO4N + fdCH4N

    return fdO2N, fdtNO3N, fpMnO2N, fpFeOH3N, fdtSO4N, fdCH4N, foxN
end

function degrade(
    dO2::Float64, dtNO3::Float64, pMnO2::Float64, pFeOH3::Float64, dtSO4::Float64,
    pfoc_kfast::Float64, psoc_kslow::Float64,
)
    fdO2, fdtNO3, fpMnO2, fpFeOH3, fdtSO4, fdCH4, fox =
        degradationfactors(dO2, dtNO3, pMnO2, pFeOH3, dtSO4)

    Rfast_dO2    = fdO2   * pfoc_kfast
    Rslow_dO2    = fdO2   * psoc_kslow
    Rfast_dtNO3  = fdtNO3 * pfoc_kfast
    Rslow_dtNO3  = fdtNO3 * psoc_kslow
    Rfast_pMnO2  = fpMnO2 * pfoc_kfast
    Rslow_pMnO2  = fpMnO2 * psoc_kslow
    Rfast_pFeOH3 = fpFeOH3* pfoc_kfast
    Rslow_pFeOH3 = fpFeOH3* psoc_kslow
    Rfast_dtSO4  = fdtSO4 * pfoc_kfast
    Rslow_dtSO4  = fdtSO4 * psoc_kslow
    Rfast_dCH4   = fdCH4  * pfoc_kfast
    Rslow_dCH4   = fdCH4  * psoc_kslow
    Rfast_total  = fox    * pfoc_kfast
    Rslow_total  = fox    * psoc_kslow

    return (Rfast_dO2, Rslow_dO2, Rfast_dtNO3, Rslow_dtNO3,
            Rfast_pMnO2, Rslow_pMnO2, Rfast_pFeOH3, Rslow_pFeOH3,
            Rfast_dtSO4, Rslow_dtSO4, Rfast_dCH4, Rslow_dCH4,
            Rfast_total, Rslow_total)
end

# =========================
# Redox kinetics (Monod-limited)
# =========================

# Helper
@inline monod(x::Float64, K::Float64) = x / (K + max(x, 0.0))

# --- Legacy mass-action k's kept only to define Vmax defaults (units: (m^3 mol^-1 a^-1)) ---
# If you prefer, you can ignore these and set VMAX_* directly below.
const kMnII_redox    = 1e2
const kFeII_redox    = 1e2
const kNH3_redox     = 1e3
const kH2S_redox     = 3e3
const kCH4_O2redox   = 1e5
const kCH4_SO4redox  = 1e5

# --- Reference concentrations to map legacy k → Vmax (mol m^-3) ---
# Choose values typical of your oxic boundary; tweak as needed.
const CREF_O2   = 0.015   # ~15 μM
const CREF_NH4  = 0.020   # 20 μM
const CREF_H2S  = 0.020   # 20 μM
const CREF_Fe   = 0.005   # 5  μM
const CREF_Mn   = 0.005   # 5  μM
const CREF_CH4  = 0.050   # 50 μM
const CREF_SO4  = 1.000   # ~1 mol m^-3 (rarely limiting above SMTZ)

# --- Monod half-saturation constants (mol m^-3) ---
# Start near the reference values; tune to your site.
const K_O2   = CREF_O2
const K_NH4  = CREF_NH4
const K_H2S  = CREF_H2S
const K_Fe   = CREF_Fe
const K_Mn   = CREF_Mn
const K_CH4  = CREF_CH4
const K_SO4  = CREF_SO4

# --- Vmax (mol m^-3 a^-1) ---
# Set to 4 * k * Cref_sub * Cref_O2 so that at [sub]=K_sub and [O2]=K_O2
# (i.e., each Monod term = 1/2), R ≈ k * Cref_sub * Cref_O2 — matching legacy at reference conditions.
const VMAX_NH3      = 4.0 * kNH3_redox    * CREF_NH4 * CREF_O2
const VMAX_H2S      = 4.0 * kH2S_redox    * CREF_H2S * CREF_O2
const VMAX_Fe       = 4.0 * kFeII_redox   * CREF_Fe  * CREF_O2
const VMAX_Mn       = 4.0 * kMnII_redox   * CREF_Mn  * CREF_O2
const VMAX_CH4_O2   = 4.0 * kCH4_O2redox  * CREF_CH4 * CREF_O2
const VMAX_CH4_SO4  = 4.0 * kCH4_SO4redox * CREF_CH4 * CREF_SO4

"Redox reaction rates (Monod-limited). Returns the same tuple order as before."
function redox(
    dO2::Float64,
    dtNH4::Float64,
    dtH2S::Float64,
    dFeII::Float64,
    dMnII::Float64,
    dCH4::Float64,
    dtSO4::Float64,
)
    rO2 = monod(dO2, K_O2)

    R_NH3_redox   = VMAX_NH3     * monod(dtNH4, K_NH4) * rO2
    R_H2S_redox   = VMAX_H2S     * monod(dtH2S, K_H2S) * rO2
    R_FeII_redox  = VMAX_Fe      * monod(dFeII, K_Fe ) * rO2
    R_MnII_redox  = VMAX_Mn      * monod(dMnII, K_Mn ) * rO2

    R_CH4_O2redox  = VMAX_CH4_O2  * monod(dCH4, K_CH4) * rO2
    R_CH4_SO4redox = VMAX_CH4_SO4 * monod(dCH4, K_CH4) * monod(dtSO4, K_SO4)

    return R_MnII_redox, R_FeII_redox, R_NH3_redox, R_H2S_redox, R_CH4_O2redox, R_CH4_SO4redox
end

# =========================
# CaCO3 kinetics (as before; optional Ω threshold hook)
# =========================

# Set a small supersaturation threshold to delay precipitation if desired.
# Keep at 0.0 to reproduce your previous behavior.
const OMEGA_PRECIP_THR = 0.0

"CaCO3 mineral dissolution and precipitation rates."
function dissolve_precipitate_CaCO3(
    pcalcite::Float64,
    paragonite::Float64,
    dCa::Float64,
    dCO3::Float64,
    KCa::Float64,
    KAr::Float64,
    T,
    diss_scheme::Integer = 1
)
    # Saturation states
    OmegaCa = dCa * dCO3 / KCa
    OmegaAr = dCa * dCO3 / KAr

    if diss_scheme == 1
        # Base RADI kinetics (Naviaux et al. 2019 for calcite, Dong et al. 2019 for aragonite)

        # Calcite dissolution (Naviaux et al. 2019)
        Rdiss_calcite =
            (0.8275 < OmegaCa <= 1.0)  ? (pcalcite * 0.00632 * (1.0 - OmegaCa)^0.11) :
            (OmegaCa <= 0.8275)        ? (pcalcite * 20.0    * (1.0 - OmegaCa)^4.7)  :
                                        0.0

        # Aragonite dissolution (Dong et al. 2019)
        Rdiss_aragonite =
            (0.835 < OmegaAr <= 1.0) ? (paragonite * 0.0038 * (1.0 - OmegaAr)^0.13) :
            (OmegaAr <= 0.835)       ? (paragonite * 0.042  * (1.0 - OmegaAr)^1.46) :
                                    0.0

    elseif diss_scheme == 2
        # Temperature dependent calcite dissolution from fitting the data within Naviaux et al. (2019) Marine Chemistry
        if 1.0 <= OmegaCa
            Rdiss_calcite = 0.0 # Supersaturation - no dissolution

        elseif TRANSITION_OMEGA < OmegaCa < 1.0
            R1_k0 = ETCH_PIT_RATE_CONSTANT * exp(0.03298*T) * (1-TRANSITION_OMEGA)^(3.8113+0.0202*T-2)
            R2_k0 = DEFECT_ASSISTED_RATE_CONSTANT * exp(0.1055*T) * (1-TRANSITION_OMEGA)^(1.672+0.02177*T-2)

            rate1 = pcalcite * R1_k0 * (1 - OmegaCa)^2 * TIME_SA_CONVERSION
            rate2 = pcalcite * R2_k0 * (1 - OmegaCa)^2 * TIME_SA_CONVERSION

            if T <= LOW_TEMP_THRESHOLD
                Rdiss_calcite = rate1
            else
                Rdiss_calcite = max(rate1, rate2)
            end

        elseif OmegaCa <= TRANSITION_OMEGA
            etch_pit_rate = pcalcite * ETCH_PIT_RATE_CONSTANT * exp(0.03298 * T) * (1-OmegaCa)^(3.8113 + 0.0202 * T) * TIME_SA_CONVERSION
            defect_assisted_rate = pcalcite * DEFECT_ASSISTED_RATE_CONSTANT * exp(0.1055 * T) * (1-OmegaCa)^(1.672 + 0.02177 * T) * TIME_SA_CONVERSION

            if T <= LOW_TEMP_THRESHOLD
                Rdiss_calcite = etch_pit_rate
            else
                Rdiss_calcite = max(etch_pit_rate, defect_assisted_rate)
            end
        end

        # Aragonite dissolution (Dong et al. 2019)
        Rdiss_aragonite =
            (0.835 < OmegaAr <= 1.0) ? (paragonite * 0.0038 * (1.0 - OmegaAr)^0.13) :
            (OmegaAr <= 0.835)       ? (paragonite * 0.042  * (1.0 - OmegaAr)^1.46) :
                                    0.0
    end

            


    # Calcite precipitation (Zuddas & Mucci 1998), normalized to 4 m^2 g^-1
    Rprec_calcite =
        (OmegaCa > 1.0 + OMEGA_PRECIP_THR) ? (0.4075 * (OmegaCa - 1.0)^1.76) : 0.0

    # Aragonite does not precipitate in this scheme
    Rprec_aragonite = 0.0

    return Rdiss_calcite, Rdiss_aragonite, Rprec_calcite, Rprec_aragonite
end

# =========================
# Glue: totals and component rates (unchanged interfaces)
# =========================

"All reaction rates."
function getreactions(
    dO2::Float64, dtNO3::Float64, pMnO2::Float64, pFeOH3::Float64, dtSO4::Float64,
    dtNH4::Float64, dtH2S::Float64, dFeII::Float64, dMnII::Float64, dCH4::Float64,
    pfoc_kfast::Float64, psoc_kslow::Float64,
    pcalcite::Float64, paragonite::Float64, dCa::Float64, dCO3::Float64, KCa::Float64, KAr::Float64, T::Float64, diss_scheme::Integer,
)
    Rfast_dO2, Rslow_dO2, Rfast_dtNO3, Rslow_dtNO3, Rfast_pMnO2, Rslow_pMnO2,
    Rfast_pFeOH3, Rslow_pFeOH3, Rfast_dtSO4, Rslow_dtSO4, Rfast_dCH4, Rslow_dCH4,
    Rfast_total, Rslow_total =
        degrade(dO2, dtNO3, pMnO2, pFeOH3, dtSO4, pfoc_kfast, psoc_kslow)

    # Monod-limited redox
    R_dMnII, R_dFeII, R_dNH3, R_dH2S, R_CH4_O2redox, R_CH4_SO4redox =
        redox(dO2, dtNH4, dtH2S, dFeII, dMnII, dCH4, dtSO4)

    Rdiss_calcite, Rdiss_aragonite, Rprec_calcite, Rprec_aragonite =
        dissolve_precipitate_CaCO3(pcalcite, paragonite, dCa, dCO3, KCa, KAr, T, diss_scheme)

    return (Rfast_dO2, Rslow_dO2, Rfast_dtNO3, Rslow_dtNO3,
            Rfast_pMnO2, Rslow_pMnO2, Rfast_pFeOH3, Rslow_pFeOH3,
            Rfast_dtSO4, Rslow_dtSO4, Rfast_dCH4, Rslow_dCH4,
            Rfast_total, Rslow_total,
            R_dMnII, R_dFeII, R_dNH3, R_dH2S, R_CH4_O2redox, R_CH4_SO4redox,
            Rdiss_calcite, Rdiss_aragonite, Rprec_calcite, Rprec_aragonite)
end

"Convert reactions to individual component rates."
function reactions2rates(
    Rfast_dO2::Float64,   Rslow_dO2::Float64,
    Rfast_dtNO3::Float64, Rslow_dtNO3::Float64,
    Rfast_pMnO2::Float64, Rslow_pMnO2::Float64,
    Rfast_pFeOH3::Float64,Rslow_pFeOH3::Float64,
    Rfast_dtSO4::Float64, Rslow_dtSO4::Float64,
    Rfast_dCH4::Float64,  Rslow_dCH4::Float64,
    Rfast_total::Float64, Rslow_total::Float64,
    R_dMnII::Float64, R_dFeII::Float64, R_dNH3::Float64, R_dH2S::Float64,
    R_CH4_O2redox::Float64, R_CH4_SO4redox::Float64,
    Rdiss_calcite::Float64, Rdiss_aragonite::Float64,
    Rprec_calcite::Float64, Rprec_aragonite::Float64,
    phiS_phi_z::Float64, RC::Float64, RN::Float64, RP::Float64,
)
    # Convenience sums
    Rdeg_dO2    = Rfast_dO2    + Rslow_dO2
    Rdeg_dtNO3  = Rfast_dtNO3  + Rslow_dtNO3
    Rdeg_dtSO4  = Rfast_dtSO4  + Rslow_dtSO4
    Rdeg_pFeOH3 = Rfast_pFeOH3 + Rslow_pFeOH3
    Rdeg_pMnO2  = Rfast_pMnO2  + Rslow_pMnO2
    Rdeg_dCH4   = Rfast_dCH4   + Rslow_dCH4
    Rdeg_total  = Rfast_total  + Rslow_total

    # Net CaCO3 dissolution
    Rdiss_CaCO3 = Rdiss_calcite + Rdiss_aragonite - Rprec_calcite - Rprec_aragonite

    # Unit conversions
    p2d = phiS_phi_z         # particulate → dissolved
    d2p = 1.0 / phiS_phi_z   # dissolved   → particulate

    # Species tendencies
    rate_dO2   = p2d * -Rdeg_dO2 -
                 (R_dFeII/4.0 + R_dMnII/2.0 + 2.0*R_dH2S + 2.0*R_dNH3) -
                 2.0*R_CH4_O2redox
    rate_dtCO2 = p2d * (RC*(Rdeg_total - 0.5*Rdeg_dCH4) + Rdiss_CaCO3) +
                 R_CH4_O2redox + R_CH4_SO4redox
    rate_dtNO3 = p2d * RC * (-0.8*Rdeg_dtNO3) + R_dNH3
    rate_dtSO4 = p2d * RC * (-0.5*Rdeg_dtSO4) + R_dH2S - R_CH4_SO4redox
    rate_dCH4  = p2d * RC * (0.5*Rdeg_dCH4) - R_CH4_O2redox - R_CH4_SO4redox
    rate_dtPO4 = p2d * RP * Rdeg_total
    rate_dtNH4 = p2d * RN * Rdeg_total - R_dNH3
    rate_dtH2S = p2d * RC * (0.5*Rdeg_dtSO4) - R_dH2S + R_CH4_SO4redox
    rate_dFeII = p2d * RC * (4.0*Rdeg_pFeOH3) - R_dFeII
    rate_dMnII = p2d * RC * (2.0*Rdeg_pMnO2) - R_dMnII
    rate_pfoc  = -Rfast_total
    rate_psoc  = -Rslow_total
    rate_pFeOH3= RC * (-4.0*Rdeg_pFeOH3) + d2p * R_dFeII
    rate_pMnO2 = RC * (-2.0*Rdeg_pMnO2) + d2p * R_dMnII

    # TA and Ca
    rate_dCa  = p2d * Rdiss_CaCO3
    rate_dalk = p2d * ( (RN - RP)            * (Rdeg_dO2 + Rdeg_dCH4) +
                        (RN - RP + 0.8RC)    *  Rdeg_dtNO3 +
                        (RN - RP + 4.0RC)    *  Rdeg_pMnO2 +
                        (RN - RP + 8.0RC)    *  Rdeg_pFeOH3 +
                        (RN - RP + 1.0RC)    *  Rdeg_dtSO4 +
                        2.0                  *  Rdiss_CaCO3 ) +
               2.0 * R_CH4_SO4redox -
               2.0 * (R_dMnII + R_dFeII + R_dNH3 + R_dH2S)

    # CaCO3 minerals
    rate_pcalcite   = -Rdiss_calcite + Rprec_calcite
    rate_paragonite = -Rdiss_aragonite

    return (rate_dO2, rate_dtCO2, rate_dtNO3, rate_dtSO4, rate_dtPO4, rate_dtNH4,
            rate_dtH2S, rate_dFeII, rate_dMnII, rate_dCH4, rate_dalk, rate_dCa,
            rate_pfoc, rate_psoc, rate_pFeOH3, rate_pMnO2, rate_pcalcite, rate_paragonite,
            Rdeg_dO2, Rdeg_dtNO3, Rdeg_dtSO4, Rdeg_pFeOH3, Rdeg_pMnO2, Rdeg_dCH4, Rdeg_total)
end

"Rates of change of each component."
function rates(
    dO2, dtNO3, pMnO2, pFeOH3, dtSO4, dtNH4, dtH2S, dFeII, dMnII, dCH4,
    pfoc_kfast, psoc_kslow, pcalcite, paragonite, dCa, dCO3, KCa, KAr,
    phiS_phi_z, RC, RN, RP, T, diss_scheme,
)
    (Rfast_dO2, Rslow_dO2, Rfast_dtNO3, Rslow_dtNO3,
     Rfast_pMnO2, Rslow_pMnO2, Rfast_pFeOH3, Rslow_pFeOH3,
     Rfast_dtSO4, Rslow_dtSO4, Rfast_dCH4, Rslow_dCH4,
     Rfast_total, Rslow_total,
     R_dMnII, R_dFeII, R_dNH3, R_dH2S, R_dCH4_O2redox, R_dCH4_SO4redox,
     Rdiss_calcite, Rdiss_aragonite, Rprec_calcite, Rprec_aragonite) =
        getreactions(dO2, dtNO3, pMnO2, pFeOH3, dtSO4, dtNH4, dtH2S, dFeII, dMnII, dCH4,
                     pfoc_kfast, psoc_kslow, pcalcite, paragonite, dCa, dCO3, KCa, KAr, T, diss_scheme)

    return reactions2rates(Rfast_dO2, Rslow_dO2, Rfast_dtNO3, Rslow_dtNO3,
                           Rfast_pMnO2, Rslow_pMnO2, Rfast_pFeOH3, Rslow_pFeOH3,
                           Rfast_dtSO4, Rslow_dtSO4, Rfast_dCH4, Rslow_dCH4,
                           Rfast_total, Rslow_total,
                           R_dMnII, R_dFeII, R_dNH3, R_dH2S, R_dCH4_O2redox, R_dCH4_SO4redox,
                           Rdiss_calcite, Rdiss_aragonite, Rprec_calcite, Rprec_aragonite,
                           phiS_phi_z, RC, RN, RP)
end

end # module React
