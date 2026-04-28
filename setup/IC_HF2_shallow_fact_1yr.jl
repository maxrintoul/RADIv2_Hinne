### Site-specific initial conditions from read_data_set_max.ipynb ###
### Site: HF2 ###

# Define model depth steps (all depths in metres)
depthSed = 0.2  # depth of the sediment
dz_top = 0.2e-3 # depth resolution at top of sediment
dz_bot = 5e-2   # depth resolution at bottom of sediment
Nz = 51         # number of depth layers

tspan = (0.0, 1.0) # in years

# +
# wave
# wave_height = 1.0 #[m]  # no corresponding value found in read_data_set_max.ipynb
# wave_period = 6.1 #[s]  # no corresponding value found in read_data_set_max.ipynb
# wavelength = (9.81*wave_period)/(2*pi) #[m]  # no corresponding value found in read_data_set_max.ipynb

depth = 25.0 #[m] average depth of enclosed bay  # no corresponding value found in read_data_set_max.ipynb
# -

permeability = 1e-13 # [m^2]  # no corresponding value found in read_data_set_max.ipynb

# Currents
U = 0.02 # m/s  # no corresponding value found in read_data_set_max.ipynb

# Define sediment porosity parameters (weighted fit with phi0 fixed to shallowest)
phiInf = 0.623115737231445
phi0 = 0.829215532415753
beta = 17.31694472932772  # converted from fitted b in cm^-1 to m^-1

# Define characteristic depths
lambda_b = 0.08  # /m  # no corresponding value found in read_data_set_max.ipynb
# lambda_f = 0.03  # /m  # no corresponding value found in read_data_set_max.ipynb
# lambda_s = 1.0   # /m  # no corresponding value found in read_data_set_max.ipynb
lambda_i = 0.05  # /m  # no corresponding value found in read_data_set_max.ipynb

# Define overlying water column properties
T_mean = 8.16577142857143 # °C
T_amp = 2.0 # °C
T_period = 1.0 # years
T_phase = 0.0 # years
T = T_mean # °C, initial temperature for chemistry calculations
S = 34.4551
P = 25.0  # dbar/depth from CTD bottom-most value
rho_sw = gsw_rho(S, T, P) # seawater density / kg/m^3

# Concentrations all in mol/kg
dO2_w = 0.000214394285714286
dtCO2_w = 0.0021287334340083655
dtNO3_w = 8.01039955872e-7
dtSO4_w = (29264.2e-6 * S / 35)  # estimated from salinity
dtNH4_w = 4.157194260107362e-6
dalk_w = 0.002325656666666667
dalk_alloch_w = copy(dalk_w)  # Tracer for allochthonous alkalinity input, e.g. from overlying water column. Assume all initial alk is allochthonous for simplicity.
dalk_aerob_w = 0.0            # Tracer for alkalinity produced from aerobic respiration
dalk_anaerob_w = 0.0          # Tracer for alkalinity produced from anaerobic respiration
dalk_carb_w = 0.0             # Tracer for alkalinity produced from carbonate dissolution

dtPO4_w = 2.2165396e-6  # no corresponding value found in read_data_set_max.ipynb
dtH2S_w = 0.0           # no corresponding value found in read_data_set_max.ipynb
dCH4_w = 0.1e-6         # no corresponding value found in read_data_set_max.ipynb
dFeII_w = 4.186537551421602e-8  # no corresponding value found in read_data_set_max.ipynb
dMnII_w = 2.1788197604572425e-7 # no corresponding value found in read_data_set_max.ipynb
dSi_w = 9.66e-6         # no corresponding value found in read_data_set_max.ipynb

# Define organic matter flux to the surface sediment
Fpom = 44.215395617325186  # no corresponding value found in read_data_set_max.ipynb
Fpom_r = 0.0               # no corresponding value found in read_data_set_max.ipynb
Fpom_s = 0.27348635644808994 # no corresponding value found in read_data_set_max.ipynb
Fpom_f = 0.7265136435519101  # no corresponding value found in read_data_set_max.ipynb
FMnO2 = 0.02 * (365.25/1000) # no corresponding value found in read_data_set_max.ipynb
FFeOH3 = 0.10 * (365.25/1000) # no corresponding value found in read_data_set_max.ipynb
Fcalcite = 0.22            # no corresponding value found in read_data_set_max.ipynb
Faragonite = 0.0           # no corresponding value found in read_data_set_max.ipynb
Fclay = (2.0 / 360.31) * 300 # no corresponding value found in read_data_set_max.ipynb
FFeS = 0.0                 # no corresponding value found in read_data_set_max.ipynb
FFeS2 = 0.0                # no corresponding value found in read_data_set_max.ipynb
FS0 = 0.0                  # no corresponding value found in read_data_set_max.ipynb
FFeOH3_PO4 = 0.5           # no corresponding value found in read_data_set_max.ipynb
rho_p = 2.65e6             # no corresponding value found in read_data_set_max.ipynb

# Define initial conditions within the sediment (scalars or arrays)
dO2_i = copy(dO2_w)
dtCO2_i = copy(dtCO2_w)
dtNO3_i = copy(dtNO3_w)
dtSO4_i = copy(dtSO4_w)
dtNH4_i = copy(dtNH4_w)
dalk_i = copy(dalk_w)

dalk_alloch_i = copy(dalk_i)    # Tracer for allochthonous alkalinity input, e.g. from overlying water column. Assume all initial alk is allochthonous for simplicity.
dalk_aerob_i = 0.0              # Tracer for alkalinity produced from aerobic respiration
dalk_anaerob_i = 0.0            # Tracer for alkalinity produced from anaerobic respiration
dalk_carb_i = 0.0               # Tracer for alkalinity produced from carbonate dissolution

dtPO4_i = copy(dtPO4_w)   # no corresponding value found in read_data_set_max.ipynb
dtH2S_i = copy(dtH2S_w)   # no corresponding value found in read_data_set_max.ipynb
dFeII_i = copy(dFeII_w)   # no corresponding value found in read_data_set_max.ipynb
dMnII_i = copy(dMnII_w)   # no corresponding value found in read_data_set_max.ipynb
dCH4_i = copy(dCH4_w)     # no corresponding value found in read_data_set_max.ipynb

pfoc_i = 0.0              # no corresponding value found in read_data_set_max.ipynb
psoc_i = 3e2              # no corresponding value found in read_data_set_max.ipynb
proc_i = 6e2              # no corresponding value found in read_data_set_max.ipynb
pFeOH3_i = 0.0            # no corresponding value found in read_data_set_max.ipynb
pMnO2_i = 0.0             # no corresponding value found in read_data_set_max.ipynb
pcalcite_i = 4.5e3        # no corresponding value found in read_data_set_max.ipynb
paragonite_i = 0.0        # no corresponding value found in read_data_set_max.ipynb
pclay_i = 1.0e3           # no corresponding value found in read_data_set_max.ipynb
pFeS_i = 0.0              # no corresponding value found in read_data_set_max.ipynb
pFeS2_i = 0.0             # no corresponding value found in read_data_set_max.ipynb
pS0_i = 0.0               # no corresponding value found in read_data_set_max.ipynb
pFeOH3_PO4_i = 0.0        # no corresponding value found in read_data_set_max.ipynb

# Define dissolution scheme for CaCO3 (1 = base RADI, 2 = temp dependent dissolution)
calcite_diss_scheme = 2   # no corresponding value found in read_data_set_max.ipynb
aragonite_diss_scheme = 2 # no corresponding value found in read_data_set_max.ipynb

# Define precipitation scheme for calcite (1 = base RADI, 2 = temp dependent precipitation)
calcite_precip_scheme = 2 # no corresponding value found in read_data_set_max.ip

# Define Q10 for secondary reactions
Q10_primary = 2.0   # no corresponding value found in read_data_set_max.ipynb
Q10_secondary = 3.8 # no corresponding value found in read_data_set_max.ipynb
Tref = 9.0          # no corresponding value found in read_data_set_max.ipynb

# Optional factorial controls
factorial_T_levels        = [10.0]
factorial_Fpom_levels     = [33.52618867924528*0.1, 33.52618867924528, 33.52618867924528*3]
factorial_U_levels        = [U]
factorial_P_levels        = [P]
factorial_Fcalcite_levels = [0.01, 0.1, 1.0]
factorial_O_levels        = [20, 100, 200, 400] * 1e-6
