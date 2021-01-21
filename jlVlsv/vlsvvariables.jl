using LaTeXStrings

# Define some units for intrinsic values
const units_predefined = Dict(
   "rhom" => "kg/m3",
   "rhoq" => "C/m3",
   "rho" => "1/m3",
   "rhobackstream" => "1/m3",
   "rhononbackstream" => "1/m3",
   "rho_v" => "1/m2s",
   "rhovbackstream" => "1/m2s",
   "rhovnonbackstream" => "1/m2s",
   "v" => "m/s",
   "vbackstream" => "m/s",
   "nNonbackstream" => "m/s",
   "b" => "T",
   "b_vol" => "T",
   "background_b" => "T",
   "perturbed_b" => "T",
   "bgb" => "T",
   "perb" => "T",
   "perb_vol" => "T",
   "e" => "V/m",
   "e_vol" => "V/m",
   "exhall_000_100" => "V/m",
   "exhall_001_101" => "V/m",
   "exhall_010_110" => "V/m",
   "exhall_011_111" => "V/m",
   "eyhall_000_010" => "V/m",
   "eyhall_001_011" => "V/m",
   "eyhall_100_110" => "V/m",
   "eyhall_101_111" => "V/m",
   "ezhall_000_001" => "V/m",
   "ezhall_010_011" => "V/m",
   "ezhall_100_101" => "V/m",
   "ezhall_110_111" => "V/m",
   "pressure" => "Pa",
   "pressure_dt2" => "Pa",
   "pressure_r" => "Pa",
   "pressure_v" => "Pa",
   "ptensordiagonal" => "Pa",
   "ptensoroffdiagonal" => "Pa",
   "ptensorbackstreamdiagonal" => "Pa",
   "ptensorbackstreamoffdiagonal" => "Pa",
   "ptensornonbackstreamdiagonal" => "Pa",
   "ptensornonbackstreamoffdiagonal" => "Pa",
   "max_v_dt" => "s",
   "max_r_dt" => "s",
   "max_fields_dt" => "s",
   "minvalue" => "s3/m6",
   "effectivesparsitythreshold" => "s3/m6",
   "rho_loss_adjust" => "1/m3",
   "energydensity" => "eV/cm3",
   "precipitationdiffflux" => "1/(cm2 sr s eV)"
)

# Define some LaTeX markup names for intrinsic values
const latex_predefined = Dict(
   "rhom" => L"$\rho_m$",
   "rhoq" => L"$\rho_q$",
   "rho" => L"$n_\mathrm{p}$",
   "rhobackstream" => L"$n_\mathrm{p,st}$",
   "rhononbackstream" => L"$n_\mathrm{p,th}$",
   "rho_v" => L"$\Gamma_\mathrm{p}$",
   "rhovbackstream" => L"$\Gamma_\mathrm{p,st}$",
   "rhovnonbackstream" => L"$\Gamma_\mathrm{p,th}$",
   "v" => L"$V$",
   "vbackstream" => L"$V_\mathrm{p,st}$",
   "vnonbackstream" => L"$V_\mathrm{p,th}$",
   "b" => L"$B$",
   "b_vol" => L"$B_\mathrm{vol}$",
   "background_b" => L"$B_\mathrm{bg}$",
   "perturbed_b" => L"$B_\mathrm{pert}$",
   "bgb" => L"$B_\mathrm{bg}$",
   "perb" => L"B_\mathrm{pert}$",
   "perb_vol" => L"B_\mathrm{vol,pert}$",
   "e" => L"$E$",
   "e_vol" => L"$E_\mathrm{vol}$",
   "exhall_000_100" => L"$E_\mathrm{Hall,000,100}$",
   "exhall_001_101" => L"$E_\mathrm{Hall,001,101}$",
   "exhall_010_110" => L"$E_\mathrm{Hall,010,110}$",
   "exhall_011_111" => L"$E_\mathrm{Hall,011,111}$",
   "eyhall_000_010" => L"$E_\mathrm{Hall,000,010}$",
   "eyhall_001_011" => L"$E_\mathrm{Hall,001,011}$",
   "eyhall_100_110" => L"$E_\mathrm{Hall,100,110}$",
   "eyhall_101_111" => L"$E_\mathrm{Hall,101,111}$",
   "ezhall_000_001" => L"$E_\mathrm{Hall,000,001}$",
   "ezhall_010_011" => L"$E_\mathrm{Hall,010,011}$",
   "ezhall_100_101" => L"$E_\mathrm{Hall,100,101}$",
   "ezhall_110_111" => L"$E_\mathrm{Hall,110,111}$",
   "pressure" => L"$P$",
   "pressure_dt2" => L"$P_{\mathrm{d}t/2}}$",
   "pressure_r" => L"$P_r$",
   "pressure_v" => L"$P_v$",
   "ptensordiagonal" => L"$\mathcal{P}_\mathrm{diag}$",
   "ptensoroffdiagonal" => L"$\mathcal{P}_\mathrm{off-diag}$",
   "ptensorbackstreamdiagonal" => L"$\mathcal{P}_\mathrm{st,diag}$",
   "ptensorbackstreamoffdiagonal" => L"$\mathcal{P}_\mathrm{st,off-diag}$",
   "ptensornonbackstreamdiagonal" => L"$\mathcal{P}_\mathrm{th,diag}$",
   "ptensornonbackstreamoffdiagonal" => L"$\mathcal{P}_\mathrm{th,off-diag}$",
   "max_v_dt" => L"$\Delta t_{\mathrm{max},v}$",
   "max_r_dt" => L"$\Delta t_{\mathrm{max},r}$",
   "max_fields_dt" => L"$\Delta t_\mathrm{max,FS}$",
   "minvalue" => L"$f_\mathrm{Min}$",
   "effectivesparsitythreshold" => L"$f_\mathrm{Min}$",
   "rho_loss_adjust" => L"$\Delta_\mathrm{loss} n_\mathrm{p}$",
)


# Define some LaTeX markup units for intrinsic values
const latexunits_predefined = Dict(
   "rhom" => L"$\mathrm{kg}\,\mathrm{m}^{-3}$",
   "rhoq" => L"$\mathrm{C}\,\mathrm{m}^{-3}$",
   "rho" => L"$\mathrm{m}^{-3}$",
   "rhobackstream" => L"$\mathrm{m}^{-3}$",
   "rhononbackstream" => L"$\mathrm{m}^{-3}$",
   "rho_v" => L"$\mathrm{m}^{-2}$s",
   "rhovbackstream" => L"$\mathrm{m}^{-2}$s",
   "rhovnonbackstream" => L"$\mathrm{m}^{-2}$s",
   "v" => L"$\mathrm{m}\,\mathrm{s}^{-1}$",
   "vbackstream" => L"$\mathrm{m}\,\mathrm{s}^{-1}$",
   "vnonbackstream" => L"$\mathrm{m}\,\mathrm{s}^{-1}$",
   "b" => L"T",
   "b_vol" => L"T",
   "background_b" => L"T",
   "perturbed_b" => L"T",
   "bgb" => L"T",
   "perb" => L"T",
   "perb_vol" => L"T",
   "e" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "e_vol" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "exhall_000_100" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "exhall_001_101" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "exhall_010_110" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "exhall_011_111" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "eyhall_000_010" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "eyhall_001_011" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "eyhall_100_110" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "eyhall_101_111" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "ezhall_000_001" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "ezhall_010_011" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "ezhall_100_101" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "ezhall_110_111" => L"$\mathrm{V}\,\mathrm{m}^{-1}$",
   "pressure" => L"Pa",
   "pressure_dt2" => L"Pa",
   "pressure_r" => L"Pa",
   "pressure_v" => L"Pa",
   "ptensordiagonal" => L"Pa",
   "ptensoroffdiagonal" => L"Pa",
   "ptensorbackstreamdiagonal" => L"Pa",
   "ptensorbackstreamoffdiagonal" => L"Pa",
   "ptensornonbackstreamdiagonal" => L"Pa",
   "ptensornonbackstreamoffdiagonal" => L"Pa",
   "max_v_dt" => L"s",
   "max_r_dt" => L"s",
   "max_fields_dt" => L"s",
   "minvalue" => L"$\mathrm{m}^{-6}\,\mathrm{s}^{3}$",
   "effectivesparsitythreshold" => L"$\mathrm{m}^{-6}\,\mathrm{s}^{3}$",
   "rho_loss_adjust" => L"$\mathrm{m}^{-3}$",
   "energydensity" => L"$\mathrm{eV}\,\mathrm{cm}^{-3}$",
   "precipitationdiffflux" => L"$\mathrm{cm}^{-2} \,\mathrm{sr}^{-1}\,\mathrm{s}^{-1}\,\mathrm{eV}^{-1}$"
)


const variables_predefined = Dict(
   "b" => meta -> sqrt.(sum(read_variable(meta, "B").^2, dims=1)),
   "e" => meta -> sqrt.(sum(read_variable(meta, "E").^2, dims=1)),
   "u" => meta -> sqrt.(sum(read_variable(meta, "vg").^2, dims=1)),
)