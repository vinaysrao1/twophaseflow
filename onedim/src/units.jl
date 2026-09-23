# Unit conversions for volumetric rates. Internally everything is m³/s.

const BARREL = 0.158987294928  # m³ (US oil barrel)

"Convert m³/h to m³/s."
m3h(q) = q / 3600
"Convert barrels per day to m³/s."
bpd(q) = q * BARREL / 86400
"Convert m³/s to m³/h."
to_m3h(q) = q * 3600
"Convert m³/s to barrels per day."
to_bpd(q) = q * 86400 / BARREL
