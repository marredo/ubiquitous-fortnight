const fs = require('fs');

//--------- 1) Thermo Data for Cyclohexane & n-Decane -----------
const cyclohexane = {
  name: 'Cyclohexane',
  // Antoine constants (°C, mmHg)
  A: 6.851,
  B: 1206,
  C: 223.1,
  // Heat of vaporization (kJ/mol)
  Hvap_kJmol: 30.8,
  // Liquid Cp (J/mol·K)
  CpL_JmolK: 156,
  // Vapor Cp polynomial (T in K): A + B*T + C*T^2 + D*T^3
  cpV_JmolK(T) {
    const A = -5.45e+01,
          B = 6.11e-01,
          C = -2.52e-04,
          D = 1.32e-08; // fixed
    return A + B*T + C*T**2 + D*T**3;
  }
};

const decane = {
  name: 'n-Decane',
  A: 6.944,
  B: 1495,
  C: 193.9,
  // Updated heat of vaporization
  Hvap_kJmol: 51.3,  
  // Liquid Cp (J/mol·K)
  CpL_JmolK: 313,
  cpV_JmolK(T) {
    const A = -7.913,
          B = 9.61e-01,
          C = -5.29e-04,
          D = 1.13e-07; // fixed
    return A + B*T + C*T**2 + D*T**3;
  }
};

//--------- 2) Antoine: T(°C) => Psat (bar) -----------
function antoineBar(T_C, comp) {
  // log10(P_mmHg) = A - B/(T_C + C)
  const log10P = comp.A - comp.B / (T_C + comp.C);
  const p_mmHg = Math.pow(10, log10P);
  return p_mmHg * 0.00133322;  // convert mmHg -> bar
}

//--------- 3) Bubble & Dew Point Calculations -----------
function bubblePoint(z1, Pbar) {
  let T_low = 50, T_high = 300, T_mid;
  for (let i = 0; i < 50; i++) {
    T_mid = 0.5 * (T_low + T_high);
    const f = z1 * antoineBar(T_mid, cyclohexane)
            + (1 - z1) * antoineBar(T_mid, decane)
            - Pbar;
    if (Math.abs(f) < 1e-6) break;
    if (f > 0) T_high = T_mid;
    else T_low = T_mid;
  }
  return T_mid;
}

function dewPoint(z1, Pbar) {
  let T_low = 50, T_high = 300, T_mid;
  for (let i = 0; i < 50; i++) {
    T_mid = 0.5 * (T_low + T_high);
    const f = z1 * (Pbar / antoineBar(T_mid, cyclohexane))
            + (1 - z1) * (Pbar / antoineBar(T_mid, decane))
            - 1;
    if (Math.abs(f) < 1e-6) break;
    if (f > 0) T_low = T_mid;
    else T_high = T_mid;
  }
  return T_mid;
}

//--------- 4) Cp Integration for Vapor -----------
function getCoeffs(comp) {
  if (comp === cyclohexane) {
    return { A: -5.45e+01, B: 6.11e-01, C: -2.52e-04, D: 1.32e-08 };
  } else {
    return { A: -7.913, B: 9.61e-01, C: -5.29e-04, D: 1.13e-07 };
  }
}

function integrateCpV(comp, Tlo, Thi) {
  const c = getCoeffs(comp);
  const iFunc = T => (
    c.A * T +
    (c.B / 2) * T**2 +
    (c.C / 3) * T**3 +
    (c.D / 4) * T**4
  );
  return iFunc(Thi) - iFunc(Tlo);
}

//--------- 5) Enthalpy Calculations -----------
function enthalpyLiquid(TK, x1) {
  const x2 = 1 - x1;
  const Cp_mix = x1 * cyclohexane.CpL_JmolK + x2 * decane.CpL_JmolK;
  const Tref = 298;
  return Cp_mix * (TK - Tref);
}

function enthalpyVapor(TK, y1) {
  const y2 = 1 - y1;
  const Tref = 298;
  
  const Hv_cy = cyclohexane.Hvap_kJmol * 1000;
  const iCy = integrateCpV(cyclohexane, Tref, TK);
  const H_cy = Hv_cy + iCy;
  
  const Hv_dec = decane.Hvap_kJmol * 1000;
  const iDec = integrateCpV(decane, Tref, TK);
  const H_dec = Hv_dec + iDec;
  
  return y1 * H_cy + y2 * H_dec;
}

//--------- 6) Feed Enthalpy via Isobaric Path -----------
function feedEnthalpyPath(T_C, z1, Pbar) {
  const n = 10;
  const T_feed_K = T_C + 273.15;
  const T_ref = 298;
  const T_b_C = bubblePoint(z1, Pbar);
  const T_b_K = T_b_C + 273.15;
  const T_d_C = dewPoint(z1, Pbar);
  const T_d_K = T_d_C + 273.15;
  
  let H;
  if (T_feed_K <= T_b_K) {
    const Cp_mix = z1 * cyclohexane.CpL_JmolK + (1 - z1) * decane.CpL_JmolK;
    H = n * Cp_mix * (T_feed_K - T_ref);
  } else if (T_feed_K >= T_d_K) {
    const Cp_liq = z1 * cyclohexane.CpL_JmolK + (1 - z1) * decane.CpL_JmolK;
    const h_liq = Cp_liq * (T_b_K - T_ref);
    
    const latent_cy = cyclohexane.Hvap_kJmol * 1000
                    + integrateCpV(cyclohexane, T_b_K, T_d_K)
                    - cyclohexane.CpL_JmolK * (T_d_K - T_b_K);
    const latent_dec = decane.Hvap_kJmol * 1000
                     + integrateCpV(decane, T_b_K, T_d_K)
                     - decane.CpL_JmolK * (T_d_K - T_b_K);
    const latent_mix = z1 * latent_cy + (1 - z1) * latent_dec;
    
    const h_vap_heat = enthalpyVapor(T_feed_K, z1) - enthalpyVapor(T_d_K, z1);
    H = n * (h_liq + latent_mix + h_vap_heat);
  } else {
    const H_liq = enthalpyLiquid(T_b_K, z1);
    const H_vap = enthalpyVapor(T_b_K, z1);
    const frac = (T_feed_K - T_b_K) / (T_d_K - T_b_K);
    const H_mix = H_liq + frac * (H_vap - H_liq);
    H = n * H_mix;
  }
  return H;
}

//--------- 7) Adiabatic Flash Solver (No Rachford–Rice) -----------
function adiabaticFlash(Tfeed_C, z1, Pbar) {
  const n = 10;
  const Tfeed_K = Tfeed_C + 273.15;
  const T_b = bubblePoint(z1, Pbar);
  const T_d = dewPoint(z1, Pbar);
  const T_b_K = T_b + 273.15;
  const T_d_K = T_d + 273.15;
  const Hfeed = feedEnthalpyPath(Tfeed_C, z1, Pbar);

  if (Tfeed_K <= T_b_K) {
    return [Tfeed_C, z1, Pbar, n, 0, Tfeed_C, z1, z1];
  } else if (Tfeed_K >= T_d_K) {
    return [Tfeed_C, z1, Pbar, 0, n, Tfeed_C, z1, z1];
  } else {
    const T_flash_C = T_b;
    const T_flash_K = T_flash_C + 273.15;
    const H_liq_sat = enthalpyLiquid(T_flash_K, z1);
    const H_vap_sat = enthalpyVapor(T_flash_K, z1);
    
    const beta = (Hfeed / n - H_liq_sat) / (H_vap_sat - H_liq_sat);
    const L = n * (1 - beta);
    const V = n * beta;
    return [Tfeed_C, z1, Pbar, L, V, T_flash_C, z1, z1];
  }
}

//--------- 8) Build Results Array -----------
const T_MIN = 120, T_MAX = 300, T_STEP = 1;
const Z_MIN = 0.0, Z_MAX = 1.0, Z_STEP = 0.01;
const P_MIN = 0.25, P_MAX = 4.0, P_STEP = 0.25;

const results = [];
for (let T = T_MIN; T <= T_MAX + 1e-9; T += T_STEP) {
  const Tf = +T.toFixed(6);
  for (let z = Z_MIN; z <= Z_MAX + 1e-9; z += Z_STEP) {
    const zf = +z.toFixed(6);
    for (let P = P_MIN; P <= P_MAX + 1e-9; P += P_STEP) {
      const Pf = +P.toFixed(6);
      results.push(adiabaticFlash(Tf, zf, Pf));
    }
  }
}

function oldStyleFormat(rows) {
  let s = "[\n";
  for (let i = 0; i < rows.length; i++) {
    s += "\t[\n";
    const row = rows[i];
    for (let j = 0; j < row.length; j++) {
      s += "\t\t" + row[j];
      if (j < row.length - 1) s += ",";
      s += "\n";
    }
    s += (i < rows.length - 1) ? "\t],\n" : "\t]\n";
  }
  s += "]";
  return s;
}

console.log("Constructed array with", results.length, "rows for Cyclohexane/n-Decane (adiabatic, no R-R)!");
const outStr = oldStyleFormat(results);
fs.writeFileSync("cyclohex_decane_solutions.json", outStr, "utf8");
console.log("Done writing cyclohex_decane_solutions.json!");
