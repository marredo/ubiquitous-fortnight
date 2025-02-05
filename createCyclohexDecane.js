/***************************************************************
 * createBenzToluene.js
 *
 * Adiabatic flash for benzene (component 1) + toluene (component 2)
 * using a 10-mole feed basis, stepping T by 1°C, and producing
 * a multi-line JSON in the "old solutions" style.
 *
 * WARNING: This can still be quite large, so ensure you
 * have enough memory/disk if you keep small steps.
 *
 ***************************************************************/

const fs = require('fs');

//-------------- 1) Thermo Data for Benzene & Toluene ---------------
const benzene = {
  name: 'Benzene',
  // Antoine constants (°C, mmHg)
  A: 6.880,
  B: 1197,
  C: 219.2,

  // Heat of vaporization (kJ/mol)
  Hvap_kJmol: 30.7,

  // Liquid Cp (J/mol·K), approx constant
  CpL_JmolK: 136,

  // Vapor Cp polynomial: A + B*T + C*T^2 + D*T^3
  // from your data: A=-3.39E+01, B=4.74E-01, C=-3.02E-04, D=7.13E-08
  cpV_JmolK(T) {
    const A=-3.39e+01, B=4.74e-01, C=-3.02e-04, D=7.13e-08;
    return A + B*T + C*T**2 + D*T**3;  // J/mol·K
  }
};

const toluene = {
  name: 'Toluene',
  A: 6.951,
  B: 1342,
  C: 219.2,

  Hvap_kJmol: 51.0,
  CpL_JmolK: 157,

  // Vapor polynomial from your snippet:
  // A=-2.44E+01, B=5.13E-01, C=-2.77E-04, D=4.91E-08
  cpV_JmolK(T) {
    const A=-2.44e+01, B=5.13e-01, C=-2.77e-04, D=4.91e-08;
    return A + B*T + C*T**2 + D*T**3;
  }
};

//------------- 2) Convert T(°C) => Psat(bar) via Antoine -------------
function antoineBar(T_C, comp) {
  // log10(P_mmHg) = A - B/(T_C + C)
  // 1 mmHg ~ 0.00133322 bar
  const log10P = comp.A - comp.B/(T_C + comp.C);
  const p_mmHg = Math.pow(10, log10P);
  return p_mmHg * 0.00133322;
}

//------------- 3) Solve for Vapor Fraction (Beta) -------------
function solveVaporFraction(T_C, Pbar, z1) {
  const K1 = antoineBar(T_C, benzene)/Pbar;
  const K2 = antoineBar(T_C, toluene)/Pbar;
  let betaLo=0, betaHi=1, beta=0.5;

  for(let i=0; i<30; i++){
    const fLo= rrFunc(z1, 1-z1, K1, K2, betaLo);
    const fHi= rrFunc(z1, 1-z1, K1, K2, betaHi);
    beta= 0.5*(betaLo + betaHi);
    const fMid= rrFunc(z1, 1-z1, K1, K2, beta);

    if(fLo*fMid<0) {
      betaHi= beta;
    } else {
      betaLo= beta;
    }
    if(Math.abs(fMid)<1e-12) break;
  }
  return { beta, K1, K2 };
}

// Rachford–Rice function
function rrFunc(z1, z2, K1, K2, beta){
  const t1= z1*(K1-1)/(1+ beta*(K1-1));
  const t2= z2*(K2-1)/(1+ beta*(K2-1));
  return t1 + t2;
}

//------------- 4) Integrate Vapor Cp Polynomials -------------
function getCoeffs(comp){
  // We'll replicate the polynomial's A,B,C,D for indefinite integration
  if(comp === benzene) {
    return { A:-3.39e+01, B:4.74e-01, C:-3.02e-04, D:7.13e-08 };
  } else {
    return { A:-2.44e+01, B:5.13e-01, C:-2.77e-04, D:4.91e-08 };
  }
}

// Indef integral => A*T + B/2*T^2 + C/3*T^3 + D/4*T^4
function integrateCpV(comp, Tlo, Thi){
  const c = getCoeffs(comp);
  function iFunc(T) {
    return (
      c.A*T
      + (c.B/2)*T*T
      + (c.C/3)*T**3
      + (c.D/4)*T**4
    );
  }
  return iFunc(Thi) - iFunc(Tlo);
}

//------------- 5) Enthalpy Calculations -------------
function enthalpyLiquid(TK, x1){
  // Weighted by x1, x2 => 1-x1
  const x2 = 1 - x1;
  const Cp_mix = x1*benzene.CpL_JmolK + x2*toluene.CpL_JmolK;
  const Tref= 298;  // K
  return Cp_mix*(TK - Tref);  // J/mol
}

function enthalpyVapor(TK, y1){
  const y2= 1- y1;
  const Tref= 298;

  // benzene
  const Hv_benz_J = benzene.Hvap_kJmol*1000; 
  const iBenz = integrateCpV(benzene, Tref, TK);
  const H_benz = Hv_benz_J + iBenz;

  // toluene
  const Hv_tol_J = toluene.Hvap_kJmol*1000;
  const iTol = integrateCpV(toluene, Tref, TK);
  const H_tol = Hv_tol_J + iTol;

  return y1*H_benz + y2*H_tol;
}

// feed enthalpy for 10 moles
function feedEnthalpy(T_C, z1, Pbar){
  // We'll do an isothermal flash at (T_C,Pbar) to see if feed is partial, etc.
  const eq= solveVaporFraction(T_C, Pbar, z1);
  if(eq.beta<0){
    // all liquid
    return 10* enthalpyLiquid(T_C+273.15, z1);
  } else if(eq.beta>1){
    // all vapor
    return 10* enthalpyVapor(T_C+273.15, z1);
  } else {
    // partial
    const Lfrac= 1- eq.beta;
    const denom= 1 + eq.beta*(eq.K1-1);
    const x1= z1/ denom;
    const y1= eq.K1*x1;
    const hL= enthalpyLiquid(T_C+273.15, x1);
    const hV= enthalpyVapor(T_C+273.15, y1);
    return 10*(Lfrac*hL + eq.beta*hV);
  }
}

//------------- 6) Adiabatic Flash Solver -------------
function adiabaticFlash(Tfeed_C, z1, Pbar){
  const Hfeed= feedEnthalpy(Tfeed_C, z1, Pbar);

  // bracket Tflash
  let Tlow= Math.min(Tfeed_C, 0), Thigh= Math.max(Tfeed_C, 500);
  let TflashC= 0.5*(Tlow + Thigh);

  let L=0, V=0, x1=z1, y1=z1;

  for(let i=0; i<30; i++){
    const eq= solveVaporFraction(TflashC, Pbar, z1);
    const Lfr= 1 - eq.beta;
    const denom= 1 + eq.beta*(eq.K1 - 1);
    const xf= z1/ denom;
    const yf= eq.K1* xf;

    const hL= enthalpyLiquid(TflashC+273.15, xf);
    const hV= enthalpyVapor(TflashC+273.15, yf);
    const Hout= 10*Lfr*hL + 10*eq.beta*hV;

    const diff= Hout - Hfeed;
    if(Math.abs(diff)<1e-2){
      L= 10*Lfr; V= 10*eq.beta; x1= xf; y1= yf;
      break;
    }
    if(diff>0){
      Thigh= TflashC;
    } else {
      Tlow= TflashC;
    }
    TflashC= 0.5*(Tlow + Thigh);

    if(i===29){
      L= 10*Lfr; V= 10*eq.beta; x1= xf; y1= yf;
    }
  }

  // clamp single-phase
  if(V<0){
    V=0; L=10; TflashC= Tfeed_C; x1=z1; y1=z1;
  } else if(V>10){
    V=10; L=0; TflashC= Tfeed_C; x1=z1; y1=z1;
  }
  // Return [Tfeed, z, P, L, V, Tflash, x1, y1]
  return [ Tfeed_C, z1, Pbar, L, V, TflashC, x1, y1 ];
}

//------------- 7) Build Large Array ( T in 1°C steps ) -------------
const T_MIN=120, T_MAX=300, T_STEP=1;        // 1°C increments
const Z_MIN=0.0, Z_MAX=1.0,  Z_STEP=0.01;    // => 101 points
const P_MIN=0.25, P_MAX=4.0, P_STEP=0.25;    // => 16 points

const results=[];

for(let T=T_MIN; T<=T_MAX+1e-9; T+=T_STEP){
  const Tf= +T.toFixed(6);  // ensure consistent float
  for(let z=Z_MIN; z<=Z_MAX+1e-9; z+=Z_STEP){
    const zf= +z.toFixed(6);
    for(let P=P_MIN; P<=P_MAX+1e-9; P+=P_STEP){
      const Pf= +P.toFixed(6);

      results.push( adiabaticFlash(Tf, zf, Pf) );
    }
  }
}

//------------- 8) Old Solutions Multi-Line Format -------------
function oldStyleFormat(rows){
  let s= "[\n";
  for(let i=0; i<rows.length; i++){
    s += "\t[\n";
    const row= rows[i];
    for(let j=0; j<row.length; j++){
      s += "\t\t" + row[j];
      if(j<row.length-1) s+=",";
      s += "\n";
    }
    s+=(i<rows.length-1) ? "\t],\n" : "\t]\n";
  }
  s+="]";
  return s;
}

// Final
console.log("Constructed array with", results.length, "rows at 1°C steps for Benzene/Toluene!");
const outStr= oldStyleFormat(results);
fs.writeFileSync("benztol_solutions.json", outStr, "utf8");
console.log("Done writing benztol_solutions.json!");
