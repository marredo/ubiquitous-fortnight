//°C
window.g = {
    cnv: undefined,
    flashDrumColor: [0, 200, 200],  // Default cyan color

    // Slider inputs
    P: .25,
    T: 271,
    moleFrac: .5,

    // Colors to be used repeatedly
    blue: [0,100,200],
    green: [0,200,0],

    // Constants
    Tref: 25,
    cpL1: .11,
    cpL2: .076,
    cpV1: .052,
    cpV2: .034,
    dHvap1: 35.3,
    dHvap2: 40.656,
    nf: 10,

    mixtures: {
        A: {
            comp1: { 
                name: "Methanol",
                mw: 32.04, 
                antoine: { A: 8.081, B: 1582.271, C: 239.726 } 
            },
            comp2: { 
                name: "Water",
                mw: 18.02, 
                antoine: { A: 8.071, B: 1730.630, C: 233.426 } 
            }
        },
        B: {
            comp1: {
                name: "n-Hexane",
                mw: 86.18,
                antoine: { A: 6.876, B: 1171, C: 224.4 }
            },
            comp2: {
                name: "n-Octane",
                mw: 114.23,
                antoine: { A: 6.919, B: 1352, C: 209.2 }
            }
        },
        C: {
            comp1: {
                name: "Benzene",
                mw: 78.11,
                antoine: { A: 6.880, B: 1197, C: 219.2 }
            },
            comp2: {
                name: "Toluene", 
                mw: 92.14,
                antoine: { A: 6.951, B: 1342, C: 219.2 }
            }
        },
        D: {
            comp1: {
                name: "Cyclohexane",
                mw: 84.16,
                antoine: { A: 6.851, B: 1206, C: 223.1 }
            },
            comp2: {
                name: "n-Decane",
                mw: 142.29,
                antoine: { A: 6.944, B: 1495, C: 193.9 }
            }
        }
    },
    showCharts: false,
}

let solutionsA;
let solutionsB;
let solutionsC;
let solutionsD;
function preload(){
    solutionsA = loadJSON("solutions.json");
    solutionsB = loadJSON("hexdec_solutions_oldstyle.json");
    solutionsC = loadJSON("benztol_solutions_big.json");
    solutionsD = loadJSON("cyclohex_decane_solutions_big.json");
}
// Moles in feed
mf = {
    lx: 50,
    rx: 200,
    by: 310,
    ty: 190,
}

// Moles in vapor
mv = {
    lx: 500,
    rx: 650,
    by: 170,
    ty: 50,
}

// Moles in liquid
ml = {
    lx: 500,
    rx: 650,
    by: 450,
    ty: 330,    
}

function setup(){
    g.cnv = createCanvas(700,500);
    g.cnv.parent("graphics-wrapper");
    document.getElementsByTagName("main")[0].remove();
}

function draw(){
    background(250);
    frame();
    mathAndDisplay();
}

// Event listeners and such
const pressure = document.getElementById("drum-pressure");
const pressureLabel = document.getElementById("drum-pressure-value");
const temperature = document.getElementById("feed-temp");
const temperatureLabel = document.getElementById("feed-temp-value");
const moleFrac = document.getElementById("methanol-feed");
const moleFracLabel = document.getElementById("methanol-feed-value");

// Add event listener for GC button
const gcButton = document.getElementById("gc-button");
gcButton.addEventListener("click", function() {
    g.showCharts = true;
    redraw();
});

// Modify existing slider event listeners
temperature.addEventListener("input", function(){
    const temp = Number(temperature.value);
    g.T = temp;
    g.showCharts = false;  // Hide charts on slider change
    temperatureLabel.innerHTML = `${temp}`;
    setTimeout(redraw(),2000);
});

moleFrac.addEventListener("input", function(){
    const temp = Number(moleFrac.value);
    g.moleFrac = temp;
    g.showCharts = false;  // Hide charts on slider change
    moleFracLabel.innerHTML = `${temp}`;
    setTimeout(redraw(),2000);
});

pressure.addEventListener("input", function(){
    const temp = Number(pressure.value);
    g.P = temp;
    g.showCharts = false;  // Hide charts on slider change
    pressureLabel.innerHTML = `${temp}`;
    setTimeout(redraw(),3000);
});

const mixtureSelect = document.getElementById("mixture-select");
const inputs = document.querySelectorAll(".mixture-data input");

mixtureSelect.addEventListener("change", function() {
    const mixture = g.mixtures[this.value];
    if (!mixture) {
        inputs.forEach(input => {
            input.value = "";
            input.disabled = true;
        });
        // Default colors
        g.green = [0, 200, 0];
        g.blue = [0, 100, 200];
        return;
    }

    // Set colors based on mixture
    switch(this.value) {
        case 'B':
            g.green = [0, 200, 200];     // cyan
            g.blue = [200, 0, 200];      // magenta
            g.flashDrumColor = [234, 211, 237];  // updated color
            break;
        case 'C':
            g.green = [255, 220, 0];     // yellow
            g.blue = [255, 130, 0];      // orange
            g.flashDrumColor = [255, 232, 181];  // updated color
            break;
        case 'D':
            g.green = [255, 60, 150];    // pink
            g.blue = [150, 0, 200];      // violet
            g.flashDrumColor = [248, 209, 241];  // updated color
            break;
        default:
            g.green = [0, 200, 0];       // default green
            g.blue = [0, 100, 200];      // default blue
            g.flashDrumColor = [0, 200, 200];  // default cyan
    }
    
    // Update values
    inputs[0].value = mixture.comp1.mw;
    inputs[1].value = mixture.comp1.antoine.A;
    inputs[2].value = mixture.comp1.antoine.B;
    inputs[3].value = mixture.comp1.antoine.C;

    inputs[4].value = mixture.comp2.mw;
    inputs[5].value = mixture.comp2.antoine.A;
    inputs[6].value = mixture.comp2.antoine.B;
    inputs[7].value = mixture.comp2.antoine.C;
});

// Trigger change event for mixture select to initialize with Mixture A
mixtureSelect.value = 'A';
mixtureSelect.dispatchEvent(new Event('change'));

function mathSolve(){
    const mixture = g.mixtures[mixtureSelect.value];
    if (!mixture) return;
    
    // Check if mixture has all required properties
    const requiredProps = ['cpL', 'cpV', 'dHvap'];
    const missingProps = {
        comp1: requiredProps.filter(prop => !mixture.comp1[prop]),
        comp2: requiredProps.filter(prop => !mixture.comp2[prop])
    };
    
    if (missingProps.comp1.length > 0 || missingProps.comp2.length > 0) {
        return ["NUMBERS", "NEEDED", "---", "---", "---"];
    }
    
    // If all properties exist, proceed with normal calculation
    let ind = Math.round((g.T-120)*1616 + 100*g.moleFrac*16 + (g.P/.25 - 1));
    let t = data[ind];
    return([t[3],t[4],t[5],t[6],t[7]]);
}