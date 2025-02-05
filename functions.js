function frame(){
    push(); fill(250); strokeWeight(1.3);
    rect(50,height/2-60,150,120); // Moles in feed
    rect(500,50,150,120); // Moles in vapor
    rect(500,330,150,120); // Moles in liquid
    pop();

    push(); // Flash drum
    beginShape(); strokeWeight(2); 
    fill(g.flashDrumColor[0], g.flashDrumColor[1], g.flashDrumColor[2], 100);
    for(let i = 0; i < flashDrum.length; i++){
        vertex(265+.9*flashDrum[i][0],120+.8*flashDrum[i][1]);
    }
    endShape();
    pop();

    push();
    noStroke(); textSize(20);
    text('moles in feed',mf.lx+12,mf.ty-8); // Box labels
    text('moles in vapor',mv.lx+11,mv.ty-8);
    text('moles in liquid',ml.lx+10,ml.ty-8);
    text('flash drum',width/2-43,height/2+6);
    text(g.P.toFixed(2)+' bar',width/2-33,height/2+30);
    text('liquid feed: '+g.T+'°C',mf.lx,mf.by+50);
    textSize(16);

    // Get the current mixture
    const mixture = g.mixtures[mixtureSelect.value];
    const comp1Name = mixture ? mixture.comp1.name : 'methanol';
    const comp2Name = mixture ? mixture.comp2.name : 'water';

    // Update component labels
    if (mixtureSelect.value === 'D') {
        // Save current text size
        let currentSize = textSize();
        textSize(14); // Smaller text for Cyclohexane
        text(comp1Name, mf.lx+5, mf.by+18);
        text(comp1Name, mv.lx+5, mv.by+18);
        text(comp1Name, ml.lx+5, ml.by+18);
        
        textSize(16); // Reset text size for n-Decane
        text(comp2Name, mf.rx-56, mf.by+18);
        text(comp2Name, mv.rx-56, mv.by+18);
        text(comp2Name, ml.rx-56, ml.by+18);
        
        // Restore original text size
        textSize(currentSize);
    } else {
        text(comp1Name, mf.lx+5, mf.by+18);
        text(comp2Name, mf.rx-56, mf.by+18);
        text(comp1Name, mv.lx+5, mv.by+18);
        text(comp2Name, mv.rx-56, mv.by+18);
        text(comp1Name, ml.lx+5, ml.by+18);
        text(comp2Name, ml.rx-56, ml.by+18);
    }
    pop();

    // Lines and arrows
    push(); strokeWeight(2);
    line(mf.rx+10,height/2,mf.rx+80,height/2); // MF to flashdrum
    arrow([mf.rx+10,height/2],[mf.rx+85,height/2],0,17,5);
    line(353,134,353,110); // Top of flashdrum to MV
    line(353,110,mv.lx-10,110);
    arrow([353,110],[mv.lx-5,110],0,17,5);
    line(353,363,353,390); // Bottom of flashdrum to ML
    line(353,390,ml.lx-10,390);
    arrow([353,390],[ml.lx-5,390],0,17,5);
    pop();
}

function mathAndDisplay(){
  // Moles in feed dipslay
  molesInFeed();
  function molesInFeed(){
      push();
      strokeWeight(.5); fill(g.green);
      rect(mf.lx+10,mf.by,60,map(g.moleFrac,0,1,0,-95));
      fill(g.blue);
      rect(mf.lx+80,mf.by,60,map(1-g.moleFrac,0,1,0,-95));
      pop();
      push();
      noStroke(); textSize(18);
      if(g.moleFrac == 1){
          text((10*g.moleFrac).toFixed(1),mf.lx+21,mf.by+map(g.moleFrac,0,1,0,-95)-5);
          text((10*(1-g.moleFrac)).toFixed(1),mf.lx+96,mf.by+map(1-g.moleFrac,0,1,0,-95)-5);
      } else if(g.moleFrac == 0){
          text((10*g.moleFrac).toFixed(1),mf.lx+26,mf.by+map(g.moleFrac,0,1,0,-95)-5);
          text((10*(1-g.moleFrac)).toFixed(1),mf.lx+91,mf.by+map(1-g.moleFrac,0,1,0,-95)-5);
      } else {
          text((10*g.moleFrac).toFixed(1),mf.lx+26,mf.by+map(g.moleFrac,0,1,0,-95)-5);
          text((10*(1-g.moleFrac)).toFixed(1),mf.lx+96,mf.by+map(1-g.moleFrac,0,1,0,-95)-5);
      }
      pop();
  }
  
  let answers = mathSolve(); // L, V, t, x, y
  if(answers[3] < 0){
    answers[3] = 0;
  } 
  if (answers[4] < 0){
    answers[4] = 0;
  }
  if(answers[0] > 10){
    answers[0] = 10;
  } else if (answers[0] < 0){
    answers[0] = 0;
  }
  if(answers[1] > 10){
    answers[1] = 10;
  } else if (answers[1] < 0){
    answers[1] = 0;
  }
  answerDisplay();
  function answerDisplay(){
    push();
    noStroke(); textSize(20);
    text(answers[2].toFixed(0)+'°C',width/2-20,height/2-15);
    
    // Only show these if showCharts is true
    if (g.showCharts) {
        text('vapor = '+answers[1].toFixed(1)+' mol',348,60);
        text('liquid = '+answers[0].toFixed(1)+' mol',348,ml.ty+85);
        text(' = '+answers[4].toFixed(2),399,90);
        text(' = '+answers[3].toFixed(2),395,ml.ty+115);

        textStyle(ITALIC);
        text('y',378,90);
        text('x',375,ml.ty+115);

        textStyle(NORMAL); textSize(16);
        text('m',388,95);
        text('m',385,ml.ty+120);
        
        // Draw the vapor flow rate text
        fill(0, 100, 200);  
        textStyle(BOLD);    
        text('Outlet Vapor Molar Flow Rate = ' + (answers[1]/10 * 100).toFixed(1) + ' mol/h', mv.lx-45, mv.by+35);
        text('Outlet Liquid Molar Flow Rate = ' + (answers[0]/10 * 100).toFixed(1) + ' mol/h', ml.lx-45, ml.by+35);
    }
    pop();

    // Only draw the charts if showCharts is true
    if (g.showCharts) {
        push();
        strokeWeight(.5); fill(g.green);
        rect(mv.lx+10,mv.by,60,map(answers[1]*answers[4],0,10,0,-95));
        rect(ml.lx+10,ml.by,60,map(answers[0]*answers[3],0,10,0,-95));
        fill(g.blue);
        rect(mv.lx+80,mv.by,60,map(answers[1]*(1-answers[4]),0,10,0,-95));
        rect(ml.lx+80,ml.by,60,map(answers[0]*(1-answers[3]),0,10,0,-95));
        pop();

        push();
        noStroke(); textSize(18);
        let t1 = answers[1]*(1-answers[4]);
        if(t1 < 0) t1 = 0;
        let t2 = answers[0]*(1-answers[3]);
        if(t2 < 0) t2 = 0;
        text((answers[1]*answers[4]).toFixed(1),mv.lx+26,mv.by+map(answers[1]*answers[4],0,10,0,-95)-5);
        text(t1.toFixed(1),mv.lx+96,mv.by+map(answers[1]*(1-answers[4]),0,10,0,-95)-5);
        text((answers[0]*answers[3]).toFixed(1),ml.lx+26,ml.by+map(answers[0]*answers[3],0,10,0,-95)-5);
        text(t2.toFixed(1),ml.lx+96,ml.by+map(answers[0]*(1-answers[3]),0,10,0,-95)-5);
        pop();
    }
  }

}


function mathSolve(){
    const mixture = g.mixtures[mixtureSelect.value];
    if (!mixture) return;
    
    // Choose which solution set to use based on mixture
    let data = solutionsA; // default to solutions.json
    if (mixtureSelect.value === 'B') {
        data = solutionsB;  // use test.json for mixture B
    } else if (mixtureSelect.value === 'C') {
        data = solutionsC;  // use test2.json for mixture C
    } else if (mixtureSelect.value === 'D') {
        data = solutionsD;  // use test3.json for mixture D
    }
    
    // Calculate index and get values
    let ind = Math.round((g.T-120)*1616 + 100*g.moleFrac*16 + (g.P/.25 - 1));
    let t = data[ind];
    
    return([t[3],t[4],t[5],t[6],t[7]]);
}

// For methanol
function PsatM(T){
  let val = (1/750.06)*10**(8.08097 - 1582.271/(T+239.726));
  return(val);
}

// For water
function PsatW(T){
  let val = (1/750.06)*10**(8.07131 - 1730.63/(T+233.426));
  return(val);
}

function calculatePsat(T, antoine) {
    let val = (1/750.06) * 10**(antoine.A - antoine.B/(T + antoine.C));
    return val;
}

// For creating arrows
function arrow(base,tip,color,arrowLength,arrowWidth){ 
  // base = [x,y] tip = [x,y]
  // let arrowLength = 20; // Length of arrow
  // let arrowWidth = 5; // width of arrow (1/2)
  let dx, dy, mag;
  let u_hat, u_perp;
  let point = new Array(2); // Point along unit vector that is base of triangle
  let vert = new Array(6); // Holds vertices of arrow
  // Need to define a unit vector
  dx = tip[0] - base[0];
  dy = tip[1] - base[1];
  mag = (dx**2 + dy**2)**(1/2);
  u_hat = [dx/mag,dy/mag];
  vert[0] = tip[0] - 2*u_hat[0]; // Shifts the arrow back some to keep the tip from going out too far
  vert[1] = tip[1] - 2*u_hat[1];
  // Perpendicular unit vector
  u_perp = [-u_hat[1],u_hat[0]];
  // Base of arrow
  point[0] = vert[0]+ -arrowLength*u_hat[0];
  point[1] = vert[1]+ -arrowLength*u_hat[1];
  
  vert[2] = point[0] + u_perp[0]*arrowWidth;
  vert[3] = point[1] + u_perp[1]*arrowWidth;
  vert[4] = point[0] + -u_perp[0]*arrowWidth;
  vert[5] = point[1] + -u_perp[1]*arrowWidth;
  push();
  stroke(color); fill(color); strokeWeight(1);
  triangle(vert[0],vert[1],vert[2],vert[3],vert[4],vert[5]);
  pop();

}

let flashDrum = [[25.3231031543052, 52.80136402387042],[25.861892583120206, 49.02983802216539],[27.20886615515772, 45.52770673486786],[28.82523444160273, 42.56436487638534],[30.710997442455245, 40.13981244671782],[32.86615515771526, 37.7152600170503],[35.290707587382784, 35.56010230179028],[37.9846547314578, 33.135549872122766],[40.67860187553283, 31.24978687127025],[43.372549019607845, 29.63341858482523],[46.33589087809037, 28.017050298380223],[49.29923273657289, 26.670076726342714],[53.34015345268543, 25.0537084398977],[56.84228473998295, 23.976129582267692],[60.34441602728048, 22.62915601023018],[64.65473145780051, 21.55157715260017],[68.15686274509804, 20.743393009377666],[71.92838874680308, 19.665814151747657],[76.2387041773231, 19.127024722932653],[80.81841432225065, 18.85763000852515],[84.85933503836317, 18.58823529411765],[88.90025575447571, 18.049445865302644],[94.28815004262574, 17.78005115089514],[98.86786018755329, 17.78005115089514],[103.44757033248082, 17.78005115089514],[107.48849104859336, 18.04944586530264],[111.26001705029839, 18.58823529411765],[115.03154305200341, 19.12702472293265],[118.80306905370844, 19.66581415174765],[122.57459505541348, 20.20460358056266],[126.34612105711851, 20.47399829497016],[129.84825234441604, 21.55157715260017],[133.08098891730606, 22.62915601023018],[136.58312020460357, 23.70673486786018],[139.81585677749362, 24.51491901108269],[143.58738277919863, 26.13128729752771],[147.62830349531117, 28.01705029838022],[151.9386189258312, 29.902813299232736],[155.44075021312875, 32.32736572890025],[158.40409207161127, 34.75191815856777],[162.1756180733163, 37.4458653026428],[164.60017050298381, 40.67860187553283],[166.75532821824382, 43.91133844842285],[167.83290707587383, 46.87468030690537],[168.64109121909635, 50.3768115942029],[168.91048593350385, 56.84228473998295],[168.91048593350385, 269.1253196930946],[168.10230179028133, 272.8968456947996],[166.48593350383632, 276.9377664109122],[164.33077578857632, 279.9011082693947],[161.6368286445013, 282.86445012787726],[159.21227621483376, 285.5583972719523],[155.97953964194375, 287.9829497016198],[152.4774083546462, 290.1381074168798],[148.97527706734869, 292.5626598465473],[145.20375106564364, 294.4484228473998],[140.8934356351236, 296.06479113384484],[136.8525149190111, 297.6811594202899],[132.27280477408357, 299.0281329923273],[128.23188405797103, 300.1057118499574],[123.92156862745098, 301.1832907075874],[118.26427962489345, 301.9914748508099],[113.6845694799659, 302.7996589940324],[108.56606990622336, 303.3384484228474],[103.98635976129583, 303.8772378516624],[98.86786018755329, 303.8772378516624],[93.47996589940324, 303.6078431372549],[88.3614663256607, 303.6078431372549],[83.24296675191816, 303.3384484228474],[79.74083546462063, 303.0690537084399],[75.4305200341006, 302.26086956521743],[71.12020460358056, 301.9914748508099],[67.34867860187553, 300.9138959931799],[63.307757885763, 299.2975277067349],[58.99744245524297, 298.21994884910487],[55.22591645353794, 296.8729752770674],[51.993179880647915, 295.2566069906223],[47.952259164535384, 293.9096334185848],[44.450127877237854, 291.4850809889173],[40.947996589940324, 289.5993179880648],[37.7152600170503, 287.17476555839727],[33.94373401534527, 284.75021312872974],[31.788576300085253, 282.3256606990622],[30.17220801364024, 280.70929241261723],[28.82523444160273, 278.5541346973572],[27.47826086956522, 275.5907928388747],[26.13128729752771, 272.62745098039215],[25.592497868712705, 268.8559249786871],[25.592497868712705, 264.2762148337596],[25.592497868712705, 56.57289002557545]];