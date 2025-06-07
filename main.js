const EPSILON = 1E-10; // precisione massima
const MAX_ITERATIONS = 1000; // numero massimo di iterazioni. Utile in caso di lenta convergenza
const SCALE_FACTOR = 100;

// inizializzatore per il metodo di Durand-Kerner
// non c'è nulla di speciale in questo numero
// eccetto per il fatto che non è un numero reale o una radice dell'unità
const initializer = math.complex(0.5, 0.7);

var funcInput;
var fdiv;
var showRootsCB;
var pauseCB;

var graphics;

let poly = "";

let hwratio = 0.75;
let w = 0;
let h = 0;

let ws = 0;
let hs = 0;

let coeff_values = [];
let roots = [];

let showRoots = false;

function setup() {
  canvas = select("#canvas");
  canvas.child(createCanvas(w = document.getElementById("canvas").clientWidth, h = w * hwratio));
  graphics = createGraphics(w, h);
  
  ws = 0.5*w/SCALE_FACTOR;
  hs = 0.5*h/SCALE_FACTOR;
  
  funcInput = createInput("1 0 0 1");
  select("#func").child(funcInput);
  btn = createButton('Aggiorna');
  select("#enter").child(btn);
  showRootsCB = createCheckbox("Mostra radici");
  select("#showRoots").child(showRootsCB);
  pauseCB = createCheckbox("Pausa");
  select("#showRoots").child(pauseCB);
  
  funcInput.input(inputChanged, false);
  btn.mousePressed(update);
  showRootsCB.changed(updateChecks);
  
  fdiv = createDiv("f(x)=");
  fdiv.position(10, 0);
  fdiv.style("font-size: large;");
  fdiv.style("text-shadow: 2px 0 0 #fff, -2px 0 0 #fff, 0 2px 0 #fff, 0 -2px 0 #fff, 1px 1px #fff, -1px -1px 0 #fff, 1px -1px 0 #fff, -1px 1px 0 #fff;");
  
  inputChanged();
  
  durandKerner();
  
  loop();
}

var it_i = 0, it_j = 0;
var px_width = 256;

function draw() {
  if(pauseCB.checked()) return;
  
  var t0 = millis();
  
  colorMode(HSB, 360, 1, 1);
  
  const segment = 360 / roots.length;
  
  if(it_i > w){
    if(px_width == 1){
      console.log("finished");
      noLoop();
      drawImageAndOverlay();
      return;
    }
    
    it_i = 0;
    px_width /= 2;
    
    console.log("drawing at pxw=" + px_width);
  }
  
  for(it_j = 0; it_j < h; it_j+=px_width){
    var re = it_i/SCALE_FACTOR;
    var im = it_j/SCALE_FACTOR;
    var z = math.complex(re - ws, im - hs);
    var root = newton(z);
    var index = 0;

    for(var k = 0; k < roots.length; k++){
      var delta = math.abs(math.subtract(root.result, roots[k]));
      if(delta < EPSILON){
        index = k;
        break;
      }
    }
    
    let c = (root.iterations == MAX_ITERATIONS) ? black : color(segment * index, 1, (root.iterations < 40) ? (root.iterations/40 + 0.5) : 1);
    graphics.stroke(c);
    graphics.fill(c);
    graphics.rect(it_i, it_j, px_width, px_width);
  }
  
  it_i+=px_width;
  
  drawImageAndOverlay();
  
  var t = millis() - t0;
}

function drawImageAndOverlay(){
  image(graphics, 0, 0);
  
  if(showRoots){
    stroke(0);
    fill(0);
    ellipseMode(CENTER);
    for(const root of roots){
      ellipse((root.re + ws) * SCALE_FACTOR, (root.im + hs) * SCALE_FACTOR, 10, 10);
    }
  }
}

/*
*  Questo metodo è utilizzato per trovare tutte le radici di una funzione
*  contemporaneamente. Utile per pre-calcolare le radici prima di applicare newton in ogni
*  punto.
*/
function durandKerner(){
  roots = [];
  
  const attractors = [];
  const moduli = [];
  const degree = coeff_values.length - 1;
  
  for(var i = 0; i < degree; i++)
    attractors[i] = fastpow(initializer, i);
  
  var delta = 100;
  var it = 0;
  while(delta > EPSILON && it < MAX_ITERATIONS){
    for(i = 0; i < degree; i++){
      var rho = 1;
      for(var j = i + 1; j < degree + i; j++){
        var index = j % degree;
        var d = math.subtract(attractors[i], attractors[index]);
        rho = math.multiply(rho, d);
      }
      
      var sigma = math.divide(f(attractors[i]), rho);
      attractors[i] = math.subtract(attractors[i], sigma);
      
      var delta0 = math.abs(attractors[i]);
      if(delta0 < delta)
        delta = delta0;
    }
    
    it++;
  }

  roots = attractors;
}

/*
*  Metodo di Newton-Raphson (tangenti) per trovare una delle radici a partire da un numero 
*  z0. Basato sulla formula di ricorrenza z_n=z_{n-1}-f(z_{n-1})/D[f(z_{n-1})].
*/
function newton(z0){
  var z = z0;
  var w = 1E16;
  var i = 0;
  
  while((math.abs(w) > EPSILON) && i < MAX_ITERATIONS){
    w = f(z);
    z = math.subtract(z, math.divide(w, df(z)));
    i++;
  }
  
  return new NewtonResult(z, i);
}

function f(z){
  var w = math.complex(0, 0);
  
  const degree = coeff_values.length;
  
  for(var i = 0; i < degree; i++){
    w = w.add(fastmulre(fastpow(z, degree-i-1), coeff_values[i]));
  }
  
  return w;
}

function df(z){
  var w = math.complex(0, 0);
  
  const len = coeff_values.length;
  
  for(var i = 0; i < len-1; i++){
    w = w.add(fastmulre(fastpow(z, len-i-2), coeff_values[i] * (len-i-1)));
  }
  
  return w;
}

class NewtonResult {
  constructor(result, iterations){
    this.result = result;
    this.iterations = iterations;
  }
}

// Eventi

function inputChanged(reset){
  noLoop();
  
  let s = funcInput.value();
  
  poly = "f(x)=";
  
  funcInput.style('color: #000');
  
  if(s == " "){
    funcInput.value('0');
    s = funcInput.value();
  }
  
  try{
    s = s.trim();

    const coefficients = s.split(" ");
    const len = coefficients.length;

    coeff_values = [];

    var i = -1;
    for(const coeff of coefficients){
      coeff_values[++i] = parseFloat(coeff);
      
      if(isNaN(coeff_values[i])) throw 'invalid string';
      
      if(coeff_values[i] != 0){
        var exp = len-i-1;
        var exps;

        switch(exp){
          case 0: 
            exps = "";
            break;
          case 1:
            exps = "x";
            break;
          default:
            exps = "x<sup>" + exp + "</sup>";
        }

        poly = poly.concat((coeff_values[i] < 0 || exp == len - 1)?"":"+", (coeff_values[i]==1 && exp!=0)?"":coeff_values[i], exps);
      }
    }
    
    if(coeff_values[0] == 0 && len == 1) poly = poly.concat('0');
    
    for(i = coeff_values.length-1; i >= 0; i--)
      coeff_values[i] /= coeff_values[0];
  
  } catch(error) {
    funcInput.style('color: #f00');
    
    if(!reset){
      poly = 'input non valido.';
    } else {
      //funcInput.value('0');
      inputChanged(false);
    }
  }
  
  fdiv.html(poly);
}

function update(){
  noLoop();
  
  let s = funcInput.value();
  
  if(s == "") funcInput.value('0');
  
  durandKerner();
  inputChanged(true);
  
  background(255);
  
  it_i = 0;
  px_width = 256;
  
  loop();
}

function updateChecks(){
  showRoots = showRootsCB.checked();
  
  draw();
}
