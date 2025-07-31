// Lightweight particle field + gentle hue shift
const canvas = document.getElementById("bg-canvas");
const ctx = canvas.getContext("2d");
let w, h, particles;

function resize(){
  w = canvas.width = window.innerWidth;
  h = canvas.height = window.innerHeight;
  particles = Array.from({length: Math.min(120, Math.floor(w*h/18000))}, () => ({
    x: Math.random()*w,
    y: Math.random()*h,
    vx: (Math.random()-0.5)*0.4,
    vy: (Math.random()-0.5)*0.4,
    r: 1.1 + Math.random()*1.9
  }));
}
window.addEventListener("resize", resize);
resize();

let hue = 210;
function tick(){
  hue = (hue + 0.02) % 360;
  ctx.clearRect(0,0,w,h);
  ctx.fillStyle = `hsla(${hue}, 70%, 60%, 0.06)`;
  ctx.fillRect(0,0,w,h);

  // draw links
  for(let i=0;i<particles.length;i++){
    for(let j=i+1;j<particles.length;j++){
      const a = particles[i], b = particles[j];
      const dx = a.x-b.x, dy = a.y-b.y, d2 = dx*dx+dy*dy;
      if (d2 < 140*140){
        const alpha = 0.35 * (1 - d2/(140*140));
        ctx.strokeStyle = `rgba(96,165,250,${alpha})`;
        ctx.beginPath(); ctx.moveTo(a.x,a.y); ctx.lineTo(b.x,b.y); ctx.stroke();
      }
    }
  }
  // draw particles
  for(const p of particles){
    p.x += p.vx; p.y += p.vy;
    if (p.x<0||p.x>w) p.vx*=-1;
    if (p.y<0||p.y>h) p.vy*=-1;
    ctx.beginPath();
    ctx.arc(p.x,p.y,p.r,0,Math.PI*2);
    ctx.fillStyle = "rgba(193, 219, 255, 0.6)";
    ctx.fill();
  }
  requestAnimationFrame(tick);
}
tick();
