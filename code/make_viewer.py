#!/usr/bin/env python3
"""Generate a self-contained web viewer for the zonated APAP model.

Precomputes the per-hepatocyte time courses for each species, gradient scheme
(assumed vs measured) and dose, and writes a single HTML file (Plotly from CDN)
that renders a position x time heatmap with a time slider and a scheme toggle.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from zonated_apap_model import run, N_HEPATOCYTES  # noqa: E402

BASE = Path(__file__).resolve().parent.parent
SPECIES = ["P", "S", "N", "G", "C"]
SPECIES_LABEL = {"P": "APAP", "S": "PAPS / APAP-S", "N": "NAPQI", "G": "Glutathione", "C": "Protein adducts"}
N_TIMES = 120


def build_data():
    data = {}
    for scheme in ("assumed", "measured"):
        data[scheme] = {}
        for dose in (4.0, 16.0):
            t, Y, _ = run(scheme, dose, t_end_days=3.0, n_points=400, rtol=1e-5, atol=1e-8)
            idx = np.linspace(0, len(t) - 1, N_TIMES).astype(int)
            times = t[idx].round(3).tolist()
            species = {}
            for si, name in enumerate(SPECIES):
                species[name] = Y[idx, :, si].round(6).tolist()  # (time, position)
            data[scheme][str(int(dose))] = {"times": times, **species}
    return data


HTML = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>APAP zonation viewer</title>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<style>
  :root { --ink:#0b0b0b; --ink2:#52514e; --muted:#898781; --grid:#e1e0d9; --blue:#2a78d6; --orange:#eb6834; --surface:#fcfcfb; }
  * { box-sizing:border-box; }
  body { margin:0; font-family:system-ui,-apple-system,"Segoe UI",sans-serif; color:var(--ink); background:var(--surface); }
  header { padding:18px 24px 8px; border-bottom:1px solid var(--grid); }
  h1 { font-size:18px; margin:0 0 4px; }
  header p { color:var(--ink2); margin:0; font-size:13px; max-width:760px; }
  .controls { display:flex; flex-wrap:wrap; gap:18px; align-items:center; padding:12px 24px; }
  .group { display:flex; align-items:center; gap:6px; }
  .group label { font-size:12px; color:var(--muted); margin-right:2px; }
  .btn { font-size:12px; padding:5px 11px; border:1px solid var(--grid); background:#fff; color:var(--ink2); border-radius:5px; cursor:pointer; }
  .btn.on { background:var(--blue); color:#fff; border-color:var(--blue); }
  .slider { width:230px; }
  .play { background:var(--orange); color:#fff; border:none; padding:6px 14px; border-radius:5px; cursor:pointer; font-size:12px; }
  #plot { height:calc(100vh - 130px); min-height:520px; }
  .note { font-size:11px; color:var(--muted); padding:0 24px 12px; }
</style>
</head>
<body>
<header>
  <h1>Zonated acetaminophen metabolism along the sinusoid</h1>
  <p>Each of 16 hepatocytes (periportal → pericentral) runs the same intracellular
  kinetics with enzyme rates modulated by a zonal gradient. The heatmap shows a
  species' concentration over <em>position</em> (x) and <em>time</em> (y); the profile
  (right) is the snapshot at the slider time. Toggle <em>assumed</em> (literature-derived)
  vs <em>measured</em> (2026 human proteomics) gradients.</p>
</header>

<div class="controls">
  <div class="group" id="species-group">
    <label>Species</label>
    <button class="btn" data-sp="P">APAP</button>
    <button class="btn" data-sp="S">PAPS</button>
    <button class="btn on" data-sp="N">NAPQI</button>
    <button class="btn" data-sp="G">GSH</button>
    <button class="btn" data-sp="C">Adducts</button>
  </div>
  <div class="group" id="scheme-group">
    <label>Gradients</label>
    <button class="btn" data-sch="assumed">Assumed</button>
    <button class="btn on" data-sch="measured">Measured</button>
  </div>
  <div class="group" id="dose-group">
    <label>Dose</label>
    <button class="btn" data-dose="4">4 g</button>
    <button class="btn on" data-dose="16">16 g</button>
  </div>
  <div class="group">
    <button class="play" id="play">Play</button>
    <input type="range" class="slider" id="time-slider" min="0" max="1" step="0.001" value="0.4"/>
  </div>
</div>

<div id="plot"></div>
<div class="note">Concentrations are nmol per hepatocyte; the adduct compartment (C) is a
proxy for injury, not necrosis itself. No intercellular transport is modelled.</div>

<script>
const DATA = __DATA__;
const POSITIONS = Array.from({length:16}, (_,i)=>i+1);
let species = 'N', scheme = 'measured', dose = '16';

const el = id => document.getElementById(id);

function timesOf(){ return DATA[scheme][dose].times; }
function valuesOf(sp){ return DATA[scheme][dose][sp]; }   // (time, position)
function nearestTimeIdx(t){ const ts = timesOf(); let lo=0, hi=ts.length-1;
  while(lo<hi){ const m=(lo+hi)>>1; if(ts[m]<t) lo=m+1; else hi=m; } return lo; }

function buildTraces(){
  const Z = valuesOf(species);          // (time, position)
  const times = timesOf();
  const tIdx = nearestTimeIdx(timeSlider());
  const profile = Z[tIdx];
  const heat = { type:'heatmap', x:POSITIONS, y:times, z:Z, zsmooth:'best',
    colorscale:[[0,'#fcfcfb'],[1,'#2a78d6']], colorbar:{title:{text:'nmol',font:{size:10}},thickness:12},
    hovertemplate:'cell %{x}<br>t=%{y:.2f} d<br>%{z:.3g}<extra></extra>',
    xaxis:'x', yaxis:'y' };
  const line = { type:'scatter', x:POSITIONS, y:profile, mode:'lines+markers',
    line:{color:'#eb6834', width:2.5}, marker:{size:6, color:'#eb6834'},
    hovertemplate:'cell %{x}<br>%{y:.3g}<extra></extra>', xaxis:'x2', yaxis:'y2' };
  return { heat, line };
}

function layout(){
  const ts = timesOf();
  const tNow = timeSlider();
  return {
    grid:{rows:1, columns:2, pattern:'independent'},
    xaxis:{title:{text:'hepatocyte (1 = periportal, 16 = pericentral)', font:{size:11}}, tickvals:[1,4,8,12,16], dtick:0},
    yaxis:{title:{text:'time (days)', font:{size:11}}},
    xaxis2:{title:{text:'hepatocyte', font:{size:11}}, tickvals:[1,4,8,12,16]},
    yaxis2:{title:{text:'nmol', font:{size:11}}},
    shapes:[{type:'line', xref:'paper', yref:'y', x0:0, x1:0.485, y0:tNow, y1:tNow,
             line:{color:'#0b0b0b', width:1, dash:'dot'}}],
    annotations:[{text:'position × time', x:0.24, y:1.0, xref:'paper', yref:'paper', showarrow:false, font:{size:11, color:'#52514e'}},
                 {text:'snapshot at slider time', x:0.75, y:1.0, xref:'paper', yref:'paper', showarrow:false, font:{size:11, color:'#52514e'}}],
    margin:{t:34,l:52,r:16,b:44},
    paper_bgcolor:'#fcfcfb', plot_bgcolor:'#fcfcfb',
    font:{family:'system-ui, "Segoe UI", sans-serif', color:'#0b0b0b'},
    width:null, height:null
  };
}

function timeSlider(){ return parseFloat(el('time-slider').value); }

function render(){
  const {heat, line} = buildTraces();
  Plotly.react(el('plot'), [heat, line], layout(), {responsive:true});
}

function setSliderFromRange(){ const ts = timesOf(); el('time-slider').max = ts[ts.length-1];
  el('time-slider').min = ts[0]; }

function wireButtons(grpId, getter, setter, apply){
  el(grpId).addEventListener('click', e=>{
    const b = e.target.closest('.btn'); if(!b) return;
    el(grpId).querySelectorAll('.btn').forEach(x=>x.classList.remove('on'));
    b.classList.add('on');
    apply(b.dataset);
  });
}

wireButtons('species-group', null, null, d=>{ species = d.sp; render(); });
wireButtons('scheme-group', null, null, d=>{ scheme = d.sch; setSliderFromRange(); render(); });
wireButtons('dose-group',  null, null, d=>{ dose  = d.dose; setSliderFromRange(); render(); });

el('time-slider').addEventListener('input', render);

let playing=false, timer=null;
el('play').addEventListener('click', ()=>{
  playing=!playing;
  el('play').textContent = playing ? 'Pause' : 'Play';
  if(playing){
    timer=setInterval(()=>{
      const s=el('time-slider'); let v=parseFloat(s.value)+0.05;
      if(v>parseFloat(s.max)) v=parseFloat(s.min);
      s.value=v; render();
    },90);
  } else { clearInterval(timer); }
});

setSliderFromRange();
el('time-slider').value = timesOf()[Math.floor(timesOf().length*0.4)];
render();
</script>
</body>
</html>
"""


def main():
    data = build_data()
    out = BASE / "viewer.html"
    out.write_text(HTML.replace("__DATA__", json.dumps(data)), encoding="utf-8")
    print(f"wrote {out.resolve()}  ({out.stat().st_size/1e3:.0f} kB)")


if __name__ == "__main__":
    main()
