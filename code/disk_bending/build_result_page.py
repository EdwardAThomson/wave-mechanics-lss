#!/usr/bin/env python3
"""Build the phase-spiral result page: inlines figures/*.png as data URIs.

Source of truth is output/phase_spiral_result.md and figures/*.png, both
committed. The generated HTML is not committed (1.4 MB of inlined images);
run this to regenerate it.
"""
import base64, pathlib

BASE = pathlib.Path("/home/user/wave-mechanics-lss/code/disk_bending/figures")
OUT = pathlib.Path("/home/user/wave-mechanics-lss/code/disk_bending/output/phase_spiral_result.html")


def uri(name):
    b = (BASE / name).read_bytes()
    return "data:image/png;base64," + base64.b64encode(b).decode()


HTML = """<title>The Vertical Phase Spiral</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
:root {
  --ground:  #EFF2F6;
  --surface: #FFFFFF;
  --sunk:    #E4E9F0;
  --ink:     #10151C;
  --ink-2:   #5A6672;
  --ink-3:   #8894A3;
  --accent:  #23478A;
  --warm:    #A03028;
  --rule:    #D3DAE3;
  --figure-adjust: none;

  --serif: "Iowan Old Style", "Charter", "Bitstream Charter", "Palatino Linotype", Palatino, Georgia, serif;
  --sans: ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
  --mono: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, "Liberation Mono", monospace;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --ground:  #0D1117;
    --surface: #161C24;
    --sunk:    #11171E;
    --ink:     #E6EBF2;
    --ink-2:   #9BA7B6;
    --ink-3:   #6E7B8B;
    --accent:  #7FA6E8;
    --warm:    #E29189;
    --rule:    #263140;
    --figure-adjust: brightness(.92);
  }
}
:root[data-theme="dark"] {
  --ground:  #0D1117;
  --surface: #161C24;
  --sunk:    #11171E;
  --ink:     #E6EBF2;
  --ink-2:   #9BA7B6;
  --ink-3:   #6E7B8B;
  --accent:  #7FA6E8;
  --warm:    #E29189;
  --rule:    #263140;
  --figure-adjust: brightness(.92);
}

* { box-sizing: border-box; }
body {
  margin: 0;
  background: var(--ground);
  color: var(--ink);
  font-family: var(--sans);
  font-size: 17px;
  line-height: 1.65;
  -webkit-font-smoothing: antialiased;
}
.wrap { max-width: 1080px; margin: 0 auto; padding: 0 24px 96px; }
.col  { max-width: 68ch; }

/* ---- masthead ---- */
header {
  border-bottom: 1px solid var(--rule);
  padding: 72px 0 40px;
  margin-bottom: 56px;
}
.eyebrow {
  font-family: var(--mono);
  font-size: 12px;
  letter-spacing: .13em;
  text-transform: uppercase;
  color: var(--accent);
  margin: 0 0 20px;
}
h1 {
  font-family: var(--serif);
  font-weight: 600;
  font-size: clamp(2.2rem, 5.2vw, 3.4rem);
  line-height: 1.08;
  letter-spacing: -.015em;
  margin: 0 0 20px;
  text-wrap: balance;
}
.standfirst {
  font-family: var(--serif);
  font-size: clamp(1.05rem, 2.2vw, 1.3rem);
  line-height: 1.5;
  color: var(--ink-2);
  margin: 0 0 32px;
  max-width: 62ch;
}
.meta {
  display: flex; flex-wrap: wrap; gap: 8px 28px;
  font-family: var(--mono); font-size: 12.5px; color: var(--ink-3);
}
.meta b { color: var(--ink-2); font-weight: 500; }

/* ---- type ---- */
h2 {
  font-family: var(--serif);
  font-weight: 600;
  font-size: 1.72rem;
  letter-spacing: -.01em;
  line-height: 1.2;
  margin: 72px 0 8px;
  text-wrap: balance;
}
h3 {
  font-family: var(--sans);
  font-weight: 650;
  font-size: 1.02rem;
  letter-spacing: .005em;
  margin: 36px 0 6px;
}
p { margin: 0 0 18px; }
a { color: var(--accent); text-decoration-thickness: 1px; text-underline-offset: 2px; }
a:focus-visible, summary:focus-visible {
  outline: 2px solid var(--accent); outline-offset: 3px; border-radius: 2px;
}
code {
  font-family: var(--mono); font-size: .88em;
  background: var(--sunk); padding: .12em .34em; border-radius: 3px;
}
strong { font-weight: 650; }
.lede { color: var(--ink-2); }

/* ---- findings: a real sequence, ordered by how surprising they were ---- */
.findings { display: grid; gap: 2px; margin: 32px 0 8px; background: var(--rule);
            border: 1px solid var(--rule); border-radius: 4px; overflow: hidden; }
.finding { background: var(--surface); padding: 22px 26px; display: grid;
           grid-template-columns: 42px 1fr; gap: 0 18px; align-items: start; }
.finding .n {
  font-family: var(--mono); font-size: 12px; color: var(--accent);
  padding-top: 4px; letter-spacing: .06em;
}
.finding h4 { margin: 0 0 4px; font-size: 1.02rem; font-weight: 650; letter-spacing: .005em; }
.finding p { margin: 0; font-size: .95rem; color: var(--ink-2); }

/* ---- figures ---- */
figure { margin: 40px 0; }
.plate {
  background: #fff; border: 1px solid var(--rule); border-radius: 4px;
  padding: 10px; overflow-x: auto;
}
.plate img { display: block; width: 100%; height: auto; min-width: 640px;
             filter: var(--figure-adjust); }
figcaption {
  font-size: .875rem; color: var(--ink-2); margin-top: 14px;
  max-width: 72ch; line-height: 1.55;
}
figcaption b { color: var(--ink); font-weight: 650; }

/* ---- tables ---- */
.tw { overflow-x: auto; margin: 26px 0; }
table { border-collapse: collapse; width: 100%; font-size: .92rem; min-width: 460px; }
th, td { text-align: left; padding: 11px 16px 11px 0; border-bottom: 1px solid var(--rule); }
thead th {
  font-family: var(--mono); font-size: 11.5px; letter-spacing: .09em;
  text-transform: uppercase; color: var(--ink-3); font-weight: 500;
  border-bottom: 1px solid var(--ink-3);
}
td.num, th.num { text-align: right; font-family: var(--mono);
                 font-variant-numeric: tabular-nums; padding-right: 0; }
tbody tr:last-child td { border-bottom: none; }
.hi { color: var(--accent); font-weight: 650; }
.hi-warm { color: var(--warm); font-weight: 650; }

/* ---- callout ---- */
.note {
  border-left: 2px solid var(--accent); padding: 2px 0 2px 20px;
  margin: 28px 0; color: var(--ink-2); font-size: .95rem;
}
.note.warm { border-left-color: var(--warm); }
.note p:last-child { margin-bottom: 0; }

pre {
  background: var(--sunk); border: 1px solid var(--rule); border-radius: 4px;
  padding: 18px 20px; overflow-x: auto; font-family: var(--mono);
  font-size: 13px; line-height: 1.6; margin: 24px 0;
}
pre code { background: none; padding: 0; font-size: inherit; }

ul { margin: 0 0 18px; padding-left: 22px; }
li { margin-bottom: 9px; }

footer {
  margin-top: 88px; padding-top: 28px; border-top: 1px solid var(--rule);
  font-family: var(--mono); font-size: 12.5px; color: var(--ink-3);
}
@media (prefers-reduced-motion: reduce) {
  * { animation: none !important; transition: none !important; }
}
</style>

<div class="wrap">
<header>
  <p class="eyebrow">wave-mechanics-lss &middot; disk bending &middot; stage 0</p>
  <h1>The Vertical Phase Spiral</h1>
  <p class="standfirst">A Schr&ouml;dinger&ndash;Poisson solver, built to simulate cosmological
  large-scale structure, pointed instead at a kicked galactic disk. The Gaia snail comes out
  with no particle noise at all &mdash; and disk self-gravity turns out to suppress it.</p>
  <div class="meta">
    <span><b>Solver</b> 1D vertical slab, warm, self-gravitating</span>
    <span><b>Grid</b> 4096 &times; 8 kpc, 320 streams</span>
    <span><b>Validation</b> 6/6 tests passing</span>
  </div>
</header>

<div class="col">
<p class="lede">Schr&ouml;dinger&ndash;Poisson is a Vlasov&ndash;Poisson solver in disguise.
Husimi-smoothing the wavefunction recovers the distribution function to
<code>O(&hbar;<sub>eff</sub>)</code>, so anything collisionless is fair game &mdash; including
vertical phase mixing in a disk. The specific pay-off is shot-noise-free phase space, which
is exactly where an N-body code has to work hardest.</p>

<p>Three findings, ordered by how surprising they were.</p>
</div>

<div class="findings">
  <div class="finding"><div class="n">01</div><div>
    <h4>The spiral forms cleanly, with no shot noise</h4>
    <p>The distribution function comes out of a single grid calculation. No particles, no
    sampling, no post-hoc smoothing kernel &mdash; contrast is limited only by the
    phase-space resolution.</p>
  </div></div>
  <div class="finding"><div class="n">02</div><div>
    <h4>Strong disk self-gravity suppresses it</h4>
    <p>Two runs identical but for surface density retain 0.4% and 65% of the corrugation after
    14 vertical periods. The self-gravitating slab sloshes as a collective mode instead of
    winding.</p>
  </div></div>
  <div class="finding"><div class="n">03</div><div>
    <h4>A harmonic vertical potential gives no spiral at all</h4>
    <p>By Kohn&rsquo;s theorem the centre-of-mass mode decouples exactly and never damps. The
    winding comes entirely from anharmonicity &mdash; an easy way to lose a week with nothing
    visibly wrong in the output.</p>
  </div></div>
</div>

<div class="col">
<h2>The spiral</h2>
<p>A warm, self-gravitating slab held by a rigid halo term is given a 15 km/s impulsive
vertical kick and evolved for 14 vertical periods. Each panel is the Husimi transform divided
by the <em>measured</em> pre-kick equilibrium, so the smooth background divides out exactly and
no circular symmetry is assumed. The vertical axis is <code>v<sub>z</sub>/&omega;<sub>z</sub></code>
in kpc, so an unperturbed orbit is a circle and a wound spiral looks like one.</p>
</div>

<figure>
  <div class="plate"><img src="__TRACER__" alt="Nine panels showing the vertical phase-space
  contrast winding from a single lobe at t=0 into a tight multi-armed spiral by t=1578 Myr."></div>
  <figcaption><b>A tracer population winds up into a multi-armed snail.</b>
  &Sigma; = 5 M<sub>&#9737;</sub>/pc&sup2;, 320 streams on a 4096 grid. The corrugation amplitude
  decays from 0.247 kpc to 0.0009 kpc, a factor of 280 &mdash; that decay <em>is</em> the winding.
  The perturbation is not lost, it is wrapped into ever finer phase-space structure until the
  mean displacement averages to nothing.</figcaption>
</figure>

<div class="col">
<h2>Self-gravity suppresses it</h2>
<p>Identical setup, ten times the stellar surface density. The sheet barely winds: its dipole
response is a collective normal mode, and because self-gravity is an internal force the whole
structure moves together rather than each orbit phase-mixing at its own frequency.</p>
</div>

<figure>
  <div class="plate"><img src="__SELFGRAV__" alt="Nine panels showing a single m=1 lobe that
  rotates but does not wind into a spiral over the same duration."></div>
  <figcaption><b>The self-gravitating case does not wind.</b>
  &Sigma; = 50 M<sub>&#9737;</sub>/pc&sup2;, everything else unchanged. The external potential is
  distinctly anharmonic in both runs and the orbital frequency varies by roughly 40% across the
  occupied range, so this is not for lack of anharmonicity.</figcaption>
</figure>

<div class="col">
<div class="tw">
<table>
  <thead><tr>
    <th>Run</th><th class="num">&Sigma; [M<sub>&#9737;</sub>/pc&sup2;]</th>
    <th class="num">h [kpc]</th><th class="num">T<sub>z</sub> [Myr]</th>
    <th class="num">Corrugation, early &rarr; late</th><th class="num">Retained</th>
  </tr></thead>
  <tbody>
    <tr><td>Tracer</td><td class="num">5</td><td class="num">0.229</td><td class="num">98.4</td>
        <td class="num">0.247 &rarr; 0.0009</td><td class="num hi">0.4%</td></tr>
    <tr><td>Self-gravitating</td><td class="num">50</td><td class="num">0.155</td><td class="num">68.8</td>
        <td class="num">0.277 &rarr; 0.180</td><td class="num hi-warm">65%</td></tr>
  </tbody>
</table>
</div>

<p>Both runs conserve energy to better than 5 parts in 10<sup>6</sup> and hold velocity-space
spill below 4 parts in 10<sup>11</sup>, so neither result is an aliasing artefact.</p>

<div class="note">
<p>The practical consequence: the observed sharpness of the Milky Way&rsquo;s phase spiral carries
information about the disk-to-halo ratio of the vertical restoring force. A simulation that gets
that ratio wrong will get the winding rate wrong even with the kick and the potential shape
otherwise right.</p>
</div>
</div>

<figure>
  <div class="plate"><img src="__MOMENTS__" alt="Three stacked panels: corrugation and mean
  vertical velocity in quadrature, oscillation envelope decaying for the tracer run and
  plateauing for the self-gravitating one, and conservation diagnostics."></div>
  <figcaption><b>The same contrast, quantitatively.</b> The middle panel is the oscillation
  envelope: the tracer run (blue) decays away, the self-gravitating run (orange) plateaus. The
  bottom panel is the conservation and aliasing check for both.</figcaption>
</figure>

<div class="col">
<h2>Validation</h2>
<p>Six tests, run in order, each catching a different class of error &mdash; about four minutes
end to end.</p>

<div class="tw">
<table>
  <thead><tr><th class="num">#</th><th>Test</th><th>Result</th></tr></thead>
  <tbody>
    <tr><td class="num">0</td><td>Isolated vertical Poisson vs analytic sech&sup2; sheet</td>
        <td>2nd order; edge force exact to 3&times;10<sup>&minus;9</sup></td></tr>
    <tr><td class="num">1</td><td>Free particle: plane-wave phase, Gaussian spreading</td>
        <td>10<sup>&minus;14</sup></td></tr>
    <tr><td class="num">2</td><td>Harmonic oscillator: spectrum and coherent state</td>
        <td>eigen-residual 10<sup>&minus;12</sup>; Strang 2nd order</td></tr>
    <tr><td class="num">3</td><td>Jeans growth and the quantum-pressure branch</td>
        <td>3&times;10<sup>&minus;9</sup> and 10<sup>&minus;11</sup></td></tr>
    <tr><td class="num">4</td><td>Warm isothermal sheet stationarity</td>
        <td>z<sub>rms</sub> drift 4&times;10<sup>&minus;6</sup> over 8 periods</td></tr>
    <tr><td class="num">5a</td><td>Hybrid (x, z) Poisson, four ways</td>
        <td>round-off (4&times;10<sup>&minus;16</sup>) end to end</td></tr>
  </tbody>
</table>
</div>

<p>Test 3 runs both signs of the dispersion relation, which pins the sign of the gravity coupling
and the kinetic operator together. Test 4 was expected to be the one that failed first; the
initial-condition construction below makes it pass by construction.</p>

<h2>Three method notes that mattered</h2>

<h3>Warm initial conditions come from eigenstates, not sampled streams</h3>
<p>A single <code>&psi; = &radic;&rho; e<sup>iS/&hbar;</sup></code> is a cold stream, and a disk
with &sigma;<sub>z</sub> = 20 km/s is not. The obvious fix &mdash; sampling streams from
<code>f(z, v<sub>z</sub>)</code> &mdash; has a trap: a constant-velocity stream is not stationary
in the sheet potential, it simply falls, and the streams that <em>are</em> stationary are
constant-energy tori needing two velocity branches with a density that diverges at the turning
point. Using eigenstates of the self-consistent vertical Hamiltonian instead makes
<code>|u<sub>n</sub>|&sup2;</code> exactly stationary by construction. Measured: the scale height
matches the analytic Spitzer value to five figures.</p>

<h3>The resolution constraint is the thesis&rsquo;s own headline lesson, reappearing</h3>
<p>The grid can only represent velocities up to
<code>v<sub>max</sub> = &hbar;<sub>eff</sub>&pi;/dz</code>. Refining
<code>&hbar;<sub>eff</sub></code> without refining <code>dz</code> makes velocity-space
diagnostics <em>worse</em>, not better &mdash; the same trap that turned out to explain the
&ldquo;quantum interference&rdquo; noise in the original 2011 cosmological runs, now in a
completely different problem.</p>

<h3>Extra streams buy variance reduction, not warmth</h3>
<p>A single wavefunction already carries the full warm distribution function. What extra streams
reduce is interference speckle, and the spurious vertical dipole that comes with it follows
<code>A &asymp; 0.70 &radic;(&hbar;<sub>eff</sub>/&omega;<sub>z</sub>) / &radic;N<sub>streams</sub></code>,
holding to a few per cent over a 64&times; range. Since the corrugation is the observable, that
noise competes directly with the signal.</p>

<h2>What is not claimed</h2>
<div class="note warm">
<ul>
  <li>This is a <b>1D vertical slab</b>. No in-plane wavevector means no bending waves, no
  bending dispersion relation and no firehose.</li>
  <li>The perturber is idealised &mdash; impulsive kicks, not a Sagittarius orbit.</li>
  <li><code>&hbar;<sub>eff</sub></code> = 0.6 kpc km/s gives a de Broglie length of 0.64 scale
  heights, so the finest arms of a real spiral would be smoothed at that scale.</li>
  <li>No comparison against an N-body run of the same setup has been made. The claim is that the
  spiral is noise-free, not that it is quantitatively closer to truth than a well-resolved
  particle calculation.</li>
</ul>
</div>

<h2>Reproducing</h2>
<pre><code>cd code/disk_bending
make &amp;&amp; make check                       # validation ladder, ~4 min

./sheet_1d --N 4096 --L 8 --hbar 0.6 --sigma 5 \\
           --nu-ext 60 --z-ext 0.5 --vkick 15 \\
           --ic multistream --periods 14 --out output/spiral_tracer

python3 plot_phase_spiral.py output/spiral_tracer</code></pre>
<p>Each production run is a few minutes on four cores. Add <code>--z-ext 50</code> to recover the
harmonic case and watch the spiral vanish.</p>
</div>

<footer>
  wave-mechanics-lss &middot; code/disk_bending &middot; Stage 0 &middot;
  full method and plan review in the directory README
</footer>
</div>
"""

HTML = (HTML.replace("__TRACER__", uri("phase_spiral_tracer.png"))
            .replace("__SELFGRAV__", uri("phase_spiral_selfgravitating.png"))
            .replace("__MOMENTS__", uri("moments_tracer_vs_selfgrav.png")))
OUT.write_text(HTML)
print(f"wrote {OUT}  ({OUT.stat().st_size/1e6:.2f} MB)")
