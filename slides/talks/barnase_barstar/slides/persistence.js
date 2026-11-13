// Every contact as a barcode over the simulation: a row per contact, inked where it is present. Drawn on a
// canvas (47 rows x 200 frames is a lot of rectangles for a vector scene) and animated here.
//
// Steps: 1 the rows are written left to right as the clock runs, most persistent on top · 2 the rows present
// in at least 95% of frames grow and get their names; the rest drop away · 3 of those, the ones that go
// through barstar Asp39.
{
const KEEP = 0.95;
const ease = (x) => (x < 0.5 ? 4 * x * x * x : 1 - (-2 * x + 2) ** 3 / 2);

Deck.slide('persistence', {
  steps: 3,
  notes: 'One row per contact (a type and a residue pair), one column per frame, 200 frames over 20 ns. Rows are '
    + 'sorted by occupancy. Most of the interface flickers: contacts form and break all the time. A small set is '
    + 'there in at least 95% of frames, and more than half of that set involves a single barstar residue, Asp39.',
  mount(frame, { palette: P, data, esc }) {
    const D = data.interface;
    const rows = D.contacts.map((c) => ({ ...c, on: [...c.present].map((ch) => ch === '1') }));
    const kept = rows.filter((r) => r.occupancy >= KEEP);
    const hub = kept.filter((r) => r.b === 'Asp39');
    frame.innerHTML = `
      <h2 class="deck-heading reveal-in in">${esc('Who holds on?')}</h2>
      <div class="ps-fig"><canvas></canvas></div>
      ${Deck.points([
        { title: `${rows.length} contacts, ${D.frames} frames`, sub: 'a row per contact, inked when it is there; most flicker', from: 1 },
        { title: `${kept.length} never let go`, sub: `present in at least ${Math.round(100 * KEEP)}% of the simulation`, from: 2 },
        { title: `${hub.length} of them meet at one residue`, sub: '<span class="ps-hub">barstar Asp39</span>', html: true, from: 3 },
      ], { className: 'ps-col' })}`;
    Deck.css(`.ps-fig { position: absolute; left: 4.2%; top: 17%; width: 62%; height: 67%; }
      .ps-fig canvas { display: block; width: 100%; height: 100%; }
      .ps-col { --col-left: 70%; --col-top: 24%; --col-gap: calc(2.6 * var(--u)); }
      .ps-hub { color: var(--plot-5); font-weight: 700; font-size: var(--t-lead); }`);

    const canvas = frame.querySelector('canvas'), ctx = canvas.getContext('2d');
    const S = { sweep: 0, focus: 0, hub: 0 };
    let gen = 0, dead = false;

    function draw() {
      const w = canvas.clientWidth, h = canvas.clientHeight, dpr = window.devicePixelRatio || 1;
      if (!w) return;
      if (canvas.width !== Math.round(w * dpr)) { canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr); }
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      ctx.clearRect(0, 0, w, h);
      const u = frame.clientWidth / 100;
      const label = 23 * u * S.focus;                    // room for names once the kept rows are named
      const x0 = label, x1 = w - 0.5 * u, top = 0, bottom = h - 3 * u;
      const all = (bottom - top) / rows.length, big = (bottom - top) / kept.length;
      const n = D.frames, bw = (x1 - x0) / n, upto = S.sweep * n;
      let k = 0;
      rows.forEach((r, i) => {
        const keep = r.occupancy >= KEEP, isHub = keep && r.b === 'Asp39';
        // kept rows grow to share the height; the others slide down and fade
        const y = keep ? top + (i * all) * (1 - S.focus) + (k * big) * S.focus : top + i * all + S.focus * 6 * u;
        const rh = keep ? all + (big - all) * S.focus : all;
        const alpha = keep ? 1 : 1 - S.focus;
        if (keep) k++;
        if (alpha <= 0.01) return;
        ctx.globalAlpha = alpha * (keep && !isHub ? 1 - 0.72 * S.hub : 1);
        ctx.fillStyle = isHub && S.hub > 0 ? P.series[4] : (keep && S.focus > 0.5 ? P.fg : P.series[0]);
        const pad = Math.min(0.14 * rh, 0.5 * u);
        for (let f = 0; f < n && f < upto; f++) {
          if (!r.on[f]) continue;
          let g = f;
          while (g + 1 < n && g + 1 < upto && r.on[g + 1]) g++;
          ctx.fillRect(x0 + f * bw, y + pad, (g - f + 1) * bw + 0.5, rh - 2 * pad);
          f = g;
        }
        if (keep && S.focus > 0.02) {
          ctx.globalAlpha = S.focus * (isHub ? 1 : 1 - 0.6 * S.hub);
          ctx.fillStyle = isHub && S.hub > 0.5 ? P.series[4] : P.fg;
          ctx.font = `600 ${1.5 * u}px Lato, sans-serif`;
          ctx.textBaseline = 'middle';
          ctx.textAlign = 'right';
          ctx.fillText(`${r.a} · ${r.b}`, label - 8.6 * u, y + rh / 2);
          ctx.font = `400 ${1.15 * u}px Lato, sans-serif`;
          ctx.fillStyle = P.muted;
          ctx.textAlign = 'left';
          ctx.fillText(r.type.replace(' / ionic', '').toLowerCase(), label - 8 * u, y + rh / 2);
        }
      });
      // the time axis, and the clock while the rows are being written
      ctx.globalAlpha = 1;
      ctx.fillStyle = P.muted;
      ctx.font = `400 ${1.4 * u}px Lato, sans-serif`;
      ctx.textBaseline = 'top';
      ctx.textAlign = 'left';
      ctx.fillText('0', x0, bottom + 0.9 * u);
      ctx.textAlign = 'right';
      ctx.fillText(`${D.ns.toFixed(0)} ns`, x1, bottom + 0.9 * u);
      if (S.sweep > 0 && S.sweep < 1) {
        ctx.strokeStyle = P.accent;
        ctx.lineWidth = 0.2 * u;
        ctx.beginPath(); ctx.moveTo(x0 + upto * bw, top); ctx.lineTo(x0 + upto * bw, bottom); ctx.stroke();
      }
    }

    const TARGET = [{ sweep: 1, focus: 0, hub: 0 }, { sweep: 1, focus: 1, hub: 0 }, { sweep: 1, focus: 1, hub: 1 }];
    function go(step, animate) {
      const id = ++gen, to = TARGET[step - 1];
      if (animate !== true) { Object.assign(S, to); draw(); return undefined; }
      const from = step === 1 ? { sweep: 0, focus: 0, hub: 0 } : { ...S }, ms = step === 1 ? 3200 : 1300, t0 = performance.now();
      return new Promise((done) => {
        const tick = (t) => {
          if (dead || id !== gen) return done();
          const x = Math.min(1, (t - t0) / ms), e = step === 1 ? x : ease(x);
          for (const key of Object.keys(S)) S[key] = from[key] + (to[key] - from[key]) * e;
          draw();
          if (x < 1) requestAnimationFrame(tick); else done();
        };
        requestAnimationFrame(tick);
      });
    }
    const resize = new ResizeObserver(draw);
    resize.observe(canvas);
    let at = 0;
    return {
      enter(step, animate) {
        Deck.show(frame, step, animate);
        const forward = animate === true && (step === at + 1 || (step === 1 && at === 0));
        at = step;
        return document.fonts.load('600 16px Lato').catch(() => {}).then(() => go(step, forward));
      },
      freeze: draw,
      destroy() { dead = true; resize.disconnect(); },
    };
  },
});
}
