// The complex in 3D, shared by every 3D slide of this talk: both chains as cartoons that follow the MD
// trajectory, the residues the story follows as sticks, and an SVG layer over the stage for residue names and
// measured distances (dashed line + Å), all placed from the same atoms every redraw.
//
//   const st = Tango.stage(host, data.interface, palette, { residues: ['Asp39', 'Arg83'] });   // sticks: all, or these
//   st.set({ apart: 9, zoom: 1 });                     // a state, at once
//   await st.to({ apart: 0 }, 1400);                   // ...or eased over 1.4 s
//   st.show({ labels: ['Asp39', 'Arg83'], pairs: [0, 1] });   // what the overlay names and measures
//   st.play(true);                                     // run the trajectory back and forth
//   return { enter, leave: st.leave, freeze: st.render, destroy: st.dispose };
//
// State keys: apart (Å each chain is pulled back along x), toonA / toonB (cartoon opacity), sticksA / sticksB
// (stick opacity), ghost (opacity of the alchemical side chain, here Tyr29 beyond Cβ), over (overlay opacity),
// zoom and at ([x, y, z] Å, what the camera turns about). Coordinates are data/interface.js's frame.
{
const SCALE = 0.1;                       // scene units per Å
const GHOST = { id: 'Tyr29', keep: ['N', 'CA', 'C', 'O', 'CB'] };   // the side chain FEP switches off

window.Tango = {
  chainColours: (P) => ({ a: P.accent, b: P.series[1] }),

  stage(host, D, P, { rock = 0.3, fps = 14, residues = null } = {}) {
    const T = window.ManimWeb.THREE;
    const view = Deck.turntable(host, { radius: 36 * SCALE, scale: SCALE, rock, tilt: 0.1 });
    const col = Tango.chainColours(P);
    const K = D.sticks, frames = K.xyz.length;

    // ---- cartoons ----
    const chains = {};
    for (const tag of ['a', 'b']) {
      const mat = new T.MeshStandardMaterial({ color: new T.Color(col[tag]), transparent: true });
      const toon = Deck.cartoon(D.chains[tag], { material: mat });
      const group = new T.Group();
      group.add(toon.mesh);
      view.model.add(group);
      chains[tag] = { mat, toon, group };
    }

    // ---- sticks: one set per chain, and the ghost's side chain on its own so it can fade alone ----
    const resChain = K.residues.map((r) => r.chain);
    const ghostRes = K.residues.findIndex((r) => r.id === GHOST.id);
    const inGhost = (i) => K.residue[i] === ghostRes && !GHOST.keep.includes(K.names[i]);
    const drawn = (i) => !residues || residues.includes(K.residues[K.residue[i]].id);
    const member = {
      a: K.names.map((_, i) => drawn(i) && resChain[K.residue[i]] === 'a'),
      b: K.names.map((_, i) => drawn(i) && resChain[K.residue[i]] === 'b' && !inGhost(i)),
      // the ghost set also holds Cβ, so the Cβ-Cγ bond has both ends
      ghost: K.names.map((n, i) => drawn(i) && K.residue[i] === ghostRes && (inGhost(i) || n === 'CB')),
    };
    const sets = {};
    for (const [name, mask] of Object.entries(member)) {
      const idx = mask.map((m, i) => (m ? i : -1)).filter((i) => i >= 0);
      const local = new Map(idx.map((i, k) => [i, k]));
      const bonds = K.bonds.filter(([i, j]) => local.has(i) && local.has(j)).map(([i, j]) => [local.get(i), local.get(j)]);
      // carbons in a tint of the chain's colour, so nitrogens and oxygens stand out from it
      const carbon = new T.Color(col[name === 'ghost' ? 'b' : name]).lerp(new T.Color(P.bg), 0.45);
      const colours = idx.map((i) => (K.elements[i] === 'C' ? carbon.clone() : new T.Color(P.cpk[K.elements[i]] || P.cpk.other)));
      const set = Deck.ballStick(idx.map((i) => K.xyz[0][i]), colours, bonds, { atom: 0.34, bond: 0.17 });
      chains[name === 'a' ? 'a' : 'b'].group.add(set.group);
      sets[name] = { set, idx };
    }

    // ---- the overlay ----
    const NS = 'http://www.w3.org/2000/svg';
    const svg = document.createElementNS(NS, 'svg');
    svg.setAttribute('class', 'tg-over');
    host.after(svg);
    Deck.css(`.deck-frame .tg-over { position: absolute; inset: 0; width: 100%; height: 100%; pointer-events: none;
        font-family: Lato, sans-serif; font-weight: 600; }
      .tg-over text { paint-order: stroke; stroke: var(--bg); stroke-linejoin: round; }`);
    const make = (tag, attrs) => {
      const el = document.createElementNS(NS, tag);
      for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, v);
      svg.append(el);
      return el;
    };
    const pairEls = D.pairs.map(() => ({
      line: make('line', { stroke: P.fg, 'stroke-linecap': 'round' }),
      text: make('text', { fill: P.fg, 'text-anchor': 'middle', 'dominant-baseline': 'central' }),
    }));
    const leadEls = K.residues.map((r) => make('line', { stroke: col[r.chain], 'stroke-opacity': 0.55, 'stroke-linecap': 'round' }));
    const labelEls = K.residues.map((r) => make('text', { fill: col[r.chain], 'text-anchor': 'middle',
      'dominant-baseline': 'central' }));
    labelEls.forEach((el, k) => { el.textContent = K.residues[k].id; });

    // ---- state ----
    const S = { apart: 0, toonA: 1, toonB: 1, sticksA: 0, sticksB: 0, ghost: 1, over: 0, zoom: 1, at: [0, 0, 0] };
    let shown = { labels: [], pairs: [] };
    let frame = 0, dir = 1, playing = false, active = false, dead = false, raf = 0, gen = 0;
    let now = K.xyz[0];                   // every stick atom at the current frame
    const lerp3 = (p, q, u) => [p[0] + (q[0] - p[0]) * u, p[1] + (q[1] - p[1]) * u, p[2] + (q[2] - p[2]) * u];

    function setFrame(f) {
      frame = Math.max(0, Math.min(frames - 1, f));
      const i = Math.min(frames - 2, Math.floor(frame)), u = frame - i;
      for (const tag of ['a', 'b']) {
        const tr = D.chains[tag].trace;
        chains[tag].toon.update(tr[i].map((p, k) => lerp3(p, tr[i + 1][k], u)));
      }
      now = K.xyz[i].map((p, k) => lerp3(p, K.xyz[i + 1][k], u));
      for (const { set, idx } of Object.values(sets)) set.put(idx.map((k) => now[k]));
      view.dirty();
    }

    function paint() {
      chains.a.group.position.x = -S.apart;
      chains.b.group.position.x = S.apart;
      for (const [tag, a] of [['a', S.toonA], ['b', S.toonB]]) {
        const c = chains[tag];
        c.mat.opacity = a; c.mat.depthWrite = a > 0.99; c.toon.mesh.visible = a > 0.001;
      }
      sets.a.set.opacity(S.sticksA);
      sets.b.set.opacity(S.sticksB);
      sets.ghost.set.opacity(S.sticksB * S.ghost);
      svg.style.opacity = S.over;
      if (api.onPaint) api.onPaint(S);
      view.dirty();
    }

    // names sit beyond their residue, pushed away from the middle of everything named; distances sit on
    // their dashed line
    view.onRender = () => {
      const u = host.clientWidth / 100;
      const i = Math.min(frames - 1, Math.round(frame));
      const centroid = (k) => {
        const pts = now.filter((_, a) => K.residue[a] === k);
        return pts.reduce((s, p) => [s[0] + p[0] / pts.length, s[1] + p[1] / pts.length, s[2] + p[2] / pts.length], [0, 0, 0]);
      };
      const at = K.residues.map((r, k) => (shown.labels.includes(r.id) ? view.project(centroid(k)) : null));
      const on = at.map((p, k) => (p ? k : -1)).filter((k) => k >= 0);
      const mid = on.length ? [on.reduce((s, k) => s + at[k][0], 0) / on.length, on.reduce((s, k) => s + at[k][1], 0) / on.length] : [0, 0];
      // names go on a ring round the residues, each joined to its residue by a leader; neighbours on the ring
      // are pushed apart until they are at least GAP radians from each other
      const ring = Math.max(14 * u, ...on.map((k) => Math.hypot(at[k][0] - mid[0], at[k][1] - mid[1]) + 8 * u));
      const order = on.map((k) => ({ k, a: Math.atan2(at[k][1] - mid[1], at[k][0] - mid[0]) })).sort((p, q) => p.a - q.a);
      const GAP = 0.42;
      for (let pass = 0; pass < 12 && order.length > 1; pass++) {
        order.forEach((p, n) => {
          const q = order[(n + 1) % order.length];
          const gap = ((q.a - p.a) % (2 * Math.PI) + 2 * Math.PI) % (2 * Math.PI);
          if (gap < GAP) { p.a -= (GAP - gap) / 2; q.a += (GAP - gap) / 2; }
        });
      }
      const angle = new Map(order.map((p) => [p.k, p.a]));
      labelEls.forEach((el, k) => {
        const lead = leadEls[k];
        el.style.display = lead.style.display = at[k] ? '' : 'none';
        if (!at[k]) return;
        const a = angle.get(k), x = mid[0] + Math.cos(a) * ring * 1.25, y = mid[1] + Math.sin(a) * ring * 0.8;
        el.setAttribute('x', x);
        el.setAttribute('y', y);
        el.setAttribute('font-size', 1.9 * u);
        el.setAttribute('stroke-width', 0.5 * u);
        const d = [x - at[k][0], y - at[k][1]], l = Math.hypot(d[0], d[1]) || 1;
        lead.setAttribute('x1', at[k][0] + d[0] / l * 1.2 * u);
        lead.setAttribute('y1', at[k][1] + d[1] / l * 1.2 * u);
        lead.setAttribute('x2', x - d[0] / l * 3.4 * u);
        lead.setAttribute('y2', y - d[1] / l * 1.6 * u);
        lead.setAttribute('stroke-width', 0.14 * u);
      });
      pairEls.forEach(({ line, text }, k) => {
        const vis = shown.pairs.includes(k);
        line.style.display = text.style.display = vis ? '' : 'none';
        if (!vis) return;
        const p = D.pairs[k];
        const [x1, y1] = view.project(now[p.i[i]]), [x2, y2] = view.project(now[p.j[i]]);
        for (const [key, v] of [['x1', x1], ['y1', y1], ['x2', x2], ['y2', y2]]) line.setAttribute(key, v);
        line.setAttribute('stroke-width', 0.28 * u);
        line.setAttribute('stroke-dasharray', `${0.5 * u} ${0.6 * u}`);
        text.textContent = `${p.d[D.kept[i]].toFixed(1)} Å`;
        // beside its line, on alternate sides so neighbouring distances do not sit on each other
        const len = Math.hypot(x2 - x1, y2 - y1) || 1, side = (k % 2 ? -1 : 1) * 2.1 * u;
        text.setAttribute('x', (x1 + x2) / 2 - (y2 - y1) / len * side);
        text.setAttribute('y', (y1 + y2) / 2 + (x2 - x1) / len * side);
        text.setAttribute('font-size', 1.55 * u);
        text.setAttribute('stroke-width', 0.55 * u);
      });
    };

    let prev = performance.now();
    const tick = (t) => {
      const dt = Math.min(0.1, (t - prev) / 1000);
      prev = t;
      if (active && playing) {
        let f = frame + dir * fps * dt;
        if (f >= frames - 1) { f = frames - 1; dir = -1; } else if (f <= 0) { f = 0; dir = 1; }
        setFrame(f);
        if (api.onFrame) api.onFrame(frame);
      }
      if (!dead) raf = requestAnimationFrame(tick);
    };
    raf = requestAnimationFrame(tick);

    const api = {
      view, data: D,
      get frame() { return frame; },
      onFrame: null,     // called with the frame while the trajectory plays
      onPaint: null,     // called with the state whenever it changes (to drive a readout)
      setFrame,
      show(what) { shown = { labels: [], pairs: [], ...what }; view.dirty(); },
      play(on) { playing = on; },
      // a state at once: any of the keys above
      set(state) {
        gen++;
        Object.assign(S, state);
        view.focus(S.at, S.zoom, 0);
        paint();
      },
      // ...or eased there over ms; a newer set()/to() takes over
      to(state, ms = 1200) {
        const id = ++gen, from = { ...S }, t0 = performance.now();
        const keys = Object.keys(state).filter((k) => k !== 'at' && k !== 'zoom');
        if ('at' in state) S.at = state.at;
        if ('zoom' in state) S.zoom = state.zoom;
        view.focus(S.at, S.zoom, ms);
        return new Promise((done) => {
          const step = (t) => {
            if (dead || id !== gen) return done();
            const x = Math.min(1, (t - t0) / ms), e = x < 0.5 ? 4 * x * x * x : 1 - (-2 * x + 2) ** 3 / 2;
            for (const k of keys) S[k] = from[k] + (state[k] - from[k]) * e;
            paint();
            if (x < 1) requestAnimationFrame(step); else done();
          };
          requestAnimationFrame(step);
        });
      },
      // the middle of a residue's atoms in frame 0, Å: somewhere for the camera to turn about
      centre(id) {
        const k = K.residues.findIndex((r) => r.id === id);
        const pts = K.xyz[0].filter((_, a) => K.residue[a] === k);
        return [0, 1, 2].map((c) => pts.reduce((s, p) => s + p[c], 0) / pts.length);
      },
      enter() { active = true; prev = performance.now(); view.start(); },
      leave() { active = false; view.stop(); },
      render() { view.render(); },
      dispose() {
        dead = true;
        cancelAnimationFrame(raf);
        for (const c of Object.values(chains)) { c.toon.dispose(); c.mat.dispose(); }
        for (const { set } of Object.values(sets)) set.dispose();
        svg.remove();
        view.dispose();
      },
    };
    setFrame(0);
    paint();
    return api;
  },
};
}
