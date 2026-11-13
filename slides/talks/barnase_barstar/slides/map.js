// From the 3D interface to PandaPPI's flat map, as a manim-web morph, then the map coloured by charge.
//
// Steps: 1 every residue that touches the other chain, as a dot where its Cα is in the view of the last slide
// · 2 the dots fall into two columns and get their names · 3 one line per residue pair in contact, drawn in
// order of how persistent it is · 4 the residues recoloured by charge; only the salt bridges stay.
{
const CX = -3.35, CY = -0.25, PER_A = 0.17;          // where the 3D view sits, and world units per Å
const XA = -4.75, XB = -1.95, TOP = 2.75, BOTTOM = -3.45;
const kind = (type) => (/salt/i.test(type) ? 'ionic' : /hydrogen/i.test(type) ? 'hbond' : 'other');
const RANK = { ionic: 2, hbond: 1, other: 0 };

Deck.slide('map', {
  steps: 4,
  notes: 'This is the PandaPPI interaction map, built from the simulation instead of one structure. A line is a '
    + 'pair of residues that are in contact in at least 20% of the 200 frames; its thickness is how often. '
    + 'Step 4: barnase brings six positive residues to the interface (Lys27, Lys39, Arg59, Arg83, Arg87, His102), '
    + 'barstar four negative ones and no positive one (Asp35, Asp39, Glu46, Glu76). His102 is counted as positive '
    + 'because it sits against a carboxylate.',
  mount(frame, { palette: P, data, esc }) {
    const D = data.interface;
    // one line per residue pair: its strongest kind and its highest occupancy
    const pairs = new Map();
    for (const c of D.contacts) {
      const key = `${c.a}|${c.b}`, k = kind(c.type), had = pairs.get(key);
      if (!had) pairs.set(key, { a: c.a, b: c.b, kind: k, occ: c.occupancy });
      else { if (RANK[k] > RANK[had.kind]) had.kind = k; had.occ = Math.max(had.occ, c.occupancy); }
    }
    const lines = [...pairs.values()].sort((p, q) => q.occ - p.occ);
    const count = (tag, sign) => D.nodes[tag].filter((n) => Math.sign(n.charge) === sign).length;
    const KIND = { ionic: P.series[4], hbond: P.series[3], other: P.series[0] };
    const chip = (colour, label) => `<span><i style="background:${colour}"></i>${esc(label)}</span>`;
    frame.innerHTML = `
      <div class="deck-stage"></div>
      <h2 class="deck-heading reveal-in in">${esc('Flatten the interface')}</h2>
      ${Deck.points([
        { title: `${D.nodes.a.length + D.nodes.b.length} residues touch the other side`, sub: `${D.nodes.a.length} on barnase, ${D.nodes.b.length} on barstar, over ${D.ns.toFixed(0)} ns`, from: 1 },
        { title: 'Unfold them into two columns', sub: 'top to bottom as they sit in the interface', from: 2 },
        { title: `${D.contacts.length} contacts, ${lines.length} residue pairs`, sub: `<span class="mp-key">${chip(KIND.ionic, 'salt bridge')}${chip(KIND.hbond, 'hydrogen bond')}${chip(KIND.other, 'other')}</span>thicker = present more of the time`, html: true, from: 3 },
        { title: 'Opposite charges face each other', sub: `barnase: <b style="color:${P.el.N}">${count('a', 1)} positive</b> · barstar: <b style="color:${P.el.O}">${count('b', -1)} negative</b>, ${count('b', 1)} positive`, html: true, from: 4 },
      ], { className: 'mp-col' })}`;
    Deck.css(`.mp-col { --col-left: 58%; --col-top: 20%; --col-gap: calc(2.1 * var(--u)); }
      .mp-key { display: flex; gap: .9em; margin-bottom: .25em; white-space: nowrap; font-size: .86em; }
      .mp-key i { display: inline-block; width: 1.5em; height: .32em; border-radius: .2em; margin-right: .4em; vertical-align: .22em; }`);

    const stage = Deck.manim(frame.querySelector('.deck-stage'), {
      build: ({ M }) => {
        const col = Tango.chainColours(P);
        const charge = (n) => (n.charge > 0 ? P.el.N : n.charge < 0 ? P.el.O : P.detail);
        const nodes = {};
        for (const tag of ['a', 'b']) {
          const list = [...D.nodes[tag]].sort((p, q) => q.at[2] - p.at[2]);    // top of the interface first
          list.forEach((n, k) => {
            const y = TOP + (BOTTOM - TOP) * k / (list.length - 1), x = tag === 'a' ? XA : XB;
            nodes[n.id] = {
              n, tag, to: [x, y, 0],
              dot: new M.Dot({ point: [CX + n.at[0] * PER_A, CY + n.at[2] * PER_A, 0], radius: 0.1, color: col[tag] }),
              column: new M.Dot({ point: [x, y, 0], radius: 0.1, color: col[tag] }),
              charged: new M.Dot({ point: [x, y, 0], radius: n.charge ? 0.15 : 0.07, color: charge(n) }),
              name: new M.Text({ text: n.id, fontSize: 19, fontFamily: 'Lato', fontWeight: 600, color: P.fg })
                .moveTo([x + (tag === 'a' ? -0.72 : 0.72), y, 0]),
            };
          });
        }
        const link = lines.map((l) => {
          const a = nodes[l.a].to, b = nodes[l.b].to;
          const line = new M.Line({ start: [a[0] + 0.16, a[1], 0], end: [b[0] - 0.16, b[1], 0], color: KIND[l.kind],
            strokeWidth: 1.2 + 4.5 * l.occ });
          line.setStrokeOpacity(0.35 + 0.6 * l.occ);
          return { line, kind: l.kind };
        });
        const head = (s, x, colour) => new M.Text({ text: s, fontSize: 24, fontFamily: 'Lato', fontWeight: 700, color: colour }).moveTo([x, TOP + 0.48, 0]);
        return { nodes: Object.values(nodes), link, heads: [head('barnase', XA - 0.4, col.a), head('barstar', XB + 0.4, col.b)] };
      },
      steps: [
        ({ nodes }, { M }) => [new M.LaggedStart(nodes.map((n) => new M.GrowFromCenter(n.dot)), { lagRatio: 0.04 })],
        async ({ nodes, heads }, { M, play }) => {
          await play(...nodes.map((n) => new M.Transform(n.dot, n.column, { duration: 1.4 })));
          await play(...nodes.map((n) => new M.FadeIn(n.name)), ...heads.map((h) => new M.FadeIn(h)));
        },
        ({ link }, { M }) => [new M.LaggedStart(link.map((l) => new M.Create(l.line, { duration: 0.6 })), { lagRatio: 0.06 })],
        ({ nodes, link }, { M }) => [
          ...nodes.map((n) => new M.Transform(n.dot, n.charged, { duration: 0.9 })),
          ...link.filter((l) => l.kind !== 'ionic').map((l) => new M.FadeOut(l.line)),
        ],
      ],
    });

    return {
      enter(step, animate) { Deck.show(frame, step, animate); return stage.go(step, animate); },
      freeze: stage.freeze,
      destroy: stage.dispose,
    };
  },
});
}
