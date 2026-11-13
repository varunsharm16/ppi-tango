// What the persistent contacts have in common, as a manim-web diagram: one aspartate, four partners, and their
// distances over the whole simulation.
//
// Steps: 1 Asp39 and the four barnase residues round it, with mean distances · 2 the four distances traced
// over 20 ns · 3 what those four residues are: barnase's catalytic site.
{
const C = [-4.3, -0.5];                 // Asp39 in the diagram
const AX = { x0: -2.5, x1: 1.15, y0: -2.9, y1: 0.9, lo: 2, hi: 6 };   // the plot's box (world) and its Å range

Deck.slide('hub', {
  steps: 3,
  notes: 'Distances are the closest N-O pair between Asp39 and each partner, per frame. Lys27, Arg83, Arg87 and '
    + 'His102 are barnase\'s active site: His102 is the general acid, Lys27 stabilises the transition state, and '
    + 'the two arginines hold the phosphate. Barstar puts a carboxylate exactly where the RNA phosphate would '
    + 'go (Buckle, Schreiber and Fersht, Biochemistry 1994). It does not just cover the active site, it '
    + 'imitates the substrate.',
  mount(frame, { palette: P, data, esc }) {
    const D = data.interface;
    const hub = D.pairs.filter((p) => p.b === 'Asp39');
    const REACH = 3.5;   // Å: the usual donor-acceptor limit for a hydrogen bond
    const together = Math.round(100 * hub[0].d.filter((_, f) => hub.every((p) => p.d[f] <= REACH)).length / hub[0].d.length);
    frame.innerHTML = `
      <div class="deck-stage"></div>
      <h2 class="deck-heading reveal-in in">${esc('One aspartate, four partners')}</h2>
      <p class="hb-axis reveal-in" data-from="2">${esc(`closest N–O distance, Å, over ${D.ns.toFixed(0)} ns`)}</p>
      ${Deck.points([
        { title: 'Four residues converge on Asp39', sub: 'two arginines, a lysine and a histidine', from: 1 },
        { title: 'They stay within hydrogen-bond reach', sub: `all four at once in ${together}% of the ${D.frames} frames`, from: 2 },
        { title: 'These four are the active site', sub: 'Asp39 sits where the RNA phosphate would. <span class="deck-accent">Barstar imitates the substrate.</span>', html: true, from: 3 },
      ], { className: 'hb-col' })}`;
    Deck.css(`.hb-col { --col-left: 62%; --col-top: 17%; --col-gap: calc(1.7 * var(--u)); }
      .hb-axis { position: absolute; left: 34.4%; top: 84%; margin: 0; font-size: var(--t-micro); color: var(--muted); }
      .hb-col .deck-point:nth-child(n) { max-width: 100%; }`);
    // the plot shares the right half with the points: it sits low, under points 1 and 2
    const px = (t) => AX.x0 + (AX.x1 - AX.x0) * t;
    const py = (d) => AX.y0 + (AX.y1 - AX.y0) * (Math.min(AX.hi, d) - AX.lo) / (AX.hi - AX.lo);

    const stage = Deck.manim(frame.querySelector('.deck-stage'), {
      build: ({ M }) => {
        const col = Tango.chainColours(P);
        const text = (s, at, colour, size = 24, weight = 600) => new M.Text({ text: s, fontSize: size, fontFamily: 'Lato', fontWeight: weight, color: colour }).moveTo([at[0], at[1], 0]);
        const tint = [P.series[2], P.series[3], P.series[1], P.series[4]];
        const spokes = hub.map((p, k) => {
          const a = Math.PI * (0.62 + 0.76 * k / (hub.length - 1)), r = 2.25;     // a fan on Asp39's left
          const at = [C[0] + r * Math.cos(a), C[1] + r * Math.sin(a) * 1.15];
          const dir = [Math.cos(a), Math.sin(a) * 1.15], len = Math.hypot(dir[0], dir[1]);
          const unit = [dir[0] / len, dir[1] / len];
          const line = new M.DashedLine({ start: [C[0] + unit[0] * 0.62, C[1] + unit[1] * 0.62, 0],
            end: [at[0] - unit[0] * 0.5, at[1] - unit[1] * 0.5, 0], color: tint[k], strokeWidth: 5, dashLength: 0.13 });
          const mid = [(C[0] + at[0]) / 2 + unit[1] * 0.42, (C[1] + at[1]) / 2 - unit[0] * 0.42];
          const trace = new M.VMobject();
          trace.setPointsAsCorners(p.d.map((d, f) => [px(f / (p.d.length - 1)), py(d), 0]));
          trace.setColor(tint[k]).setStrokeWidth(2.5).setFillOpacity(0);
          return {
            node: new M.Circle({ radius: 0.48, color: col.a, fillOpacity: 0.2, strokeWidth: 5 }).moveTo([at[0], at[1], 0]),
            name: text(p.a, at, col.a, 21),
            line, dist: text(`${p.mean.toFixed(1)} Å`, mid, tint[k], 20), trace,
          };
        });
        const axis = (a, b) => new M.Line({ start: [a[0], a[1], 0], end: [b[0], b[1], 0], color: P.muted, strokeWidth: 2 });
        const ticks = [2, 3, 4, 5, 6].map((d) => text(String(d), [AX.x0 - 0.25, py(d)], P.muted, 17, 400));
        const reach = new M.DashedLine({ start: [AX.x0, py(REACH), 0], end: [AX.x1, py(REACH), 0], color: P.detail, strokeWidth: 2, dashLength: 0.1 });
        return {
          asp: new M.Circle({ radius: 0.6, color: P.el.O, fillOpacity: 0.22, strokeWidth: 6 }).moveTo([C[0], C[1], 0]),
          aspName: text('Asp39', C, P.el.O, 23),
          minus: text('−', [C[0] + 0.5, C[1] + 0.52], P.el.O, 34, 700),
          spokes,
          axes: [axis([AX.x0, AX.y0], [AX.x1, AX.y0]), axis([AX.x0, AX.y0], [AX.x0, AX.y1]), reach, ...ticks,
            text('3.5 Å', [AX.x1 + 0.4, py(REACH)], P.detail, 17, 400)],
        };
      },
      steps: [
        async ({ asp, aspName, minus, spokes }, { M, play }) => {
          await play(new M.GrowFromCenter(asp), new M.FadeIn(aspName), new M.FadeIn(minus));
          await play(new M.LaggedStart(spokes.map((s) => new M.AnimationGroup([new M.GrowFromCenter(s.node), new M.FadeIn(s.name)])), { lagRatio: 0.25 }));
          await play(...spokes.flatMap((s) => [new M.Create(s.line, { duration: 0.7 }), new M.FadeIn(s.dist)]));
        },
        async ({ spokes, axes }, { M, play }) => {
          await play(...axes.map((a) => new M.FadeIn(a)));
          await play(...spokes.map((s) => new M.Create(s.trace, { duration: 2.6 })));
        },
        ({ spokes }, { M }) => [new M.LaggedStart(spokes.map((s) => new M.Indicate(s.node, { scaleFactor: 1.12 })), { lagRatio: 0.2 })],
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
