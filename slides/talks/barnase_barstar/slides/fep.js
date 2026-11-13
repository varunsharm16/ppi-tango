// The free energy calculation as a manim-web graph: the cost of switching the Tyr29 side chain off, in the
// complex and in barstar on its own, and the gap between them.
//
// Steps: 1 the axes and both curves, drawn window by window · 2 the gap at the end: ΔΔG · 3 what that means
// for binding, and the measured value.
{
const BOX = { x0: -6.5, x1: -0.9, y0: -2.6, y1: 2.5 };
const RT = 0.5925;                                      // kcal/mol at 298 K

Deck.slide('fep', {
  steps: 3,
  notes: 'Alchemical free energy perturbation with OpenMM: the Tyr29 side chain beyond C-beta is decoupled in 12 '
    + 'windows (charges first, then Lennard-Jones), 1 ns each, in the complex and in free barstar, analysed with '
    + 'MBAR. The difference between the two legs is the side chain\'s contribution to binding. Experiment for '
    + 'Tyr29 to Ala is 3.4 kcal/mol (Schreiber and Fersht). Caveats if asked: short sampling, statistical error '
    + 'only, and a real mutation to alanine would use a hybrid topology.',
  mount(frame, { palette: P, data, esc }) {
    const F = data.interface.fep, n = F.window.length;
    const ddg = F.dG[n - 1], err = F.dG_err[n - 1];
    const fold = Math.exp(ddg / RT);
    const round = (v) => { const m = 10 ** Math.floor(Math.log10(v)); return Math.round(v / m) * m; };
    const split = F.lambda_electrostatics.findIndex((l) => l === 0);     // the window where the charges are gone
    frame.innerHTML = `
      <div class="deck-stage"></div>
      <h2 class="deck-heading reveal-in in">${esc('What is one side chain worth?')}</h2>
      ${Deck.points([
        { title: 'Switch it off twice', sub: `in the complex and in barstar alone, in ${n} small steps`, from: 1 },
        { title: `${ddg.toFixed(1)} ± ${err.toFixed(1)} kcal/mol harder in the complex`, sub: 'the gap is what the side chain gives to binding', from: 2 },
        { title: `About ${round(fold).toLocaleString('en')}× weaker without it`, sub: 'measured for Tyr29→Ala: 3.4 kcal/mol <span class="deck-detail">(Schreiber &amp; Fersht)</span>', html: true, from: 3 },
      ], { className: 'fe-col' })}
      <p class="deck-caption deck-detail reveal-in" data-from="1">${esc('Free energy to switch off the Tyr29 side chain, kcal/mol · 1 ns per step · statistical error only')}</p>`;
    Deck.css(`.fe-col { --col-left: 60%; --col-top: 22%; --col-gap: calc(2.4 * var(--u)); }`);

    const top = Math.ceil(Math.max(...F.dG_complex) / 10) * 10;
    const px = (w) => BOX.x0 + (BOX.x1 - BOX.x0) * w / (n - 1);
    const py = (g) => BOX.y0 + (BOX.y1 - BOX.y0) * g / top;

    const stage = Deck.manim(frame.querySelector('.deck-stage'), {
      build: ({ M }) => {
        const text = (s, at, colour, size = 20, weight = 400) => new M.Text({ text: s, fontSize: size, fontFamily: 'Lato', fontWeight: weight, color: colour }).moveTo([at[0], at[1], 0]);
        const line = (a, b, colour, w = 2) => new M.Line({ start: [a[0], a[1], 0], end: [b[0], b[1], 0], color: colour, strokeWidth: w });
        const curve = (ys, colour) => {
          const v = new M.VMobject();
          v.setPointsAsCorners(ys.map((g, w) => [px(w), py(g), 0]));
          v.setColor(colour).setStrokeWidth(5).setFillOpacity(0);
          return v;
        };
        const cBound = P.series[4], cFree = P.series[0];
        const dots = (ys, colour) => ys.map((g, w) => new M.Dot({ point: [px(w), py(g), 0], radius: 0.07, color: colour }));
        const ticks = [];
        for (let g = 0; g <= top; g += 10) ticks.push(text(String(g), [BOX.x0 - 0.3, py(g)], P.muted, 18));
        const xr = px(n - 1) + 0.3;
        return {
          axes: [line([BOX.x0, BOX.y0], [BOX.x1 + 0.1, BOX.y0], P.muted), line([BOX.x0, BOX.y0], [BOX.x0, BOX.y1], P.muted), ...ticks,
            new M.DashedLine({ start: [px(split), BOX.y0, 0], end: [px(split), BOX.y1, 0], color: P.detail, strokeWidth: 2, dashLength: 0.1 }),
            text('charges off', [(px(0) + px(split)) / 2, BOX.y0 - 0.33], P.muted, 19),
            text('atoms off', [(px(split) + px(n - 1)) / 2, BOX.y0 - 0.33], P.muted, 19)],
          bound: curve(F.dG_complex, cBound), free: curve(F.dG_free, cFree),
          marks: [...dots(F.dG_complex, cBound), ...dots(F.dG_free, cFree)],
          names: [text('in the complex', [px(n - 4), py(F.dG_complex[n - 4]) + 0.42], cBound, 21, 600),
            text('barstar alone', [px(n - 3), py(F.dG_free[n - 3]) - 0.5], cFree, 21, 600)],
          gap: line([xr, py(F.dG_free[n - 1])], [xr, py(F.dG_complex[n - 1])], P.fg, 6),
          caps: [line([xr - 0.12, py(F.dG_free[n - 1])], [xr + 0.12, py(F.dG_free[n - 1])], P.fg, 4),
            line([xr - 0.12, py(F.dG_complex[n - 1])], [xr + 0.12, py(F.dG_complex[n - 1])], P.fg, 4)],
          ddg: text(`ΔΔG = ${ddg.toFixed(1)}`, [xr + 0.95, (py(F.dG_free[n - 1]) + py(F.dG_complex[n - 1])) / 2], P.fg, 24, 700),
        };
      },
      steps: [
        async ({ axes, bound, free, marks, names }, { M, play }) => {
          await play(...axes.map((a) => new M.FadeIn(a)));
          await play(new M.Create(bound, { duration: 2.2 }), new M.Create(free, { duration: 2.2 }));
          await play(...marks.map((m) => new M.FadeIn(m)), ...names.map((t) => new M.FadeIn(t)));
        },
        async ({ gap, caps, ddg }, { M, play }) => {
          await play(new M.Create(gap, { duration: 0.6 }), ...caps.map((c) => new M.FadeIn(c)));
          await play(new M.Write(ddg, { duration: 0.9 }));
        },
        ({ ddg }, { M }) => [new M.Indicate(ddg, { scaleFactor: 1.15 })],
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
