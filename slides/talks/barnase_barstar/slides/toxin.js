// Why the pair exists, as a manim-web cartoon: barnase cuts RNA; barstar plugs the site; the pair is among the
// tightest known.
//
// Steps: 1 barnase, and an RNA strand that runs through its active site and is cut · 2 barstar docks over the
// site and the next strand passes untouched · 3 the affinity, and the question the talk answers.
{
const A = [-4.5, -0.3], RA = 1.45, RB = 1.2;       // barnase's centre and the two radii, world units
const SITE = A[0] + RA;                            // the active site: barnase's right-hand edge
const wave = (x0, y0, y1, n = 60) => Array.from({ length: n + 1 }, (_, k) => {
  const y = y0 + (y1 - y0) * k / n;
  return [x0 + 0.16 * Math.sin(3.1 * y), y, 0];
});

Deck.slide('toxin', {
  steps: 3,
  notes: 'Barnase is a small ribonuclease that Bacillus amyloliquefaciens secretes. Inside the cell that makes '
    + 'it, it would destroy the cell\'s own RNA, so the same cell makes barstar, which binds barnase and blocks '
    + 'its active site until it is exported. Kd is about 1e-14 M (Schreiber and Fersht; Hartley gives 1e-13): '
    + 'among the tightest protein-protein complexes known. The rest of the talk asks which contacts do that.',
  mount(frame, { esc }) {
    frame.innerHTML = `
      <div class="deck-stage"></div>
      <h2 class="deck-heading reveal-in in">${esc('A toxin and its antidote')}</h2>
      ${Deck.points([
        { title: 'Barnase cuts RNA', sub: 'a ribonuclease a soil bacterium secretes; inside its own cell it would be lethal', from: 1 },
        { title: 'Barstar switches it off', sub: 'made by the same cell, it sits on the active site until barnase is exported', from: 2 },
        { title: '<b class="tx-kd">K<sub>d</sub> ≈ 10<sup>−14</sup> M</b>', sub: 'among the tightest protein pairs known. <span class="deck-accent">Which contacts do that?</span>', html: true, from: 3 },
      ], { className: 'tx-col' })}`;
    Deck.css(`.tx-kd { font-weight: inherit; display: inline-block; } .tx-kd sub { vertical-align: sub; font-size: .65em; } .tx-kd sup { vertical-align: super; font-size: .65em; }
      .tx-col { --col-left: 60%; --col-top: 24%; --col-gap: calc(2.6 * var(--u)); }`);

    const stage = Deck.manim(frame.querySelector('.deck-stage'), {
      build: ({ M, P }) => {
        const col = Tango.chainColours(P);
        const strand = (pts, colour) => {
          const v = new M.VMobject();
          v.setPointsAsCorners(pts);
          v.setColor(colour).setStrokeWidth(6).setFillOpacity(0);
          return v;
        };
        const text = (s, at, colour, size = 30) => new M.Text({ text: s, fontSize: size, fontFamily: 'Lato', fontWeight: 600, color: colour }).moveTo([at[0], at[1], 0]);
        return {
          barnase: new M.Circle({ radius: RA, color: col.a, fillOpacity: 0.2, strokeWidth: 6 }).moveTo([A[0], A[1], 0]),
          nameA: text('barnase', [A[0], A[1]], col.a),
          top: strand(wave(SITE + 0.12, A[1], 2.9), P.fg),
          bottom: strand(wave(SITE + 0.12, -3.5, A[1]), P.fg),
          rna: text('RNA', [SITE + 0.95, 2.6], P.muted, 26),
          // barstar starts off to the right and slides in
          barstar: new M.Circle({ radius: RB, color: col.b, fillOpacity: 0.2, strokeWidth: 6 }).moveTo([SITE + RB - 0.12 + 4, A[1], 0]),
          nameB: text('barstar', [SITE + RB - 0.12 + 4, A[1]], col.b),
          spared: strand(wave(SITE + 2 * RB + 0.6, -3.5, 2.9), P.fg),
          intact: text('intact', [SITE + 2 * RB + 1.35, 2.6], P.muted, 26),
        };
      },
      steps: [
        async ({ barnase, nameA, top, bottom, rna }, { M, play, wait, fast }) => {
          await play(new M.Create(barnase, { duration: 0.9 }), new M.FadeIn(nameA));
          await play(new M.Create(bottom, { duration: 0.7 }), new M.Create(top, { duration: 0.7 }), new M.FadeIn(rna));
          await wait(0.3);
          if (!fast) await play(new M.Flash(new M.Dot({ point: [SITE + 0.12, A[1], 0], radius: 0.01 }), { color: P_ACCENT(), flashRadius: 0.55, duration: 0.5 }));
          await play(new M.Shift(top, { direction: [0.55, 0.45, 0], duration: 0.7 }), new M.Shift(rna, { direction: [0.55, 0.45, 0], duration: 0.7 }),
            new M.Shift(bottom, { direction: [0.45, -0.4, 0], duration: 0.7 }));
        },
        async ({ top, bottom, rna, barstar, nameB, spared, intact }, { M, play }) => {
          await play(new M.FadeOut(top), new M.FadeOut(bottom), new M.FadeOut(rna), new M.FadeIn(barstar), new M.FadeIn(nameB));
          await play(new M.Shift(barstar, { direction: [-4, 0, 0], duration: 1.1 }), new M.Shift(nameB, { direction: [-4, 0, 0], duration: 1.1 }));
          await play(new M.Create(spared, { duration: 1.0 }), new M.FadeIn(intact));
        },
        ({ barnase, barstar }, { M }) => [new M.Indicate(barnase, { scaleFactor: 1.04 }), new M.Indicate(barstar, { scaleFactor: 1.04 })],
      ],
    });
    const P_ACCENT = () => Deck.palette.series[4];

    return {
      enter(step, animate) { Deck.show(frame, step, animate); return stage.go(step, animate); },
      freeze: stage.freeze,
      destroy: stage.dispose,
    };
  },
});
}
