// The simulation, live: the cartoons and side chains follow the 20 ns trajectory, back and forth.
//
// Steps: 1 the whole complex moving, with the clock · 2 into barstar Asp39 and the four barnase residues round
// it, each with its distance measured in every frame.
{
Deck.slide('dynamics', {
  steps: 2,
  interactive: true,
  notes: 'A 20 ns explicit-solvent simulation (OpenMM, Amber14SB, TIP3P, 150 mM NaCl, 298 K), one frame every '
    + '0.1 ns. The drawn coordinates are a 5-frame rolling mean so the ribbons do not shiver; the distances are '
    + 'measured on the raw frames. Step 2: the four dashed lines are the closest N-O pair between Asp39 and '
    + 'each partner. Watch the numbers: Arg83, Arg87 and Lys27 hardly move; His102 breathes a little.',
  mount(frame, { palette: P, data, esc }) {
    const D = data.interface;
    const hub = D.pairs.map((p, k) => (p.b === 'Asp39' ? k : -1)).filter((k) => k >= 0);
    const names = ['Asp39', ...hub.map((k) => D.pairs[k].a)];
    frame.innerHTML = `
      <div class="deck-stage turnable"></div>
      <h2 class="deck-heading reveal-in in">${esc(`${D.ns.toFixed(0)} nanoseconds of motion`)}</h2>
      <p class="dy-clock reveal-in in"><b></b> ns</p>
      <p class="deck-caption reveal-in" data-show="1">${esc('Everything jiggles. Which contacts survive it?')}</p>
      <p class="deck-caption reveal-in" data-show="2">${esc('Barstar Asp39 and its four partners on barnase: distances measured live, in Å.')}</p>`;
    Deck.css(`.dy-clock { position: absolute; right: 4.2%; top: 6.2%; margin: 0; font-size: var(--t-lead); color: var(--muted);
        font-variant-numeric: tabular-nums; }
      .dy-clock b { color: var(--fg); }`);
    const clock = frame.querySelector('.dy-clock b');
    const st = Tango.stage(frame.querySelector('.deck-stage'), D, P, { rock: 0.22, residues: names });
    const tick = (f) => { clock.textContent = D.time[Math.round(f)].toFixed(1); };
    st.onFrame = tick;
    tick(0);
    const STATES = [
      { zoom: 1.4, at: [0, 0, 0], toonA: 1, toonB: 1, sticksA: 0, sticksB: 0, over: 0, ghost: 0 },
      { zoom: 3.4, at: st.centre('Asp39'), toonA: 0.16, toonB: 0.16, sticksA: 1, sticksB: 1, over: 1, ghost: 0 },
    ];
    st.show({ labels: names, pairs: hub });
    let at = 0;
    return {
      enter(step, animate) {
        st.enter();
        st.play(true);
        Deck.show(frame, step, animate);
        const from = at;
        at = step;
        if (animate === true && step === from + 1 && step > 1) return st.to(STATES[step - 1], 1800);
        st.set(STATES[step - 1]);
      },
      leave() { st.play(false); st.leave(); },
      freeze: st.render,
      destroy: st.dispose,
    };
  },
});
}
