// The other residue the story follows: barstar Tyr29, uncharged, and what a free energy calculation does to it.
//
// Steps: 1 the camera moves to Tyr29 and its hydrogen bond, the simulation still running · 2 the simulation
// stops and the side chain beyond Cβ is switched off, as the FEP calculation does it.
{
Deck.slide('tyrosine', {
  steps: 2,
  interactive: true,
  notes: 'Tyr29 carries no charge. Its hydroxyl hydrogen-bonds to barnase Arg83 (the closest N-O pair is shown; '
    + 'present in about 93% of frames) and its ring packs against the barnase surface. Step 2 shows what the '
    + 'free energy calculation does: the side-chain atoms beyond C-beta stop interacting with everything else, '
    + 'first their charges, then their size, in 12 steps. The atoms never leave the coordinates; here they fade '
    + 'to show how much of them the rest of the system still feels.',
  mount(frame, { palette: P, data, esc }) {
    const D = data.interface;
    const pair = D.pairs.findIndex((p) => p.b === 'Tyr29');
    const hbond = D.contacts.find((c) => c.b === 'Tyr29' && /hydrogen/i.test(c.type));
    frame.innerHTML = `
      <div class="deck-stage turnable"></div>
      <h2 class="deck-heading reveal-in in">${esc('The quiet one: Tyr29')}</h2>
      <div class="ty-meter reveal-in" data-from="2"><span>${esc('Tyr29 side chain, felt by its surroundings')}</span><div><i></i></div><b></b></div>
      <p class="deck-caption reveal-in" data-show="1">${esc(`No charge: one hydrogen bond (there ${Math.round(100 * hbond.occupancy)}% of the time) and a ring packed against barnase.`)}</p>
      <p class="deck-caption reveal-in" data-show="2">${esc('What is it worth? Switch the side chain off, in the computer, and count the cost.')}</p>`;
    Deck.css(`.ty-meter { position: absolute; right: 4.2%; top: 17%; width: 27%; font-size: var(--t-small); color: var(--muted); }
      .ty-meter div { height: calc(.8 * var(--u)); margin-top: .5em; border-radius: 1em; background: color-mix(in srgb, var(--muted) 22%, transparent); overflow: hidden; }
      .ty-meter i { display: block; height: 100%; width: 100%; background: var(--plot-2); border-radius: 1em; }
      .ty-meter b { display: block; margin-top: .3em; font-size: var(--t-lead); color: var(--fg); font-variant-numeric: tabular-nums; }`);
    const bar = frame.querySelector('.ty-meter i'), pct = frame.querySelector('.ty-meter b');
    const st = Tango.stage(frame.querySelector('.deck-stage'), D, P, { rock: 0.3, residues: ['Tyr29', 'Arg83'] });
    st.onPaint = (S) => { bar.style.width = `${100 * S.ghost}%`; pct.textContent = `${Math.round(100 * S.ghost)}%`; };
    const t = st.centre('Tyr29'), r = st.centre('Arg83');
    const here = [0, 1, 2].map((c) => 0.65 * t[c] + 0.35 * r[c]);    // between the two, nearer the tyrosine
    const STATES = [
      { zoom: 4.8, at: here, toonA: 0.16, toonB: 0.16, sticksA: 1, sticksB: 1, over: 1, ghost: 1 },
      { zoom: 4.8, at: here, toonA: 0.16, toonB: 0.16, sticksA: 1, sticksB: 1, over: 1, ghost: 0 },
    ];
    let at = 0;
    return {
      enter(step, animate) {
        st.enter();
        st.play(step === 1);
        st.show({ labels: ['Tyr29', 'Arg83'], pairs: step === 1 ? [pair] : [] });   // nothing left to measure to
        Deck.show(frame, step, animate);
        const from = at;
        at = step;
        if (animate === true && step === 2 && from === 1) return st.to(STATES[1], 3600);
        st.set(STATES[step - 1]);
      },
      leave() { st.play(false); st.leave(); },
      freeze: st.render,
      destroy: st.dispose,
    };
  },
});
}
