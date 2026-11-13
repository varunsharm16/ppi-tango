// The real pair, in 3D, from the first frame of the simulation.
//
// Steps: 1 the two proteins apart · 2 they dock · 3 the camera moves into the interface, the cartoons thin out
// and the six residues the story follows appear as sticks, named.
{
const STATES = [
  { apart: 9, zoom: 1.2, at: [0, 0, 0], toonA: 1, toonB: 1, sticksA: 0, sticksB: 0, over: 0 },
  { apart: 0, zoom: 1.4, at: [0, 0, 0], toonA: 1, toonB: 1, sticksA: 0, sticksB: 0, over: 0 },
  { apart: 0, zoom: 2.7, at: [0, -1, 1], toonA: 0.22, toonB: 0.22, sticksA: 1, sticksB: 1, over: 1 },
];

Deck.slide('complex', {
  steps: 3,
  interactive: true,
  notes: 'Barnase in blue (110 residues, 108 of them in the structure), barstar in gold (89, 87 in the structure). The structure is PDB 1BRS, chains A and D, after '
    + 'equilibration. On step 3 the six residues are: barnase Lys27, Arg83, Arg87 and His102, and barstar Asp39 '
    + 'and Tyr29. Drag to turn it if someone asks; the next arrow key carries on.',
  mount(frame, { palette: P, data, esc }) {
    const D = data.interface;
    const n = (tag) => D.chains[tag].residues.length;
    frame.innerHTML = `
      <div class="deck-stage turnable"></div>
      <h2 class="deck-heading reveal-in in">${esc('Two small proteins, one interface')}</h2>
      <p class="cx-tag cx-a reveal-in" data-show="1">barnase<span>${n('a')} residues</span></p>
      <p class="cx-tag cx-b reveal-in" data-show="1">barstar<span>${n('b')} residues</span></p>
      <p class="deck-caption reveal-in" data-show="2">${esc('Docked: barstar covers the pocket where barnase would hold RNA.')}</p>
      <p class="deck-caption reveal-in" data-show="3">${esc('Six residues at the centre of the interface. We will follow these.')}</p>`;
    Deck.css(`.cx-tag { position: absolute; top: 77%; width: 30%; margin: 0; text-align: center; font-size: var(--t-lead);
        font-weight: 700; }
      .cx-tag span { display: block; font-size: var(--t-small); font-weight: 400; color: var(--muted); }
      .cx-a { left: 10%; color: var(--accent); } .cx-b { left: 60%; color: var(--plot-2); }`);
    const st = Tango.stage(frame.querySelector('.deck-stage'), D, P);
    st.show({ labels: D.sticks.residues.map((r) => r.id) });
    let at = 0;
    return {
      enter(step, animate) {
        st.enter();
        Deck.show(frame, step, animate);
        const from = at;
        at = step;
        if (animate === true && step === from + 1) return st.to(STATES[step - 1], step === 3 ? 1800 : 1400);
        st.set(STATES[step - 1]);
      },
      leave: st.leave,
      freeze: st.render,
      destroy: st.dispose,
    };
  },
});
}
