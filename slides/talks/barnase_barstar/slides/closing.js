// The take-home, and how the talk was made. Out of every section, so no footer.
//
// Steps: 1-3 the three points, one per click · 4 the two commands and the credits.
{
const KEEP = 0.95, RT = 0.5925;

Deck.slide('closing', {
  steps: 4,
  notes: 'Three things to remember. Then, if asked how this was made: ppi-tango runs PandaPPI on every frame of '
    + 'the trajectory; one script turns that into the data file these slides draw. The MD and FEP runs are in '
    + 'the repository under examples/simulations.',
  mount(frame, { data, esc }) {
    const D = data.interface, F = D.fep;
    const kept = D.contacts.filter((c) => c.occupancy >= KEEP);
    const ddg = F.dG[F.dG.length - 1];
    const fold = Math.exp(ddg / RT), mag = 10 ** Math.floor(Math.log10(fold));
    const MADE = [
      'ppitango md top.pdb traj.dcd --chains A D --dt 0.1',
      'python slides/tools/prep_interface.py top.pdb traj.dcd --chains A D --focus D:39 D:29 --dg dg.csv data/interface.json',
    ];
    frame.innerHTML = `
      <h2 class="deck-heading reveal-in in">${esc('What holds the handshake')}</h2>
      ${Deck.points([
        { title: 'A few contacts carry the interface', sub: `${kept.length} of ${D.contacts.length} are always there; ${kept.filter((c) => c.b === 'Asp39').length} of those go through barstar Asp39`, from: 1 },
        { title: 'Barstar wins by imitation', sub: 'Asp39 takes the place of the RNA phosphate in barnase\'s active site', from: 2 },
        { title: 'Not all of it is charge', sub: `the Tyr29 side chain alone is worth ${ddg.toFixed(1)} kcal/mol, about ${(Math.round(fold / mag) * mag).toLocaleString('en')}× in affinity`, from: 3 },
      ], { className: 'cl-col' })}
      <div class="cl-made reveal-in" data-from="4">${MADE.map((m) => `<code class="deck-code">${esc(m)}</code>`).join('')}</div>
      <p class="cl-credits deck-detail reveal-in" data-from="4">Contacts: <a href="https://github.com/pritampanda15/PandaPPI">PandaPPI</a> ·
        slides: <a href="https://github.com/ramithuh/protein-design-interactive-slides">live slides template</a>,
        <a href="https://github.com/maloyan/manim-web">manim-web</a>, <a href="https://threejs.org">three.js</a>,
        <a href="https://revealjs.com">Reveal.js</a> · simulations: <a href="https://openmm.org">OpenMM</a> ·
        structure: PDB 1BRS</p>`;
    Deck.css(`.cl-col { --col-left: 8%; --col-top: 20%; --col-gap: calc(2.2 * var(--u)); --col-title: var(--t-lead); }
      .cl-made { position: absolute; left: 8%; right: 4.2%; top: 68%; display: flex; flex-direction: column; gap: .6em;
        font-size: var(--t-micro); }
      .cl-made code { white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
      .deck-frame .cl-credits { position: absolute; left: 8%; right: 4.2%; top: 84%; margin: 0; font-size: var(--t-micro); }
      .deck-frame .cl-credits a { color: inherit; text-decoration: underline; text-decoration-color: color-mix(in srgb, currentColor 40%, transparent); }`);
    return { enter: (step, animate) => Deck.show(frame, step, animate) };
  },
});
}
