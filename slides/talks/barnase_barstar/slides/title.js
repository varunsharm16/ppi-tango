// The title slide. One step; the keys line is for whoever opens the deck, not part of the talk.
Deck.slide('title', {
  steps: 1,
  notes: 'A ten-minute story about one protein pair, told from a 20 ns simulation and one free energy '
    + 'calculation. Everything is driven by the arrow keys. The 3D slides can be dragged and zoomed if someone '
    + 'asks to see something, and the next arrow key carries on.',
  mount(frame, { css, esc }) {
    const KEYS = [['→', 'next'], ['←', 'back'], ['F', 'fullscreen'], ['S', 'speaker notes'], ['O', 'overview'], ['T', 'dark / light'], ['L', 'laser']];
    frame.innerHTML = `
      <h1 class="deck-title reveal-in" data-from="1">${esc('The tightest handshake')}</h1>
      <p class="deck-subtitle reveal-in" data-from="1" style="--i:1">${esc('How barstar disarms barnase, one contact at a time')}</p>
      <div class="deck-byline reveal-in" data-from="1" style="--i:2">${esc('Barnase–barstar, PDB 1BRS')}<span>${esc('20 ns of molecular dynamics and one free energy calculation · made with ppi-tango')}</span></div>
      <p class="ti-keys deck-detail reveal-in" data-from="1" style="--i:3">${KEYS.map(([k, what]) =>
        `<span><kbd>${esc(k)}</kbd>${esc(what)}</span>`).join('')}</p>`;
    css(`.ti-keys { position: absolute; left: 8%; right: 4.2%; top: 86.5%; margin: 0; display: flex; flex-wrap: wrap;
        gap: .4em 1.3em; font-size: var(--t-caption); }
      .ti-keys kbd { margin-right: .45em; }`);
    return { enter: (step, animate) => Deck.show(frame, step, animate) };
  },
});
