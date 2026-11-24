"""The animated 2D interface report: PandaPPI's interaction map, alive.

Contact lines fade in and out with their (smoothed) occupancy, distances tick, and the side panels
follow the contact counts and buried surface along the mode / trajectory / lambda axis.
"""
import subprocess

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation
from matplotlib.lines import Line2D
from pandappi.visualization.interaction_map import (BURIED_MIN, CHAIN_COLORS, CLASS_COLOR, GREY, INK, NATURE_RC,
                                                    RED, RESIDUE_CLASS, TYPE_STYLE, _style)

from .analysis import Scene

RC = {**NATURE_RC, 'font.size': 11, 'axes.labelsize': 11, 'axes.titlesize': 12, 'xtick.labelsize': 9.5,
      'ytick.labelsize': 10, 'legend.fontsize': 9.5, 'axes.linewidth': 0.8, 'xtick.major.width': 0.8,
      'ytick.major.width': 0.8, 'xtick.major.size': 3, 'ytick.major.size': 3}
AMINO = {r for r in RESIDUE_CLASS if len(r) == 3}


def _short(res: str) -> str:
    """'ASP 100A' -> 'Asp 100A': amino acids in mixed case, everything else (nucleotides, ligands) as it is."""
    name, _, rest = res.partition(' ')
    return f'{name.capitalize() if name in AMINO else name} {rest}'


DIRECTIONAL = ('Hydrogen bond', 'Salt bridge / ionic', 'Disulfide')  # the contacts that get a distance label


def sample(arr: np.ndarray, u: float, periodic: bool) -> np.ndarray:
    """Value of a per-frame array (frames on the last axis) at fractional frame u."""
    n = arr.shape[-1]
    i = int(np.floor(u))
    f = u - i
    j = (i + 1) % n if periodic else min(i + 1, n - 1)
    i = i % n if periodic else min(i, n - 1)
    return arr[..., i] * (1 - f) + arr[..., j] * f


def position_text(scene: Scene, u: float) -> str:
    ens = scene.ensemble
    x = sample(ens.x, u, False) if not ens.periodic else (u / ens.n_frames) % 1
    if ens.kind == 'md':
        return f"t = {x:.1f} {ens.xlabel.split('(')[-1].rstrip(')')}"
    if ens.kind == 'fep':
        return f'λ = {x:.2f}'
    return f'phase {x:.2f}'


class Report:
    """Builds the figure once, then update(u) moves every artist to fractional frame u."""

    def __init__(self, scene: Scene, title: str = '', images=None, labels=None):
        """With `images` (one rendered 3D picture per movie frame) the layout becomes a split screen:
        the models on the left, the map and two small timelines on the right. `labels` then holds, per
        movie frame, the (x, y, text, colour) distance labels to write on the picture, in pixels."""
        self.s, ens = scene, scene.ensemble
        self.periodic = ens.periodic
        a, b = ens.chains
        with plt.rc_context(RC):
            self.fig = fig = plt.figure(figsize=(12.8, 7.2), dpi=150)
            fig.text(0.025, 0.955, title or ens.title, fontsize=17, fontweight='bold', va='center')
            self.status = fig.text(0.025, 0.915, '', fontsize=11.5, va='center')
            fig.text(0.975, 0.012, 'ppi-tango · contacts by PandaPPI' + ''.join(f' · {n}' for n in ens.notes),
                     fontsize=8, ha='right', color='#8C8C8C')
            self.type_bars, self.type_series, self.area_bars, self.area_series = [], [], [], []
            if images is None:
                self._map(fig.add_axes([0.02, 0.13, 0.47, 0.74]), a, b, legend=(0.03, 0.012), note=(0.255, 0.118))
                self._timeline(fig.add_axes([0.585, 0.70, 0.38, 0.16]), fig.add_axes([0.585, 0.46, 0.38, 0.16]))
                self._bars(fig.add_axes([0.655, 0.115, 0.10, 0.20]), fig.add_axes([0.855, 0.115, 0.11, 0.20]), a, b)
            else:
                ax = fig.add_axes([0.005, 0.035, 0.49, 0.85])
                ax.axis('off')
                self.picture = ax.imshow(plt.imread(images[0]), interpolation='lanczos')
                self.images, self.labels = images, labels or [[]] * len(images)
                self.tags = [ax.text(0, 0, '', fontsize=9.5, fontweight='bold', ha='center', va='center', visible=False,
                                     bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.8))
                             for _ in range(max(map(len, self.labels), default=0))]
                self._map(fig.add_axes([0.50, 0.375, 0.49, 0.515]), a, b, legend=(0.515, 0.25), note=(0.745, 0.362),
                          ncol=3)
                self._timeline(fig.add_axes([0.545, 0.095, 0.185, 0.10]), fig.add_axes([0.79, 0.095, 0.175, 0.10]),
                               compact=True)
        self.update(0)

    # ---- panel a: the interaction map
    def _map(self, ax, a, b, legend, note, ncol=4):
        s = self.s
        n = max(len(s.left), len(s.right), 2)
        spread = lambda res: {r: n - 1 - i * (n - 1) / max(len(res) - 1, 1) if len(res) > 1 else (n - 1) / 2
                              for i, r in enumerate(res)}
        ly, ry = spread(s.left), spread(s.right)
        t = np.linspace(0, 1, 60)
        rank = list(TYPE_STYLE)
        pairs = {}
        for k in s.drawn:
            c = s.contacts[k]
            pairs.setdefault((c['residue_a'], c['residue_b']), []).append(k)
        self.edges = []
        for i, ((l, r), ks) in enumerate(pairs.items()):
            ks.sort(key=lambda k: rank.index(s.contacts[k]['type']) if s.contacts[k]['type'] in rank else 99)
            for m, k in enumerate(ks):
                bend = (m - (len(ks) - 1) / 2) * 0.5  # fan out several contact types between one pair
                y0, y1 = ly[l], ry[r]
                ys = (1 - t) ** 2 * y0 + 2 * (1 - t) * t * ((y0 + y1) / 2 + bend) + t ** 2 * y1
                kind = s.contacts[k]['type']
                color, ls, badge = _style(kind)
                line, = ax.plot(t, ys, color=color, ls=ls, lw=1.6, zorder=1, solid_capstyle='round')
                j = int(len(t) * (0.2 + 0.55 * ((i * 0.618 + m * 0.31) % 1)))  # stagger badges along the line
                btxt = ax.text(t[j], ys[j], badge, fontsize=8, color=color, ha='center', va='center', zorder=3,
                               bbox=dict(boxstyle='circle,pad=0.2', fc='white', ec=color, lw=0.8))
                dtxt = None
                if kind in DIRECTIONAL:
                    j2 = min(j + 4, len(t) - 1)
                    dtxt = ax.text(t[j2], ys[j2], '', fontsize=8, color=color, zorder=2, ha='left', va='bottom',
                                   bbox=dict(boxstyle='round,pad=0.08', fc='white', ec='none'))
                self.edges.append((k, line, btxt, dtxt))
        self.nodes = []
        for chain, pos, x, ha, side in ((a, ly, 0.0, 'right', 'residue_a'), (b, ry, 1.0, 'left', 'residue_b')):
            for res, y in pos.items():
                txt = ax.text(x, y, _short(res), ha=ha, va='center', fontsize=10, fontweight='bold', zorder=4,
                              bbox=dict(boxstyle='round,pad=0.3', ec=INK, lw=0.6,
                                        fc=CLASS_COLOR[RESIDUE_CLASS.get(res.split()[0], 'Other')]))
                mine = [k for k in s.drawn if s.contacts[k][side] == res]
                self.nodes.append((txt, s.occ[mine].max(0), s.area_of(chain, res)))
        for name, x, ha, color in ((s.ensemble.sides[0], 0.0, 'right', CHAIN_COLORS[0]),
                                   (s.ensemble.sides[1], 1.0, 'left', CHAIN_COLORS[1])):
            ax.text(x, n + 0.1, name, ha=ha, va='bottom', fontsize=13, fontweight='bold', color=color)
        ax.set_xlim(-0.32, 1.32)
        ax.set_ylim(-0.8, n + 0.9)
        ax.axis('off')
        shown = [k for k in TYPE_STYLE if any(s.contacts[i]['type'] == k for i in s.drawn)]
        handles = [Line2D([], [], color=_style(k)[0], ls=_style(k)[1], lw=1.6, label=f'{_style(k)[2]}  {k}')
                   for k in shown]
        ax.figure.legend(handles=handles, loc='lower left', bbox_to_anchor=legend, ncol=ncol,
                         handlelength=2.2, columnspacing=1.6)
        text = f'heavy outline: buries ≥ {BURIED_MIN:g} Å²'
        if s.hidden:
            text += f' · {s.hidden} more contact{"s" * (s.hidden > 1)} on less busy residues not drawn'
        ax.figure.text(*note, text, fontsize=8.5, ha='center', color='#8C8C8C')

    # ---- timeline: contacts by type, buried surface (and ΔG)
    def _timeline(self, top, bottom, compact=False):
        s, ens = self.s, self.s.ensemble
        x = ens.x
        # a mode oscillation closes on itself: repeat the first frame at the end of the axis
        wrap = (lambda v: np.append(v, v[0])) if ens.periodic else (lambda v: v)
        x, bsa = (np.append(x, 1.0) if ens.periodic else x), wrap(s.bsa)
        top.plot(x, wrap(s.total), color=INK, lw=2, label='All contacts')
        for kind in s.types[:3]:
            color, ls, _ = _style(kind)
            top.plot(x, wrap(s.type_counts[kind]), color=color, ls=ls, lw=1.5, label=kind)
        top.set_ylim(0, max(s.total.max(), 1) * 1.15)
        if compact:
            top.set_title('Contacts', loc='left', fontweight='bold')
            top.set_xlabel(ens.xlabel)
        else:
            top.legend(loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=4, handlelength=1.6, columnspacing=1.0,
                       borderaxespad=0.2)
            top.set_xticklabels([])
        bottom.plot(x, bsa, color=INK, lw=2)
        bottom.fill_between(x, bsa, s.bsa.min() - 0.15 * (np.ptp(s.bsa) + 1), color=GREY, alpha=0.25, lw=0)
        bottom.set_ylim(s.bsa.min() - 0.15 * (np.ptp(s.bsa) + 1), s.bsa.max() + 0.15 * (np.ptp(s.bsa) + 1))
        bottom.set_title('Buried surface (Å²)', loc='left', fontweight='bold')
        bottom.set_xlabel(ens.xlabel)
        self.cursors = [ax.axvline(x[0], color=INK, lw=1, alpha=0.6) for ax in (top, bottom)]
        self.dots = [top.plot([x[0]], [s.total[0]], 'o', color=INK, ms=6)[0],
                     bottom.plot([x[0]], [s.bsa[0]], 'o', color=INK, ms=6)[0]]
        self.curve_dot = None
        if ens.curve is not None:
            twin = bottom.twinx()
            twin.spines['right'].set_visible(True)
            twin.spines['right'].set_color(RED)
            twin.plot(ens.x, ens.curve, color=RED, lw=2)
            twin.set_ylabel(ens.curve_label, color=RED)
            twin.tick_params(axis='y', colors=RED)
            self.curve_dot, = twin.plot([x[0]], [ens.curve[0]], 'o', color=RED, ms=6)
        for ax in (top, bottom):
            ax.set_xlim(x[0], x[-1])

    # ---- live bars: contacts by type, most buried residues
    def _bars(self, left, right, a, b):
        s = self.s
        kinds = s.types[:7][::-1]
        self.type_bars = left.barh(range(len(kinds)), [0] * len(kinds), height=0.65,
                                   color=[_style(k)[0] for k in kinds])
        self.type_series = [s.type_counts[k] for k in kinds]
        left.set_yticks(range(len(kinds)), kinds)
        left.set_xlim(0, max(max(v.max() for v in self.type_series), 1) * 1.1)
        left.set_title('Contacts by type', loc='left', fontweight='bold', x=-1.45)
        left.set_xlabel('contacts (n)')
        top = np.argsort(-s.area.mean(1))[:7][::-1]
        self.area_bars = right.barh(range(len(top)), [0] * len(top), height=0.65,
                                    color=[CHAIN_COLORS[0] if s.residues[i][0] == a else CHAIN_COLORS[1] for i in top])
        self.area_series = s.area[top]
        # the chain in front of each residue, unless its name already carries one (a side of several chains)
        ens = s.ensemble
        front = {c: side[6:] + ':' if side.startswith('Chain ') else '' for c, side in zip(ens.chains, ens.sides)}
        right.set_yticks(range(len(top)), [f'{front[s.residues[i][0]]}{_short(s.residues[i][1])}' for i in top])
        right.set_xlim(0, max(self.area_series.max(), 1) * 1.1)
        right.set_title('Most buried residues', loc='left', fontweight='bold', x=-0.95)
        right.set_xlabel('buried (Å²)')

    def picture_frame(self, i: int):
        """Split screen only: show rendered picture i with its distance labels."""
        self.picture.set_data(plt.imread(self.images[i]))
        for tag in self.tags:
            tag.set_visible(False)
        for tag, (x, y, text, color) in zip(self.tags, self.labels[i]):
            tag.set(position=(x, y), text=text, color=color, visible=True)

    def update(self, u: float):
        s, ens, per = self.s, self.s.ensemble, self.periodic
        for k, line, badge, dist in self.edges:
            o = float(sample(s.occ[k], u, per))
            line.set_alpha(o)
            line.set_linewidth(0.6 + 1.4 * o)
            badge.set_alpha(o)
            badge.get_bbox_patch().set_alpha(o)
            if dist is not None:
                dist.set_text(f'{sample(s.dist[k], u, per):.1f} Å')
                dist.set_alpha(o)
                dist.get_bbox_patch().set_alpha(0.85 * o)
        for txt, active, area in self.nodes:
            o = 0.3 + 0.7 * float(sample(active, u, per))
            txt.set_alpha(o)
            patch = txt.get_bbox_patch()
            patch.set_alpha(o)
            patch.set_linewidth(1.8 if sample(area, u, per) >= BURIED_MIN else 0.6)
        x = (u / ens.n_frames) % 1 if per else float(sample(ens.x, u, False))
        total, bsa = float(sample(s.total, u, per)), float(sample(s.bsa, u, per))
        for c in self.cursors:
            c.set_xdata([x, x])
        self.dots[0].set_data([x], [total])
        self.dots[1].set_data([x], [bsa])
        if self.curve_dot is not None:
            self.curve_dot.set_data([x], [float(sample(ens.curve, u, False))])
        for bar, series in zip(self.type_bars, self.type_series):
            bar.set_width(float(sample(series, u, per)))
        for bar, series in zip(self.area_bars, self.area_series):
            bar.set_width(float(sample(series, u, per)))
        hb = float(sample(s.type_counts.get('Hydrogen bond', np.zeros(ens.n_frames)), u, per))
        text = f'{position_text(s, u)}   ·   {total:.0f} contacts ({hb:.0f} hydrogen bonds)   ·   ' \
               f'buried surface {bsa:,.0f} Å²'
        if ens.curve is not None:
            text += f"   ·   {ens.curve_label.split(' (')[0]} {sample(ens.curve, u, False):+.1f} kcal/mol"
        self.status.set_text(text)


def positions(scene: Scene, n: int) -> np.ndarray:
    """Fractional frame shown by each of the n movie frames."""
    ens = scene.ensemble
    span = ens.n_frames if ens.periodic else ens.n_frames - 1
    return np.arange(n) / (n if ens.periodic else n - 1) * span


def render(scene: Scene, output: str, seconds: float = 12.0, fps: int = 30, title: str = '', dpi: int = 150,
           images=None, labels=None) -> str:
    """Write the animated report to `output` (.mp4). `images`, `labels`: see Report."""
    import imageio_ffmpeg
    report = Report(scene, title, images, labels)
    n = len(images) if images else max(int(seconds * fps), 2)
    us = positions(scene, n)

    def draw(i):
        if images:
            report.picture_frame(i)
        report.update(us[i])
    anim = FuncAnimation(report.fig, draw, frames=n)
    plt.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()
    # yuv420p + even dimensions: plays in PowerPoint, Keynote, QuickTime and browsers
    writer = FFMpegWriter(fps=fps, codec='libx264', extra_args=['-pix_fmt', 'yuv420p', '-crf', '18'])
    try:
        with plt.rc_context(RC):
            anim.save(output, writer=writer, dpi=dpi)
    except (subprocess.CalledProcessError, OSError):   # ffmpeg gave up; its own words are printed by matplotlib
        raise _encoder_stopped(output) from None
    finally:
        plt.close(report.fig)
    return output


def _encoder_stopped(output: str) -> OSError:
    return OSError(f'could not write {output}: the video encoder (ffmpeg) stopped. Is the disk full, or the file '
                   'not writable?')


def to_gif(mp4: str, width: int = 960, fps: int = 12) -> str:
    """Convert a movie to a GIF next to it (for READMEs and chat; slides should use the MP4)."""
    import imageio_ffmpeg
    gif = mp4[:-4] + '.gif'
    chain = (f'fps={fps},scale={width}:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=96[p];'
             f'[b][p]paletteuse=dither=none')
    try:
        subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(), '-y', '-loglevel', 'error', '-i', mp4, '-vf', chain, gif],
                       check=True)
    except (subprocess.CalledProcessError, OSError):
        raise _encoder_stopped(gif) from None
    return gif


def snapshot(scene: Scene, output: str, frame: float = 0, title: str = '') -> str:
    """A single still of the report (PNG/PDF), at fractional frame `frame`."""
    report = Report(scene, title)
    report.update(frame)
    with plt.rc_context(RC):
        report.fig.savefig(output, dpi=150)
    plt.close(report.fig)
    return output
