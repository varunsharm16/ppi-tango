"""Pack a talk into one self-contained HTML file that can be sent to someone (chat attachment, email, a web host).

    python tools/bundle.py talks/barnase_barstar                    # writes talks/barnase_barstar/barnase_barstar.html
    python tools/bundle.py talks/barnase_barstar -o ~/Desktop/talk.html

Every <script src> and <link rel="stylesheet"> of the talk's index.html is inlined, and fonts the stylesheets
point to become data URIs, so the file needs nothing next to it and no network. The result is written without
<html>/<head>/<body> tags (browsers add them), which also lets hosts that wrap a page in their own skeleton
take it as is. Check it like any deck: python tools/check_deck.py is for folders, so open the file and step
through, or see the README.
"""

import argparse
import base64
import mimetypes
import re
from pathlib import Path


def inline_css(path: Path) -> str:
    css = path.read_text(encoding="utf-8")

    def data_uri(m):
        ref = m.group(2)
        if ref.startswith(("data:", "http:", "https:", "#")):
            return m.group(0)
        target = (path.parent / ref.split("?")[0].split("#")[0]).resolve()
        mime = mimetypes.guess_type(target.name)[0] or "application/octet-stream"
        return f'url("data:{mime};base64,{base64.b64encode(target.read_bytes()).decode()}")'

    return re.sub(r"url\((['\"]?)([^)'\"]+)\1\)", data_uri, css)


def bundle(talk: Path) -> str:
    html = (talk / "index.html").read_text(encoding="utf-8")

    def link(m):
        return f"<style>\n{inline_css(talk / m.group(1))}\n</style>"

    def script(m):
        js = (talk / m.group(1)).read_text(encoding="utf-8")
        # a literal "</script" inside the code would end the tag early
        return "<script>\n" + re.sub(r"</(script)", r"<\\/\1", js, flags=re.I) + "\n</script>"

    tags = r'<link rel="stylesheet" href="([^"]+)">', r'<script src="([^"]+)"></script>'
    left = re.findall(r'(?:src|href)="(?!https?:|data:|#)[^"]+"', re.sub("|".join(tags), "", html))
    if left:  # a local file this script does not know how to inline
        raise SystemExit(f"not inlined: {left}")
    html = re.sub(r"</?(?:html|head|body)\b[^>]*>\n?", "", html)  # before inlining: code may hold such tags in strings
    html = re.sub(tags[0], link, html)
    return re.sub(tags[1], script, html)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("talk", type=Path, help="a talk folder holding index.html")
    ap.add_argument("-o", "--out", type=Path, help="output file (default: <talk>/<talk name>.html)")
    args = ap.parse_args()
    out = args.out or args.talk / f"{args.talk.resolve().name}.html"
    out.write_text(bundle(args.talk), encoding="utf-8")
    print(f"{out}  {out.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
