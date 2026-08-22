# ruff: noqa: E501

"""
Make annotated X-ray image review page for a set of obsids.

Reads processed outputs from the astromon pipeline work directory and
the astromon HDF5 database, then writes a self-contained HTML file with:
  - log1p-scaled X-ray images (broad-band thresholded counts)
  - circle overlays for each celldetect source, colour-coded by flags
  - acis_streak polygon overlays
  - hover tooltips with per-source details

Usage
-----
    python -m astromon.scripts.analysis.review_xray_images \\
        --obsids 30270 29490 30282 \\
        --workdir /Volumes/Black/data/astromon/work \\
        --db     /Volumes/Black/data/astromon/astromon.h5 \\
        --out    xray_review.html

Image scaling
-------------
Block-sum (4x4 pixels), then clip at the 99.5th percentile of nonzero
blocks, then log1p.  This puts single-photon blocks at ~8-20% of the
colour-map range (dark purple in inferno) while avoiding the very bright
source from flattening the scale for the rest of the field.

Coordinate system
-----------------
Chandra celldetect source X,Y are physical sky pixel coordinates.
The image FITS header provides LTV1/LTV2 so that
    img_col = X_phys + LTV1
    img_row = Y_phys + LTV2
The image array is stored top-down (row 0 = lowest Y in sky), so source
positions use img_row/H directly (no Y flip) and the canvas renders rows
0..H-1 top to bottom.  The same transform is applied to streak polygon
vertices.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import tables as tb
from astropy.io import fits

# ---------------------------------------------------------------------------
# Data extraction
# ---------------------------------------------------------------------------


def load_db_sources(obsid: int, db_path: Path) -> dict[int, dict]:
    """Return {id: {snr, acis_streak, grating_arm}} for celldetect sources."""
    result: dict[int, dict] = {}
    with tb.open_file(str(db_path)) as h5:
        for row in h5.root.astromon_xray_src.where(f"obsid == {obsid}"):
            method = row["detect_method"]
            if isinstance(method, bytes):
                method = method.decode()
            if method != "celldetect":
                continue
            src_id = int(row["id"])
            result[src_id] = {
                "snr": float(row["snr"]),
                "acis_streak": bool(row["acis_streak"]),
                "grating_arm": bool(row["grating_arm"]),
            }
    return result


def load_cat_src_rfc(obsid: int, db_path: Path) -> set[int]:
    """Return set of celldetect x_ids that have an RFC entry in cat_src."""
    rfc_ids: set[int] = set()
    with tb.open_file(str(db_path)) as h5:
        for row in h5.root.astromon_cat_src.where(f"obsid == {obsid}"):
            cat = row["catalog"]
            if isinstance(cat, bytes):
                cat = cat.decode()
            if cat.strip() in ("RFC", "ICRS"):
                rfc_ids.add(int(row["x_id"]))
    return rfc_ids


def read_image(obs_dir: Path) -> tuple[np.ndarray, float, float, int, int]:
    """Return (data, ltv1, ltv2, W, H) from the broad_thresh.img."""
    img_path = next(obs_dir.glob("images/*_broad_thresh.img"))
    with fits.open(img_path) as hdul:
        hdr = hdul[0].header
        data = hdul[0].data.astype(np.float32)
    return (
        data,
        float(hdr["LTV1"]),
        float(hdr["LTV2"]),
        int(hdr["NAXIS1"]),
        int(hdr["NAXIS2"]),
    )


def downsample(data: np.ndarray, factor: int = 4, clip_pct: float = 99.5) -> list:
    """Block-sum, clip at nonzero percentile, log1p-normalise to 0-255 bytes.

    Sum pooling preserves single-event pixels.  Clipping at the 99.5th
    percentile of nonzero blocks prevents a single bright source from
    dominating the scale; log1p then makes individual photons visible.
    """
    H, W = data.shape
    H2, W2 = H // factor, W // factor
    ds = (
        data[: H2 * factor, : W2 * factor]
        .reshape(H2, factor, W2, factor)
        .sum(axis=(1, 3))
    )
    nz = ds[ds > 0]
    if len(nz) == 0:
        return ds.astype(np.uint8).tolist()
    clip_val = float(np.percentile(nz, clip_pct))
    ds_clipped = np.clip(ds, 0, clip_val)
    ds_log = np.log1p(ds_clipped)
    vmax = ds_log.max()
    if vmax > 0:
        ds_log = ds_log / vmax * 255
    return ds_log.astype(np.uint8).tolist()


def phys_to_norm(
    x_phys: float, y_phys: float, ltv1: float, ltv2: float, W: int, H: int
) -> tuple[float, float]:
    """Physical Chandra sky pixel → normalised canvas [0, 1].

    No Y-flip: image rows are stored with row 0 at the top of the
    display, matching the canvas rendering direction.
    """
    img_col = x_phys + ltv1
    img_row = y_phys + ltv2
    return float(img_col / W), float(img_row / H)


def read_sources(
    obs_dir: Path,
    db_srcs: dict[int, dict],
    rfc_ids: set[int],
    ltv1: float,
    ltv2: float,
    W: int,
    H: int,
) -> list[dict]:
    """Source positions + flags in normalised canvas coordinates."""
    src_file = next(obs_dir.glob("sources/*_psf_size_celldetect.fits"))
    with fits.open(src_file) as hdul:
        rows = hdul[1].data
    results = []
    for i, row in enumerate(rows):
        src_id = i + 1
        db_info = db_srcs.get(
            src_id, {"snr": 1.0, "acis_streak": False, "grating_arm": False}
        )
        xn, yn = phys_to_norm(float(row["X"]), float(row["Y"]), ltv1, ltv2, W, H)
        results.append(
            {
                "id": src_id,
                "x": round(xn, 5),
                "y": round(yn, 5),
                "snr": round(db_info["snr"], 1),
                "streak": db_info["acis_streak"],
                "arm": db_info["grating_arm"],
                "rfc": src_id in rfc_ids,
            }
        )
    return results


def read_streak_polys(
    obs_dir: Path, ltv1: float, ltv2: float, W: int, H: int
) -> list[list[dict]]:
    """Streak polygon vertex lists in normalised canvas coordinates."""
    streak_file = next(obs_dir.glob("images/*_acis_streaks.fits"), None)
    if streak_file is None:
        return []
    polys = []
    with fits.open(streak_file) as hdul:
        for ext in hdul[1:]:
            if ext.data is None:
                continue
            for row in ext.data:
                xs = np.asarray(row["X"], dtype=float)
                ys = np.asarray(row["Y"], dtype=float)
                verts = []
                for xv, yv in zip(xs, ys, strict=True):
                    xn, yn = phys_to_norm(xv, yv, ltv1, ltv2, W, H)
                    verts.append({"x": round(xn, 5), "y": round(yn, 5)})
                if verts:
                    polys.append(verts)
    return polys


# ---------------------------------------------------------------------------
# HTML generation (self-contained)
# ---------------------------------------------------------------------------

HTML_TEMPLATE = r"""<title>Chandra Source Detection Review</title>
<style>
*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
:root {
  --bg: #090b0f; --surface: #111520; --border: #1e2d40; --text: #c4d4e8; --muted: #4e6a85;
  --c-streak: #e05050; --c-arm: #c89020; --c-rfc: #38c87a; --c-normal: #4fa8d8;
}
html, body { background: var(--bg); color: var(--text);
  font-family: ui-monospace,'Cascadia Code','SF Mono',Menlo,monospace;
  font-size: 13px; line-height: 1.5; min-height: 100vh; }
.page { max-width: 1160px; margin: 0 auto; padding: 28px 18px 60px; }
header { margin-bottom: 22px; border-bottom: 1px solid var(--border); padding-bottom: 16px; }
h1 { font-size: 17px; font-weight: 600; }
header p { margin-top: 7px; color: var(--muted); font-size: 12px; max-width: 660px;
  font-family: ui-sans-serif,system-ui,sans-serif; line-height: 1.65; }
.legend { display: flex; flex-wrap: wrap; gap: 5px 18px; margin-bottom: 22px;
  padding: 10px 14px; background: var(--surface); border: 1px solid var(--border); border-radius: 3px; }
.li { display: flex; align-items: center; gap: 7px; font-size: 11px; }
.dot { width:10px; height:10px; border-radius:50%; flex-shrink:0; }
.ring { width:10px; height:10px; border-radius:50%; border:2px solid #38c87a; flex-shrink:0; }
.sw { width:16px; height:8px; border-radius:2px; flex-shrink:0; }
.panels { display: grid; grid-template-columns: repeat(auto-fill, minmax(320px, 1fr)); gap: 18px; }
.panel { background: var(--surface); border: 1px solid var(--border); border-radius: 4px; overflow: hidden; }
.ph { padding: 9px 13px 8px; border-bottom: 1px solid var(--border);
  display: flex; align-items: baseline; gap: 9px; }
.oid { font-size: 14px; font-weight: 700; }
.lbl { font-size: 11px; color: var(--muted); font-family: ui-sans-serif,system-ui,sans-serif; }
.cw { position:relative; line-height:0; cursor:crosshair; }
canvas { display:block; width:100%; image-rendering:pixelated; }
.ps { padding: 9px 13px; border-top: 1px solid var(--border); font-size: 11px; color: var(--muted);
  font-family: ui-sans-serif,system-ui,sans-serif; display:flex; flex-wrap:wrap; gap: 3px 14px; }
.sv { color:var(--text); font-family:ui-monospace,monospace; }
.note { margin-top:18px; padding:11px 15px; background:var(--surface);
  border:1px solid var(--border); border-left:3px solid #3a7fd4;
  font-family:ui-sans-serif,system-ui,sans-serif; font-size:12px; color:var(--muted); line-height:1.65; }
.note strong { color:var(--text); }
#tt { position:fixed; background:#0c1422ee; border:1px solid #1e2d40; border-radius:3px;
  padding:7px 11px; font-size:11px; color:#c4d4e8; line-height:1.7;
  pointer-events:none; z-index:100; white-space:nowrap; display:none; }
.fs{color:#e05050} .fa{color:#c89020} .fr{color:#38c87a} .fo{color:#4e6a85}
</style>
<div class="page">
  <header>
    <h1>Chandra Source Detection — acis_streak review</h1>
    <p>TITLE_PLACEHOLDER</p>
  </header>
  <div class="legend">
    <div class="li"><span class="dot" style="background:#4fa8d8"></span>normal</div>
    <div class="li"><span class="dot" style="background:#e05050"></span>acis_streak=True</div>
    <div class="li"><span class="dot" style="background:#c89020"></span>grating_arm=True</div>
    <div class="li"><span class="ring"></span>RFC in cat_src</div>
    <div class="li"><span class="sw" style="background:rgba(220,60,60,.22);border:1px solid rgba(220,60,60,.65)"></span>streak polygon</div>
  </div>
  <div class="panels" id="panels"></div>
  <div class="note" id="note">Loading…</div>
</div>
<div id="tt"></div>
<script>
const VIZ = DATA_PLACEHOLDER;
const CMAP=[[[0,4,4,12],[.12,40,14,68],[.30,120,24,100],[.55,200,80,24],[.78,244,168,18],[1,252,250,210]]];
function cm(t){t=Math.max(0,Math.min(1,t));const s=CMAP[0];let i=0;while(i<s.length-2&&t>s[i+1][0])i++;const f=(t-s[i][0])/(s[i+1][0]-s[i][0]);return s[i].slice(1).map((v,j)=>Math.round(v+f*(s[i+1][j+1]-v)));}
function sc(s){if(s.streak)return'#e05050';if(s.arm)return'#c89020';return'#4fa8d8';}
function build(d){
  const{obsid,label,image,img_shape,sources,streak_polys}=d;
  const[dsH,dsW]=img_shape,S=360,sx=S/dsW,sy=S/dsH;
  const panel=document.createElement('div');panel.className='panel';
  const ph=document.createElement('div');ph.className='ph';
  ph.innerHTML=`<span class="oid">${obsid}</span><span class="lbl">${label}</span>`;
  panel.appendChild(ph);
  const cw=document.createElement('div');cw.className='cw';
  const cv=document.createElement('canvas');cv.width=S;cv.height=S;
  cw.appendChild(cv);panel.appendChild(cw);
  const ctx=cv.getContext('2d');
  const id=ctx.createImageData(S,S);
  for(let r=0;r<dsH;r++)for(let c=0;c<dsW;c++){
    const[R,G,B]=cm(image[r][c]/255);
    const x0=Math.round(c*sx),y0=Math.round(r*sy),x1=Math.round((c+1)*sx),y1=Math.round((r+1)*sy);
    for(let py=y0;py<y1;py++)for(let px=x0;px<x1;px++){const i=(py*S+px)*4;id.data[i]=R;id.data[i+1]=G;id.data[i+2]=B;id.data[i+3]=255;}
  }
  ctx.putImageData(id,0,0);
  streak_polys.forEach(poly=>{
    if(!poly.length)return;
    ctx.save();ctx.beginPath();ctx.moveTo(poly[0].x*S,poly[0].y*S);
    for(let i=1;i<poly.length;i++)ctx.lineTo(poly[i].x*S,poly[i].y*S);
    ctx.closePath();ctx.fillStyle='rgba(220,60,60,0.18)';ctx.fill();
    ctx.strokeStyle='rgba(220,60,60,0.7)';ctx.lineWidth=1;ctx.setLineDash([3,3]);ctx.stroke();ctx.restore();
  });
  [...sources].sort((a,b)=>a.snr-b.snr).forEach(src=>{
    const cx=src.x*S,cy=src.y*S,r=Math.max(6,Math.min(22,4+Math.sqrt(src.snr)*1.1)),col=sc(src);
    ctx.save();
    ctx.beginPath();ctx.arc(cx,cy,r,0,2*Math.PI);ctx.fillStyle=col+'22';ctx.fill();
    ctx.beginPath();ctx.arc(cx,cy,r,0,2*Math.PI);ctx.strokeStyle=col;ctx.lineWidth=src.snr>30?2.5:1.5;ctx.stroke();
    if(src.rfc){ctx.beginPath();ctx.arc(cx,cy,r+5,0,2*Math.PI);ctx.strokeStyle='#38c87a';ctx.lineWidth=1.5;ctx.setLineDash([3,2]);ctx.stroke();ctx.setLineDash([]);}
    if(src.snr>=4.5||src.rfc){
      const lbl=src.snr>=10?`${src.id} ${src.snr.toFixed(0)}`:`${src.id}`;
      ctx.font=(src.snr>40?'700 ':'')+`11px ui-monospace,monospace`;
      const tw=ctx.measureText(lbl).width,lx=cx+r+4,ly=cy;
      ctx.fillStyle='rgba(9,11,15,0.7)';ctx.fillRect(lx-1,ly-7,tw+4,14);
      ctx.fillStyle=src.snr>40?'#fff8d0':'#8ab4d0';ctx.textBaseline='middle';ctx.textAlign='left';ctx.fillText(lbl,lx+1,ly);
    }
    ctx.restore();
  });
  const ns=sources.filter(s=>s.streak).length,na=sources.filter(s=>s.arm).length,nr=sources.filter(s=>s.rfc).length;
  const br=sources.reduce((a,b)=>b.snr>a.snr?b:a);
  const ps=document.createElement('div');ps.className='ps';
  ps.innerHTML=`<span>n=<span class="sv">${sources.length}</span></span><span>streak <span class="sv" style="color:#e05050">${ns}</span></span><span>arm <span class="sv" style="color:#c89020">${na}</span></span><span>RFC <span class="sv" style="color:#38c87a">${nr}</span></span><span>peak SNR <span class="sv">${br.snr}&times;</span> id ${br.id} <span style="color:${br.streak?'#e05050':'#4fa8d8'}">${br.streak?'streak':'ok'}</span></span>`;
  panel.appendChild(ps);
  const tt=document.getElementById('tt');
  cw.addEventListener('mousemove',e=>{
    const rect=cv.getBoundingClientRect(),sc2=S/rect.width,mx=(e.clientX-rect.left)*sc2,my=(e.clientY-rect.top)*sc2;
    let best=null,bestD=30;
    sources.forEach(s=>{const d=Math.hypot(s.x*S-mx,s.y*S-my);if(d<bestD){bestD=d;best=s;}});
    if(best){
      const f=[];if(best.streak)f.push('<span class="fs">acis_streak</span>');if(best.arm)f.push('<span class="fa">grating_arm</span>');if(best.rfc)f.push('<span class="fr">RFC in cat_src</span>');if(!f.length)f.push('<span class="fo">no flags</span>');
      tt.innerHTML=`<b>id ${best.id}</b>  SNR ${best.snr.toFixed(1)}<br>${f.join(' · ')}`;
      tt.style.display='block';tt.style.left=(e.clientX+14)+'px';tt.style.top=(e.clientY-10)+'px';
    }else{tt.style.display='none';}
  });
  cw.addEventListener('mouseleave',()=>{tt.style.display='none';});
  return panel;
}
document.getElementById('panels').replaceChildren(...VIZ.map(build));
document.getElementById('note').innerHTML=NOTE_PLACEHOLDER;
</script>"""


def build_html(viz_data: list[dict], title_note: str, result_note: str) -> str:
    compact = json.dumps(viz_data, separators=(",", ":"))
    html = HTML_TEMPLATE
    html = html.replace("TITLE_PLACEHOLDER", title_note)
    html = html.replace("DATA_PLACEHOLDER", compact)
    html = html.replace("NOTE_PLACEHOLDER", json.dumps(result_note))
    return html


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def find_obs_dir(workdir: Path, obsid: int) -> Path:
    """Locate the obsid work directory under workdir (any obs?? subdirectory)."""
    candidates = list(workdir.glob(f"obs*/{obsid}"))
    if not candidates:
        candidates = list(workdir.glob(f"{obsid}"))
    if not candidates:
        raise FileNotFoundError(
            f"No work directory found for obsid {obsid} under {workdir}"
        )
    return candidates[0]


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--obsids", nargs="+", type=int, required=True, help="Obsids to include"
    )
    ap.add_argument(
        "--workdir",
        type=Path,
        default=Path("/Volumes/Black/data/astromon/work"),
        help="Root of pipeline work directory",
    )
    ap.add_argument(
        "--db",
        type=Path,
        default=Path("/Volumes/Black/data/astromon/astromon.h5"),
        help="astromon HDF5 database",
    )
    ap.add_argument(
        "--out", type=Path, default=Path("xray_review.html"), help="Output HTML file"
    )
    ap.add_argument(
        "--factor",
        type=int,
        default=4,
        help="Downsampling factor (default 4, i.e. 720→180)",
    )
    ap.add_argument(
        "--clip-pct",
        type=float,
        default=99.5,
        help="Percentile clip for image scale (default 99.5)",
    )
    args = ap.parse_args()

    viz_data = []
    for obsid in args.obsids:
        obs_dir = find_obs_dir(args.workdir, obsid)
        print(f"obsid {obsid}: {obs_dir}")

        db_srcs = load_db_sources(obsid, args.db)
        rfc_ids = load_cat_src_rfc(obsid, args.db)
        data, ltv1, ltv2, W, H = read_image(obs_dir)

        img_shape, img_ds = (
            (H // args.factor, W // args.factor),
            downsample(data, factor=args.factor, clip_pct=args.clip_pct),
        )
        sources = read_sources(obs_dir, db_srcs, rfc_ids, ltv1, ltv2, W, H)
        streak_polys = read_streak_polys(obs_dir, ltv1, ltv2, W, H)

        bright = max(sources, key=lambda s: s["snr"]) if sources else {}
        print(
            f"  {len(sources)} sources, brightest id={bright.get('id')} "
            f"snr={bright.get('snr')} streak={bright.get('streak')}"
        )

        viz_data.append(
            {
                "obsid": obsid,
                "label": str(obs_dir.parent.name + "/" + obs_dir.name),
                "image": img_ds,
                "img_shape": list(img_shape),
                "sources": sources,
                "streak_polys": streak_polys,
            }
        )

    title_note = (
        "Image: log&#8321;&#8199;(block sum) clipped at 99.5th percentile. "
        "Circle radius ∝ &radic;SNR. Red = <code>acis_streak</code>, "
        "yellow = <code>grating_arm</code>, dashed green ring = RFC in cat_src."
    )
    result_note = (
        "<strong>Reading the display:</strong> "
        "Image uses log&#8321;&#8199;-scaled 4&times;4 block sums clipped at 99.5th percentile "
        "so individual photons are visible. "
        "Source circles are colour-coded: "
        "<span style='color:#e05050'>red = acis_streak</span>, "
        "<span style='color:#c89020'>yellow = grating_arm</span>, "
        "blue = no exclusion flags. "
        "Dashed <span style='color:#38c87a'>green ring</span> = RFC entry in cat_src. "
        "Red shaded region = acis_streak_map polygon. "
        "Hover any source circle for details."
    )
    html = build_html(viz_data, title_note, result_note)
    args.out.write_text(html)
    print(f"Wrote {args.out} ({args.out.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
