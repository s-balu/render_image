# render_sequence — batch interpolated animation driver

## Overview

`render_sequence.sh` wraps `render_image` to produce a **single, self-consistently
numbered sequence of PNG frames** across an arbitrary number of snapshots.  You
point it at a plain-text list of snapshot roots and it handles:

- assigning `snap_a` / `snap_b` / `snap_prev` / `snap_next` per interval
- running `render_image` once per snapshot interval with `-itmax N`
- renaming each interval's frames into one global `frame.0000.png … frame.NNNN.png`
  sequence (no duplicate endpoint frames)
- optionally running multiple intervals in parallel with `-j`

The result is ready to pass straight to `ffmpeg`.

---

## Quick start

```bash
chmod +x render_sequence.sh

# 2-snapshot Hermite, 8 interpolated + 1 endpoint = 9 frames per gap
./render_sequence.sh \
    --binary ./render_image \
    --config base.yml \
    --outdir ./frames \
    --nframes 9 \
    --mode hermite \
    snapshots.txt

# Then assemble:
ffmpeg -framerate 24 -i './frames/frame.%04d.png' \
       -c:v libx264 -crf 18 -pix_fmt yuv420p animation.mp4
```

---

## Snapshot list format

One snapshot root per line; `#` comments and blank lines are skipped:

```
# snapshots.txt
/data/HiSURFS_0089/HiSURFS_0089
/data/HiSURFS_0091/HiSURFS_0091
/data/HiSURFS_0093/HiSURFS_0093
/data/HiSURFS_0095/HiSURFS_0095
```

---

## Options

| Flag | Long form | Default | Description |
|------|-----------|---------|-------------|
| `-b` | `--binary` | _(required)_ | Path to the `render_image` binary |
| `-c` | `--config` | _(required)_ | Base YAML config (all shared settings go here) |
| `-o` | `--outdir` | _(required)_ | Output directory for PNG frames |
| `-n` | `--nframes` | `9` | Frames per interval **including both endpoint frames** |
| `-m` | `--mode` | `hermite` | `hermite` / `catmull` / `catmull4` |
| `-p` | `--prefix` | `frame` | Filename prefix (e.g. `frame` → `frame.0000.png`) |
| `-j` | `--jobs` | `1` | Parallel render_image processes |
| `-d` | `--dry-run` | off | Print commands only, write nothing |
| `-v` | `--verbose` | off | Extra progress output |

---

## Interpolation modes and snapshot roles

### `hermite` — 2-snapshot cubic Hermite (your current default)

Uses `snap_a` and `snap_b`.  The script assigns:

```
interval i:  snap_a = snaps[i]   snap_b = snaps[i+1]
```

Needs at least **2** snapshots.  K = N_snaps − 1 intervals.

### `catmull` — 3-snapshot Catmull-Rom

Uses `snap_prev`, `snap_a`, `snap_b`.  Smoother tangents at interval boundaries.

```
interval i:  snap_prev = snaps[i]   snap_a = snaps[i+1]   snap_b = snaps[i+2]
```

Needs at least **3** snapshots.  K = N_snaps − 2 intervals.

### `catmull4` — 4-snapshot Catmull-Rom

Adds `snap_next` for a fully centred tangent at the right endpoint too.

```
interval i:  snap_prev = snaps[i]   snap_a = snaps[i+1]
             snap_b    = snaps[i+2]  snap_next = snaps[i+3]
```

Needs at least **4** snapshots.  K = N_snaps − 3 intervals.

---

## Frame count arithmetic

```
Total frames = N_intervals × (N − 1) + 1
```

where N = `--nframes`.  The last frame of each interval is identical to the
first frame of the next, so it is written only once.

| snapshots | N (--nframes) | mode | intervals | total frames |
|-----------|--------------|------|-----------|-------------|
| 4 | 9 | hermite | 3 | 25 |
| 4 | 9 | catmull | 2 | 17 |
| 5 | 9 | catmull4 | 2 | 17 |
| 10 | 25 | hermite | 9 | 217 |

---

## Config file strategy

Put **all** rendering settings in `base.yml` (colormap, opacity, view centre,
`isHDF5`, `dark_matter`, `lock_levels`, etc.).  The script overrides only
`snap_a`, `snap_b`, `snap_prev`, `snap_next`, `output`, and `itmax` per
interval on the command line.

> **Tip:** keep `lock_levels: true` in the config so the colour scale is
> auto-set on the very first frame and then frozen for all subsequent frames.
> This ensures consistent brightness across the entire animation.

---

## Parallel execution

Each snapshot interval is completely independent, so `-j 4` will run 4
intervals simultaneously.  Memory use scales with the number of jobs — for
large simulations keep jobs × memory-per-render below available RAM.

---

## ffmpeg assembly

```bash
ffmpeg -framerate 24 \
       -i './frames/frame.%04d.png' \
       -c:v libx264 -crf 18 \
       -pix_fmt yuv420p \
       animation.mp4
```

Adjust `-framerate` to taste.  For smoother playback with many interpolated
frames, 30 or 60 fps works well.
