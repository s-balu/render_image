#!/usr/bin/env bash
# =============================================================================
# render_sequence.sh — batch driver for imbas_renderer interpolated animations
#
# Reads a snapshot list file and renders interpolated frames between consecutive
# pairs (2-snapshot Hermite) or triplets (3/4-snapshot Catmull-Rom), placing
# all frames into a single output directory with self-consistent sequential
# numbering ready for ffmpeg.
#
# SNAPSHOT LIST FILE FORMAT  (one entry per line, # comments allowed):
#   /path/to/snap_0091/snap_0091
#   /path/to/snap_0092/snap_0092
#   /path/to/snap_0093/snap_0093
#   ...
#
# USAGE:
#   render_sequence.sh [OPTIONS] <snapshot_list_file>
#
# OPTIONS:
#   -b, --binary    <path>   Path to imbas_renderer binary        [required]
#   -c, --config    <path>   Base YAML config file                [required]
#   -o, --outdir    <path>   Output directory for PNG frames      [required]
#   -n, --nframes   <N>      Interpolated frames per interval     [default: 9]
#                            (includes both endpoint frames, so N=9 gives
#                             7 interpolated + 2 endpoint frames per gap,
#                             but endpoints are shared between intervals
#                             so only N-1 unique frames are written per gap)
#   -m, --mode      <mode>   Interpolation mode:
#                              hermite   — 2-snapshot cubic Hermite  [default]
#                              catmull   — 3-snapshot Catmull-Rom
#                              catmull4  — 4-snapshot Catmull-Rom
#   -p, --prefix    <str>    Frame filename prefix                [default: frame]
#   -j, --jobs      <N>      Parallel imbas_renderer invocations   [default: 1]
#   -l, --levels  vmin,vmax  Skip probe; use fixed log10 density bounds
#                            e.g. --levels -1.0,4.5
#   -d, --dry-run            Print commands without executing
#   -v, --verbose            Print progress messages
#   -h, --help               Show this help
#
# FRAME NUMBERING:
#   For K snapshot intervals each producing N frames (shared endpoints):
#     Total unique frames = K * (N - 1) + 1
#   Frames are numbered frame.0000.png, frame.0001.png, ...
#
# EXAMPLE — 2-snapshot Hermite, 8 interpolated frames per interval:
#   render_sequence.sh \
#     --binary ./imbas_renderer \
#     --config base.yml \
#     --outdir ./frames \
#     --nframes 9 \
#     --mode hermite \
#     snapshots.txt
#
# EXAMPLE — Catmull-Rom (3-snap), 24 frames per interval:
#   render_sequence.sh -b ./imbas_renderer -c base.yml -o ./frames \
#     -n 25 -m catmull snapshots.txt
#
# THEN ASSEMBLE WITH FFMPEG:
#   ffmpeg -framerate 24 -i ./frames/frame.%04d.png \
#          -c:v libx264 -crf 18 -pix_fmt yuv420p animation.mp4
#
# =============================================================================

set -euo pipefail

# ─── Runtime compatibility check ─────────────────────────────────────────────
BASH_MAJOR="${BASH_VERSINFO[0]}"
BASH_MINOR="${BASH_VERSINFO[1]}"
if [[ $BASH_MAJOR -lt 3 || ( $BASH_MAJOR -eq 3 && $BASH_MINOR -lt 2 ) ]]; then
    echo "ERROR: bash 3.2 or newer required (you have $BASH_VERSION)" >&2
    exit 1
fi
if [[ $BASH_MAJOR -lt 4 ]]; then
    echo "NOTE: Running on bash $BASH_VERSION (macOS system shell)." \
         "Parallel -j mode requires bash 4.3+ for best results." >&2
fi

# ─── Defaults ────────────────────────────────────────────────────────────────
BINARY=""
CONFIG=""
OUTDIR=""
NFRAMES=9
MODE="hermite"
PREFIX="frame"
JOBS=1
DRY_RUN=0
VERBOSE=0
LEVELS=""       # empty = run probe pass; "vmin,vmax" = skip probe

# ─── Argument parsing ────────────────────────────────────────────────────────
usage() {
    sed -n '/#.*USAGE/,/^# ====/p' "$0" | grep '^#' | sed 's/^# \?//'
    exit 0
}

die() { echo "ERROR: $*" >&2; exit 1; }
info() { [[ $VERBOSE -eq 1 ]] && echo "[INFO] $*" || true; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        -b|--binary)   BINARY="$2";  shift 2 ;;
        -c|--config)   CONFIG="$2";  shift 2 ;;
        -o|--outdir)   OUTDIR="$2";  shift 2 ;;
        -n|--nframes)  NFRAMES="$2"; shift 2 ;;
        -m|--mode)     MODE="$2";    shift 2 ;;
        -p|--prefix)   PREFIX="$2";  shift 2 ;;
        -j|--jobs)     JOBS="$2";    shift 2 ;;
        -l|--levels)   LEVELS="$2"; shift 2 ;;
        -d|--dry-run)  DRY_RUN=1;   shift   ;;
        -v|--verbose)  VERBOSE=1;   shift   ;;
        -h|--help)     usage ;;
        -*)            die "Unknown option: $1" ;;
        *)             SNAP_LIST="$1"; shift ;;
    esac
done

# ─── Validate inputs ─────────────────────────────────────────────────────────
[[ -z "${SNAP_LIST:-}" ]] && die "No snapshot list file provided."
[[ -z "$BINARY" ]]        && die "--binary is required."
[[ -z "$CONFIG" ]]        && die "--config is required."
[[ -z "$OUTDIR" ]]        && die "--outdir is required."

[[ -f "$SNAP_LIST" ]] || die "Snapshot list not found: $SNAP_LIST"
[[ -f "$BINARY"    ]] || die "imbas_renderer binary not found: $BINARY"
[[ -x "$BINARY"    ]] || die "imbas_renderer binary is not executable: $BINARY"
[[ -f "$CONFIG"    ]] || die "Config file not found: $CONFIG"

[[ "$NFRAMES" -ge 2 ]] || die "--nframes must be >= 2."
[[ "$JOBS"    -ge 1 ]] || die "--jobs must be >= 1."

case "$MODE" in
    hermite|catmull|catmull4) ;;
    *) die "Unknown mode '$MODE'. Choose: hermite | catmull | catmull4." ;;
esac

# ─── Read snapshot list (strip comments and blank lines) ─────────────────────
# Use a while-read loop instead of mapfile for bash 3.2 compatibility (macOS).
SNAPS=()
while IFS= read -r line; do
    SNAPS+=("$line")
done < <(grep -v '^\s*#' "$SNAP_LIST" | grep -v '^\s*$')
N_SNAPS=${#SNAPS[@]}

info "Snapshot list: $SNAP_LIST  ($N_SNAPS entries)"

# Minimum snapshots needed per mode
case "$MODE" in
    hermite)  MIN_SNAPS=2 ;;
    catmull)  MIN_SNAPS=3 ;;
    catmull4) MIN_SNAPS=4 ;;
esac

[[ $N_SNAPS -ge $MIN_SNAPS ]] || \
    die "Mode '$MODE' needs at least $MIN_SNAPS snapshots; list has $N_SNAPS."

# ─── Prepare output directory ────────────────────────────────────────────────
if [[ $DRY_RUN -eq 0 ]]; then
    mkdir -p "$OUTDIR"
fi

# ─── Compute interval structure ──────────────────────────────────────────────
# Number of snapshot intervals (gaps between consecutive pair endpoints):
#   hermite:  K = N_SNAPS - 1   (each interval uses snaps [i, i+1])
#   catmull:  K = N_SNAPS - 2   (each interval uses snaps [i-1, i, i+1])
#   catmull4: K = N_SNAPS - 3   (each interval uses snaps [i-1, i, i+1, i+2])
case "$MODE" in
    hermite)  N_INTERVALS=$(( N_SNAPS - 1 )) ;;
    catmull)  N_INTERVALS=$(( N_SNAPS - 2 )) ;;
    catmull4) N_INTERVALS=$(( N_SNAPS - 3 )) ;;
esac

# Total unique frames:
#   Each interval contributes (NFRAMES - 1) unique frames
#   (the last frame of interval i == the first of interval i+1, so shared).
#   Plus the very first frame of the entire sequence.
TOTAL_FRAMES=$(( N_INTERVALS * (NFRAMES - 1) + 1 ))

echo "Mode:            $MODE"
echo "Snapshot count:  $N_SNAPS"
echo "Intervals:       $N_INTERVALS"
echo "Frames/interval: $NFRAMES  (including shared endpoints)"
echo "Total frames:    $TOTAL_FRAMES"
echo "Output dir:      $OUTDIR"
echo ""

# ─── Temporary per-interval directory helper ─────────────────────────────────
# imbas_renderer writes frames as <root>.NNNN.png with its own counter
# starting at 0000 for each invocation.  We redirect each interval to
# a temp sub-directory, then rename/move the PNGs into OUTDIR using the
# global frame counter.

TMPBASE=$(mktemp -d "${OUTDIR}/.rtmp_XXXXXX")
trap 'rm -rf "$TMPBASE"' EXIT

# ─── Level probe pass ────────────────────────────────────────────────────────
# Run imbas_renderer once with -itmax 1 on the first interval to let it
# auto-level on frame 0, then capture the vmin/vmax it prints.
# All subsequent interval renders use -no_auto_levels -vmin X -vmax Y
# so the colour scale is identical across the entire sequence.
#
# Skip if the user supplied --levels vmin,vmax explicitly.

LEVEL_FLAGS=""   # will be set to "-no_auto_levels -vmin X -vmax Y"

if [[ -n "$LEVELS" ]]; then
    # User-supplied levels — parse and use directly
    VMIN="${LEVELS%%,*}"
    VMAX="${LEVELS##*,}"
    LEVEL_FLAGS="-no_auto_levels -vmin $VMIN -vmax $VMAX"
    echo "Levels (user-supplied): vmin=$VMIN  vmax=$VMAX"
else
    echo "Running level probe on first interval..."

    # Build the probe command using the same snapshot flags as interval 0
    case "$MODE" in
        hermite)
            PROBE_FLAGS="-snap_a \"${SNAPS[0]}\" -snap_b \"${SNAPS[1]}\""
            ;;
        catmull)
            PROBE_FLAGS="-snap_prev \"${SNAPS[0]}\" -snap_a \"${SNAPS[1]}\" -snap_b \"${SNAPS[2]}\""
            ;;
        catmull4)
            PROBE_FLAGS="-snap_prev \"${SNAPS[0]}\" -snap_a \"${SNAPS[1]}\" -snap_b \"${SNAPS[2]}\" -snap_next \"${SNAPS[3]}\""
            ;;
    esac

    PROBE_TMPDIR="${TMPBASE}/probe"
    mkdir -p "$PROBE_TMPDIR"
    PROBE_CMD="\"$BINARY\" -config \"$CONFIG\" $PROBE_FLAGS -output \"${PROBE_TMPDIR}/probe\" -itmax 1"

    info "Probe command: $PROBE_CMD"

    if [[ $DRY_RUN -eq 0 ]]; then
        # Capture stdout; imbas_renderer prints the locked levels as:
        #   "Levels locked (pre-inpaint): vmin=X.XXX vmax=X.XXX"
        PROBE_LOG="${PROBE_TMPDIR}/probe.log"
        eval "$PROBE_CMD" 2>&1 | tee "$PROBE_LOG"

        LEVELS_LINE=$(grep "Levels locked" "$PROBE_LOG" | tail -1)
        if [[ -z "$LEVELS_LINE" ]]; then
            echo "WARNING: Could not find 'Levels locked' in probe output." >&2
            echo "         Falling back to per-interval auto-levels (flicker may occur)." >&2
            echo "         Check $PROBE_LOG for details." >&2
        else
            VMIN=$(echo "$LEVELS_LINE" | grep -oE 'vmin=[^ ]+' | cut -d= -f2)
            VMAX=$(echo "$LEVELS_LINE" | grep -oE 'vmax=[^ ]+' | cut -d= -f2)
            LEVEL_FLAGS="-no_auto_levels -vmin $VMIN -vmax $VMAX"
            echo "Levels locked globally: vmin=$VMIN  vmax=$VMAX"
        fi
    else
        echo "[DRY-RUN] $PROBE_CMD"
        echo "[DRY-RUN] Would parse 'Levels locked' line and pin vmin/vmax globally."
        LEVEL_FLAGS="-no_auto_levels -vmin <probed> -vmax <probed>"
    fi
fi
echo ""
# Each entry: "<interval_index> <first_global_frame> <snap_prev> <snap_a> <snap_b> <snap_next>"
# We emit them in order and later run them (possibly in parallel with xargs).

CMDS_FILE=$(mktemp)

GLOBAL_FRAME=0  # next frame index to be written

for (( i=0; i<N_INTERVALS; i++ )); do

    # Frame range for this interval:
    #   first frame: GLOBAL_FRAME
    #   last  frame: GLOBAL_FRAME + NFRAMES - 1  (but last == first of next interval)
    # We write frames [GLOBAL_FRAME .. GLOBAL_FRAME + NFRAMES - 2] for this interval
    # (skip the last one to avoid duplication), EXCEPT for the very last interval
    # which writes all NFRAMES frames including the true endpoint.

    INTERVAL_TMPDIR="${TMPBASE}/interval_$(printf '%04d' $i)"
    INTERVAL_ROOT="${INTERVAL_TMPDIR}/frame"

    # Assign snapshot roles per mode
    case "$MODE" in
        hermite)
            SNAP_A="${SNAPS[$i]}"
            SNAP_B="${SNAPS[$((i+1))]}"
            EXTRA_FLAGS="-snap_a \"$SNAP_A\" -snap_b \"$SNAP_B\""
            ;;
        catmull)
            # P=i, A=i+1, B=i+2
            SNAP_PREV="${SNAPS[$i]}"
            SNAP_A="${SNAPS[$((i+1))]}"
            SNAP_B="${SNAPS[$((i+2))]}"
            EXTRA_FLAGS="-snap_prev \"$SNAP_PREV\" -snap_a \"$SNAP_A\" -snap_b \"$SNAP_B\""
            ;;
        catmull4)
            # P=i, A=i+1, B=i+2, N=i+3
            SNAP_PREV="${SNAPS[$i]}"
            SNAP_A="${SNAPS[$((i+1))]}"
            SNAP_B="${SNAPS[$((i+2))]}"
            SNAP_NEXT="${SNAPS[$((i+3))]}"
            EXTRA_FLAGS="-snap_prev \"$SNAP_PREV\" -snap_a \"$SNAP_A\" -snap_b \"$SNAP_B\" -snap_next \"$SNAP_NEXT\""
            ;;
    esac

    # Frames to keep from this render: all if last interval, else skip last
    if [[ $i -lt $((N_INTERVALS - 1)) ]]; then
        KEEP_FRAMES=$(( NFRAMES - 1 ))
    else
        KEEP_FRAMES=$NFRAMES
    fi

    # Write command to the command file:
    #   Fields: interval_dir | global_start | keep_frames | render_command
    CMD="\"$BINARY\" -config \"$CONFIG\" $EXTRA_FLAGS $LEVEL_FLAGS -output \"$INTERVAL_ROOT\" -itmax $NFRAMES"
    echo "${INTERVAL_TMPDIR}|${GLOBAL_FRAME}|${KEEP_FRAMES}|${CMD}" >> "$CMDS_FILE"

    GLOBAL_FRAME=$(( GLOBAL_FRAME + KEEP_FRAMES ))
done

# ─── Execute (serial or parallel) ────────────────────────────────────────────
run_interval() {
    local line="$1"
    IFS='|' read -r INTERVAL_DIR GSTART KEEP RENDER_CMD <<< "$line"

    mkdir -p "$INTERVAL_DIR"

    info "Running: $RENDER_CMD"

    if [[ $DRY_RUN -eq 0 ]]; then
        eval "$RENDER_CMD" || { echo "ERROR: imbas_renderer failed for $INTERVAL_DIR" >&2; return 1; }
    else
        echo "[DRY-RUN] $RENDER_CMD"
    fi

    # Rename interval frames → global frame numbers
    for (( f=0; f<KEEP; f++ )); do
        LOCAL_NAME=$(printf '%s.%04d.png' "$INTERVAL_DIR/frame" "$f")
        GLOBAL_NAME=$(printf '%s/%s.%04d.png' "$OUTDIR" "$PREFIX" "$(( GSTART + f ))")
        if [[ $DRY_RUN -eq 0 ]]; then
            [[ -f "$LOCAL_NAME" ]] || { echo "ERROR: expected frame not found: $LOCAL_NAME" >&2; return 1; }
            mv "$LOCAL_NAME" "$GLOBAL_NAME"
        else
            echo "[DRY-RUN] mv $LOCAL_NAME -> $GLOBAL_NAME"
        fi
    done
}

if [[ $JOBS -le 1 ]]; then
    while IFS= read -r line; do
        run_interval "$line"
    done < "$CMDS_FILE"
else
    # Portable parallel execution using background jobs + a semaphore counter.
    # Avoids export -f / xargs -P which behave differently on macOS bash 3.2.
    echo "[INFO] Running up to $JOBS parallel jobs..."
    RUNNING=0
    while IFS= read -r line; do
        run_interval "$line" &
        RUNNING=$(( RUNNING + 1 ))
        if [[ $RUNNING -ge $JOBS ]]; then
            wait -n 2>/dev/null || wait   # wait -n needs bash 4.3+; fall back
            RUNNING=$(( RUNNING - 1 ))
        fi
    done < "$CMDS_FILE"
    wait   # drain any remaining background jobs
fi

rm -f "$CMDS_FILE"

# ─── Summary ─────────────────────────────────────────────────────────────────
if [[ $DRY_RUN -eq 0 ]]; then
    ACTUAL=$(ls "$OUTDIR"/${PREFIX}.*.png 2>/dev/null | wc -l)
    echo ""
    echo "Done. $ACTUAL frames written to $OUTDIR/"
    echo ""
    echo "Assemble with ffmpeg:"
    echo "  ffmpeg -framerate 24 -i '${OUTDIR}/${PREFIX}.%04d.png' \\"
    echo "         -c:v libx264 -crf 18 -pix_fmt yuv420p animation.mp4"
else
    echo ""
    echo "[DRY-RUN complete — no files written]"
fi
