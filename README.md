# Real-time afterdischarge (AD) detection — prelim data pipeline

Self-contained sub-project for the stim grant prelim data with Mike Beauchamp.
Goal: show that 50 Hz stim sessions can be processed **causally / in real time**
to detect afterdischarges, and quantify accuracy (sensitivity, specificity,
false-positive rate) against human ground truth, plus AD onset-timing error.

It reuses the algorithms from the original `find_ad_fcn.m` but refactors the
detector into a true streaming (chunk-at-a-time) form, and adds a
detector-independent ground-truth GUI so we can compute sensitivity/specificity
(not just PPV). This is a standalone project, separate from the prior work.

## Project layout

The **git repo is this `scripts/` folder** (code + `hfs_sessions.csv` + README).
Generated data lives in the parent workspace, *outside* the repo:

```
<workspace>/                 (not tracked)
├── scripts/                 <- THE GIT REPO
│   ├── *.m, helpers/
│   ├── hfs_sessions.csv     input list of ieeg datasets / time windows
│   ├── rt_local_config.m    per-machine paths (you create this; gitignored)
│   ├── .gitignore
│   └── README.md
├── clips/                   per-stim clips — the working data (generated)
├── results/                 detections, evaluation output (generated)
├── ground_truth/            human AD annotations: clip_gt.csv (generated)
└── data/                    optional full-session .mat
```

**Clips are the unit of everything.** Rather than storing whole sessions (all
~130 channels for thousands of seconds ≈ multi-GB each), the pipeline saves one
short clip per stimulation containing only the stim electrode's contacts
(~7 MB each). Every downstream step — annotation, detection, evaluation —
operates on these clips.

`rt_paths.m` locates itself, treats `scripts/` as the repo root and its parent
as the workspace where generated data is written, and reads per-machine paths
from `rt_local_config.m` (falling back to `seizure_termination_paths` on the
original machine).

## Files (in scripts/)

| File | Role |
|------|------|
| `rt_paths.m` | Self-locating paths/config. Generated output goes to the project's `clips/`, `results/`, `ground_truth/`. |
| `test_stim_detection.m` | **Run first.** Streams 5-min chunks of one session, detects the first few stims, checks the detected pair against the "Closed relay" annotations, and plots them. |
| `export_stim_clips.m` | **Main data step.** Streams to find every stim (no whole-session load), then saves one short clip per stimulation (stim-electrode channels only) to `clips/` as `.mat` + `.edf`, plus `clip_index.csv`. |
| `find_stim_events.m` | 50 Hz stim detection on a data chunk → per-event channel pair + on/off times. |
| `find_session_stims.m` *(helper)* | Streams a whole session in chunks, runs `find_stim_events`, merges boundary splits/re-detections. |
| `annotate_ground_truth.m` | GUI: per **clip**, mark AD yes/no and click onset/offset. Writes `ground_truth/clip_gt.csv`. |
| `rt_ad_detector_init.m` / `rt_ad_detector_step.m` | Streaming detector. Feed it blocks of 100–1000 ms; it detects stim + ADs causally with persistent state. |
| `run_realtime_sim.m` | Replay each **clip** block-by-block; record detections + per-block CPU latency + real-time factor. |
| `evaluate_detector.m` | Match per-clip detections to per-clip ground truth → sensitivity/specificity/FPR/PPV + onset error. |
| `sweep_block_sizes.m` | Run sim+eval across block sizes (100/250/500/1000 ms) for a single summary table. |
| `run_pipeline.m` | Convenience driver for the machine steps (sim + eval). |
| `export_sessions.m` *(optional)* | Pull a whole session to `data/` as `.mat`. Not part of the main flow; only if you want full sessions. |
| `helpers/` | Shared functions (bipolar pairs, high-pass filter, power, channel selection, label cleanup). |

## How to run (in MATLAB)

```matlab
cd .../seizure_termination/ad_realtime_prelim/scripts

% 0) Sanity-check stim detection on the first few stims of session 1
test_stim_detection(1)

% 1) Main data step: save one clip per stimulation (stim-electrode channels)
export_stim_clips          % all sessions; export_stim_clips(1) for just session 1

% 2) Human ground truth (detector-independent), one clip at a time
annotate_ground_truth      % mark AD y/n + onset/offset per clip

% 3) Streaming detection + timing over clips (pick a block size in ms)
run_realtime_sim(250)

% 4) Accuracy stats vs ground truth
evaluate_detector(250)

% Or sweep block sizes for the prelim table:
sweep_block_sizes([100 250 500 1000])
```

## Running on another machine / server

The repo is the `scripts/` folder and is self-contained and portable. To set it
up on a new server, clone it into a workspace folder (generated data will be
written one level up, alongside the clone):

```bash
mkdir ad_realtime_prelim && cd ad_realtime_prelim
git clone <your-repo-url> scripts
```

Then in MATLAB:

```matlab
cd ad_realtime_prelim/scripts
copyfile rt_local_config_template.m rt_local_config.m   % then edit it
```

Edit `rt_local_config.m` with three values for that machine: the IEEG MATLAB
toolbox folder, your ieeg.org `.bin` password file, and your ieeg.org login.
That file is gitignored, so credentials are never committed. `rt_paths` picks
it up automatically (and falls back to the original lab paths on the original
machine if no local config exists).

What is and isn't in the repo:

- **In git:** all code (including the bundled `download_ieeg_data`),
  `hfs_sessions.csv`, this README, `.gitignore`, the config *template*.
- **Not in git (regenerated or copied):** `clips/`, `results/`, `ground_truth/`,
  `data/` (all in the parent workspace), your `rt_local_config.m`, and any `.bin`.

So on a new server you either (a) re-create the clips with `export_stim_clips`
(needs ieeg.org access there), or (b) copy the `clips/` and `ground_truth/`
folders over manually if you just want to run detection/evaluation offline.

**Requirements:** MATLAB R2021a+ (for `edfwrite`; older versions still work if
you skip EDF and use the `.mat` clips), Signal Processing Toolbox (`bandpower`,
`bandstop`, `edfwrite`), Statistics Toolbox (`prctile`), and the IEEG MATLAB
toolbox + ieeg.org credentials (only needed for the download steps; detection
and evaluation on existing clips need neither).

## Notes / design choices

- **Trial unit = clip = one stimulation.** Sensitivity = fraction of true-AD
  clips the detector catches; specificity = fraction of no-AD clips it
  correctly leaves alone. This is what lets us report spec/FPR, which a
  PPV-only review of detections cannot.
- **Clips hold only the stim electrode.** That is exactly the set of channels
  the AD detector inspects (contacts on the stimulated electrode), so nothing
  is lost, and each clip is tiny. The detector still auto-detects the stim pair
  and on/off inside the clip (the bipolar pairs are present).
- **Truly causal detector.** `rt_ad_detector_step` carries a residual sample
  tail and all filter/buffer/baseline state across blocks, so results are
  identical regardless of block size — only *latency* changes with block size.
  Worst-case per-block CPU time vs. the block budget is reported to argue
  real-time feasibility.
- **Each clip is saved as both `.mat` (lossless) and `.edf` (portable).** The
  `.mat` keeps exact labels, fs, NaNs and absolute times; the EDF is zero-filled
  for NaN and padded to whole seconds.
- Parameters are the values tuned in the original `find_ad_fcn`; override any
  by passing a `params` struct to `rt_ad_detector_init` / `find_stim_events` /
  `export_stim_clips`.
- One intentional fix vs. the original: excluded channels are masked
  column-wise (`dataChunk(:,exclude)=nan`) rather than via linear indexing.
