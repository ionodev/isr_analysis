# To do

Open work in this codebase, in priority order. Most items come from `tbd`/`TBD`
comments left in the source; the file and line of each is given so the comment
can be read in context.

Items 18-20 come from comparing the measured system temperature against a
440 MHz sky model; see memo 7 and the `tsys-sky-noise` branch.

---

## 1. Antenna selection was interpolated as if it were continuous — DONE

`millstone_radar_state.py:get_antenna_select`

This item was originally written as "the antenna control metadata and the power
meter contradict each other, so nothing is analysed at all". An audit of the
full 49 hour eclipse2024 recording shows that was wrong on both counts. The two
sources agree at 96 to 97 per cent of pulses, and the pipeline was already
accepting 93 per cent of them. The original claim came from one sample at the
very start of the recording, generalised without checking.

The real defect was in the interpolation of the antenna state, which is a step
function switching roughly every five minutes, 595 times in the recording:

- `interp1d` with its default linear kind makes the value sweep continuously
  between -1 and +1 across every switch, so `tx_ant(t) <= -0.99` and
  `>= 0.99` are both false for a slice of each transition. Those pulses were
  discarded: 93.1 per cent accepted instead of 96.3.
- The first sample was moved 24 hours earlier and the last 24 hours later to
  extend the range. Linearly interpolated, that ramps from the first value
  towards the second across the whole 24 hours, so times before the first event
  returned approximately the *second* event's value: at the start of
  eclipse2024, zenith while MISA was transmitting at 1.14 MW.
- `get_antenna_select` returned `(rx_sel, tx_sel)` while all three callers
  unpack `tx_ant, rx_ant`. Harmless as long as both are tested for the same
  value, but wrong for anything that checks them separately.

Fixed with `kind="previous"`, `fill_value` holding the first and last value,
and the return order corrected. Every sampled pulse now evaluates to exactly
+-1, acceptance goes from 93.1 to 96.3 per cent, all of the gain on `misa-l`,
and the start of a recording is attributed to the right antenna.

No cross-source arbitration is needed: where the step model is definite, the
power meter confirms it at 96.0 per cent (MISA) and 97.3 per cent (zenith).

## 2. Uncertainties do not match the parameters they are attached to — DONE

Closed on branch `plasma-parameter-uncertainties`, see memo 6. Two indexing
defects in the covariance were repaired in `fit_lpi.py` and `fit_ionline.py`,
the covariance is now returned and stored per range gate, `dTe` is propagated
and written, `dne` is corrected, `plot_pp.py` consumes `dTe`, and a Monte Carlo
gives a median empirical/reported ratio of 0.989. The same work was carried
into `fit_lp.py`.

Left open by it, and listed below: the variance convention of item 2b, and the
space-object threshold of item 2c.

## 2b. The variance convention between the inversion and the fit — DONE

`outlier_lpi.py:559`, `fit_lpi.py:298`, `fit_ionline.py:293`

`outlier_lpi.py` reports `acfs_var` from a complex least squares, so it is
`E|xhat-x|^2`, a total complex variance. `fit_acf` builds its weight matrix for
the real and imaginary components separately, each carrying half of that, but
used `var` directly, assuming twice the variance it should and overstating every
uncertainty by sqrt(2). Fixed: `std = sqrt(var_scale*var/2)` at all six sites
across the two files.

Measured rather than argued. The scatter of the ACF estimate between
consecutive 10 s files, detrended for ionospheric drift, against the `acfs_var`
those same files report:

| filtering | acfs_var / observed | reported uncertainty vs observed scatter |
|---|---|---|
| production, 18 kHz, filter 100 | 1.67 | 2.58x |
| wide, 100 kHz, filter 20 | 1.21 | 2.20x |

The calibration of `acfs_var` depends on the processing parameters, so no single
number for it is valid across configurations.

`var_scale = 2` is kept, and the evidence says it should be. It is not
compensating for the convention error: it stands in for lag-to-lag correlation
that the diagonal weight matrix ignores. Measured correlation at unit lag
offset is 0.66 for the production filtering, an effective independent fraction
of 0.34, justifying an inflation near 2.9; at 100 kHz with filter_len 20 it is
0.06 and 1.4 would do. Two sits between the two and under-compensates in
production. Removing it would have made the uncertainties worse.

With both, reported uncertainties land within about 7 per cent of the observed
scatter in the production configuration, against 2.6 times before.

Two hypotheses were tested and rejected along the way. The global median floor
on `sigma_lp_est` makes no measurable difference: run against the same window,
the ratio is 1.214 with it and 1.204 with a per-measurement floor. And the
1.67 measured on the archive is a processing-parameter effect, not a code
defect; the same code on the same window with a wider pass band gives 1.21.

Still open: the principled fix is a non-diagonal measurement covariance rather
than a scalar inflation, which would remove the configuration dependence of
`var_scale`.

## 2c. The space-object threshold was tuned against an inflated uncertainty — DONE

`fit_lpi.py:709`, `fit_ionline.py:631`

The detector flagged a gate when `gres[0] < 300.0 and gsigma[0] < 100`: a
Doppler width too narrow for thermal plasma, and an uncertainty small enough to
believe it. The first cut is physical, the second is a confidence cut written as
an absolute number, so it moved when the covariance repair shrank the
uncertainties. The detector silently became more sensitive.

Measured over 30 files and 184 range gates, 5520 fits with both the original
and the repaired `fit_gaussian`:

| cut | original cov | repaired cov | overlap |
|---|---|---|---|
| `width<300 and dwidth<100` | 1 | 2 | 1 |
| `width + 1*dwidth < 300` | 13 | 37 | 13 |
| `width + 2*dwidth < 300` | 3 | 3 | 3 |
| `width + 3*dwidth < 300` | 1 | 1 | 1 |

Replaced with a significance cut, `gres[0] + k*gsigma[0] < 300.0`, which asks
that the width be below the physical threshold at k sigma. k is 3: it was 2
when this was first set, and the variance convention fix of item 2b later
shrank `gsigma` by sqrt(2), so 2*sqrt(2) preserves the strictness it was tuned
at and selects the same bright gates. It needs no retuning when the error
bars change, which matters because item 2b may change them again.

Why two sigma: it is the loosest cut that is stable across the repair, every
gate it selects is bright (3.7 to 6.8 times the profile median, against 0.9 for
all gates), and it contains everything the previous test found. One sigma
admits a gate with negative power. The median fitted width is 493 m/s with a
median uncertainty of 1315 m/s, so the width alone is undetermined at most
gates and the confidence cut does nearly all the work, which is why it should
not have been an absolute constant.

End to end on one 300 s period: 2 detections, at 881 and 890 km, two range
gates apart and both bright, the same object smeared across neighbouring gates.

Not validated against an independent catalogue: "bright and spectrally narrow"
is strong evidence of a hard target, but confirming the objects would mean
running the CFAR detector over the same window.

## 3. Fourth moments of lagged products — MEASURED, implementation deferred

`outlier_lpi.py:543`

The comment doubts that `<(m_t m*_{t+tau})(m_t m*_{t+tau})>` vanishes, "as it
might not be zero when snr is high!!!". For a Gaussian process that pseudo
variance is `R(tau)^2`, so the error on the ACF estimate is not circular: its
real and imaginary parts have unequal variances, correlated, along a direction
set by the ACF phase. One `acfs_var` per range and lag cannot describe that,
and `fit_acf` weights both components equally.

Measured, on 20 consecutive files, 8602 range-lag cells. A per-cell test cannot
see it: with 18 degrees of freedom the estimator's own floor is 0.236 and the
median circularity is 0.303. Pooling the cells with the phase alignment the
`R^2` prediction implies pulls it out:

| cells | pseudo variance / variance | significance |
|---|---|---|
| all (8602) | +0.0236 +- 0.0028 | 8.4 sigma |
| weak signal (4301) | +0.0045 +- 0.0038 | 1.2 sigma |
| strong (3440) | +0.0141 +- 0.0042 | 3.3 sigma |
| strongest (861) | +0.1568 +- 0.0117 | 13.4 sigma |

The imaginary part is +0.002, consistent with zero, as the prediction requires.
The effect is real, it grows with signal exactly as predicted, and it is
confined to the bright gates: 16 per cent of the variance at the strongest
cells, an error ellipse eccentric by about 17 per cent, under 1.5 per cent at
the median cell.

Not implemented. The fix is a widely linear (augmented complex) least squares
in the inversion, storing the pseudo variance beside `acfs_var`, and a
non-circular weighting in `fit_acf`: a rewrite of both the inversion and the
fit, which is what the comment means by "a slightly different method". For an
effect that perturbs a minority of gates by 17 per cent in one direction of an
error ellipse the fit does not currently represent at all.

This is the first item where the measurement argues against doing the work. It
is recorded so that the judgement can be revisited if a science case needs the
bright gates' error ellipses, or if the inversion is rewritten for another
reason.

## 4. Jammer notch is hardcoded to one dataset — DONE

`fit_lp.py:380`

`fit_spectra` now takes `notch_bands_hz` as `[[low,high],...]` in Hz and
`fit_bandwidth_hz`, both passed through `run_analysis.py`, and logs how many
frequency bins survive on every run so the exclusion is visible rather than
silent.

The default keeps the historical 20-28 kHz band rather than defaulting to
empty, against the original plan, because the measurement argued for it. A
first look said there was no jammer in the eclipse data: the band sits at 1.11
times the surrounding median on `misa-l`, 0.99 on `zenith-l`. That test was too
blunt. Fitting with and without the notch on the same 25 files shows the band
matters:

| | with notch | without |
|---|---|---|
| median dTe | 28.31 | 33.61 |
| median dTi | 41.24 | 50.21 |
| median dvi | 6.46 | 8.23 |
| bins fitted | 376 | 409 |

Every uncertainty is 24 per cent larger without the notch, and vi moves by 5
per cent, because the residual those bins contribute is what `fit_spec` uses to
estimate sigma. Looking again bin by bin, there is a broad feature over roughly
24-28 kHz peaking at 1.11 times the 30-50 kHz baseline at 26 kHz: weak in a
median, strong where it is strong.

Interference remains a property of a recording, not of the instrument, so this
should be reviewed per campaign. That is now possible without editing source.

## 4b. `fit_lp.py` ran analyses on import — DONE

`fit_lp.py:690`

The driver block at the bottom, `dirs=[...]` followed by three `fit_spectra`
calls against hardcoded `/media/j/...` paths, had no `__main__` guard, so it
ran whenever the module was imported. `run_analysis.py` imports it to reach
`fit_spectra`, so every long-pulse run began by attempting three analyses
nobody asked for. They failed fast only because those paths do not exist on
this machine. Guarded; import now takes 1.4 s and does nothing.

## 5. Space object removal window was hardcoded — DONE

`fit_lp.py:space_object_halfwidth`

A hard target's echo spreads over the transmit pulse, so the contaminated
region is one pulse length either side of the detection. That was
`for j in range(-17,17)`, correct only for the 30 us gate and 480 us pulse it
was written for. It is now `ceil(pulse_length_us / range_gate_us)`, with the
gate size and the mode read from the file's own `rg` and `mode` fields and a
`pulse_length_us` override passed through `run_analysis.py`.

| | derived | was |
|---|---|---|
| mode 300, 30 us gates | ±16 | 17 |
| mode 800, 30 us gates | ±67 | 17 |
| mode 300, 60 us gates | ±8 | 17 |

Mode 300 reproduces the old constant, which is the check that it was right for
the case it was written for. Mode 800 is the payoff: the low-elevation scan
transmits a 2 ms pulse, so the fixed 17 under-blanked by a factor of four and
left contamination in about 100 gates around every detected object. At a 60 us
gate it over-blanked by double.

The pulse lengths are measured from `tx-h` at half maximum rather than inferred
from the transmit gates, which carry guard time: mode 300 rises at sample 103
and falls at 581, giving 479 us against a 569 us gate; mode 800 at 103 and
2100, giving 1998 us against 2102. The 479 matches the 480 us that appears
throughout the codebase.

## 6. Output metadata written as constants — DONE

`fit_lp.py:713`

`range_avg_limits_km` was written as `[0,1500]` and `range_avg_window_km` from a
fixed 480 us pulse, whatever the run actually did. Both now describe the run:

| | mode 300 | mode 800 | was |
|---|---|---|---|
| range resolution | 71.8 km | 299.5 km | 72.0 km |

`range_avg_window_km` is the pulse length expressed in range, reusing the
measured pulse lengths of item 5, and `range_avg_limits_km` is the span `ridx`
selects. The long pulse does no range averaging, so the "window" is its
inherent range resolution; the comment in the source says so.

For mode 300 the old constant was nearly right, 72.0 against 71.8 km, because
it was written for that mode. For mode 800 it was wrong by a factor of four.
The same pattern as items 4 and 5: correct for its own case, silently wrong for
the other.

Not decorative. `plot_pp.py:49` reads both and passes them into the merged
`ppar-*.h5`, so anything downstream describing the resolution of a long-pulse
product has been reading 72 km regardless of mode.

## 7. Integration periods have no end time — DONE

`outlier_lpi.py:748`, `fit_lpi.py:666`

`outlier_lpi.py` stored only `i0`, so `fit_lpi.py` had to take the end of a
group from the last file's *start*, and said so in a comment. Every recorded
span was short by one file length.

`outlier_lpi.py` now writes `t0` and `t1`, the nominal integration window,
alongside `i0`, which is kept because existing files and readers use it.
`fit_lpi.py` reads `t1` when present and falls back to `i0` otherwise, so older
files still load with their old, short span.

Verified on freshly written files: each carries `t1 - t0 = 10.0 s`, matching
`avg_dur`, and for a group of six the span is now 60.0 s where the old code
reported 50.0 s.

The span is the nominal window rather than the timestamp of the last surviving
pulse, which would vary with how many pulses selection rejected and would make
adjacent periods appear to have gaps between them.

## 7b. `fit_lp.py` and `fit_ionline.py` do not import — DONE

`fit_lp.py:18`, `fit_ionline.py:11`

Both carried `import isr_spec.il_interp as il`, which treats a flat module as a
package, plus two further calls into an `il_interp` API that no longer exists:
`ilint(fname=...)` and `getspec(mol_frac=...)`. Neither module could be
imported, so the `fit_lp` pipeline step would have crashed immediately if it
had been enabled, and `fit_ionline.py` was dead code.

Fixed on the same branch: both now use the flat
import and the lazy `_init_tables` pattern of `fit_lpi.py`, and the renamed
`ion1_frac` keyword. `run_analysis.py` passes `radar_freq_hz` and `table_dir`
through to `fit_lp` as it already did for `fit_lpi`.

Also fixed: `fit_spectra` wrote its `pp-*.h5` to `dirname` while checking for
existing output under `output_base`, so running with a separate output
directory wrote results into the raw data archive.

## 8. `read_vector_c81d` is removed in digital_rf 3 — DONE

`outlier_lpi.py`, `avg_range_doppler_spec.py`, `plot_raw_voltage.py`,
`tx_delay.py`

17 call sites replaced with
`read_vector_1d(...).astype("c8", casting="unsafe", copy=False)`, the form the
library's own deprecation notice recommends.

Checked before replacing rather than assumed: on this data `read_vector_1d`
already returns `complex64`, so the values are bit identical
(`n.array_equal` true, max difference 0.0) and the `astype` is a no-op. 20
reads took 0.011 s against 0.015 s for the old call, so nothing is lost.

Verified afterwards by running the voltage path with
`-W error::FutureWarning`, which now completes instead of raising.

`plasma_line_clicker.py` was listed here in error: it makes no `read_vector`
calls at all.

## 9. Transmit pulse timing tables are duplicated and experiment specific — DONE

`radar_timing.py` (new), `outlier_lpi.py:27`, `avg_range_doppler_spec.py:19`,
`tx_delay.py:38`, `fit_lp.py:137`

The `tmm` dictionary was defined three times and the copies had drifted:

| | coded sweepids 1-32 | mode 300 | mode 800 | `read_length` |
|---|---|---|---|---|
| `outlier_lpi.py` | yes | yes | no | no |
| `avg_range_doppler_spec.py` | no | yes | yes | yes |
| `tx_delay.py` | yes | yes | yes | no |

Checked field by field before merging rather than assumed: every value the
copies shared was identical, so the union is unambiguous. One table now lives
in `radar_timing.py` with `T_INJECTION`, `CODED_SWEEPIDS` and the measured
`TX_PULSE_LENGTH_US` of item 5, which `fit_lp.py` had its own copy of.

Merging is safe even though `avg_range_doppler_spec.py` uses membership of the
table as a filter (`if sid[key] not in tmm.keys()`, line 302): the pulse loop
already skips on `sid[key] != mode` at line 298, so the coded entries the table
gains are never reached.

Verified against the pre-change code rather than by inspection:

- all five modules import; the merged table matches the old values at all 34
  keys, field by field.
- `tx_delay.py` reproduces its recorded delays exactly: zenith-l 11.880 us,
  misa-l 12.340 us on eclipse2024.
- a 20 s LPI run on misa-l writes files bit identical to the same run on the
  pre-change checkout, all 20 datasets in each of the two files.
- a 20 s range-doppler run on misa-l mode 300 likewise, both files, 16
  datasets each.

The 16 `LinAlgError: Singular matrix` messages that run prints at the start of
eclipse2024 appear identically in both checkouts, so they are pre-existing and
not from this change. Worth a look on their own, but the code catches them and
moves on.

Not done: loading the table from the config or the metadata. Nothing in the
recording carries these numbers, so there is nothing to read them from; what is
fixed is that one experiment's values now sit in one place. Left on the list
below as a smaller item.

## 10. "halting" message that does not halt — DONE

`outlier_lpi.py:100`, `outlier_lpi.py:334`, `outlier_lpi.py:510`,
`avg_range_doppler_spec.py:44`, `avg_range_doppler_spec.py:303`

Three messages that described something other than what the code did:

1. The pulse loop printed "unknown pulse code %d encountered, halting", then
   `continue`d, then had an unreachable `exit(0)`. It now counts the skipped
   pulses per code and prints one line per integration period after the loop:
   "pulse code %d is not in the timing table, skipped %d pulses". Printing at
   the point of the skip would give one line per pulse, hundreds per period,
   which is why the old message mattered: a log full of "halting" from a run
   that was still going.

2. `estimate_dc`, in both files, printed "unknown pulse, ignoring" and then did
   not ignore — it fell through to `tmm[sid[key]]["last_echo"]` and would have
   raised `KeyError` on the very pulse it claimed to skip. The missing
   `continue` is added. Both copies are uncalled, found while merging the
   timing tables in item 9, so this was latent rather than live.

3. `avg_range_doppler_spec.py`'s own membership check had its message commented
   out. It is unreachable — the `sweepid != mode` test above it has already
   selected a mode that is in the table — so the comment now says that instead
   of leaving a dead `print`.

Verified: a 20 s LPI run on misa-l writes files bit identical to the previous
commit, and the same run with sweepid 5 deleted from the table exercises the
new path, reporting one line per integration period, "pulse code 5 is not in
the timing table, skipped 22 pulses" and then 24 for the second period, where
the old code would have printed 46 claims of halting.

## 11. Avoidable matrix work in the lag profile inversion — DONE

`outlier_lpi.py:683`

The comment asked for `AA=n.dot(AA,Sinv)` first, to save the `Sinv` dot
products, and said there was no time to test it. Done and tested.

The inversion now forms `B = Sinv A` once and uses it on both sides:

| | was | now |
|---|---|---|
| Fisher information | `(A^H Sinv) . (Sinv A)` | `B^H B` |
| right hand side | `(A^H Sinv) . mm` | `B^H mm` |

Algebraically identical: `Sinv` is real and diagonal, so `B^H B = A^H Sinv^2 A`
either way, and `mm` already carries one factor of `1/sigma`, so `B^H mm` is
still `A^H Sinv^2 m`. One sparse diagonal product per lag is saved.

Measured on a 20 s misa-l run, the two runs back to back on an otherwise idle
machine, 74 solves of a 121068 x 117 design matrix:

| | baseline | reordered |
|---|---|---|
| solve time | 5.46 s | 4.56 s |
| per solve | 0.0738 s | 0.0616 s |

16 % off the solve. Worth keeping in proportion: the same run spends 45 s
building ambiguity functions, so this is about 2 % of the work, and the comment
was right that it was never going to be dramatic.

Output is bit identical, all 20 datasets in both files, so the estimates and
their errors are unchanged, not merely close.

## 12. Range gates of a single size only — MEASURED, decision pending

`outlier_lpi.py:149`

"tbd: add range gates of different sizes". Range resolution is uniform across
the profile, though the useful resolution is not: at 30 us gates the zero-lag
SNR per gate per file is 1.3-1.4 below 400 km, 0.4 to 600 km and 0.1 above, so
two thirds of the gates carry almost nothing on their own.

Measured before deciding, because the case for changing the inversion rested on
an error argument that turned out to be backwards.

**Gate-to-gate correlation.** The inversion computes the whole `Sigma` and
stores only `n.diag(Sigma)`. Instrumented to keep all of it, 44 lags of one
period:

| separation | correlation |
|---|---|
| 1 gate (4.5 km) | -0.338 |
| 2 gates | +0.007 |
| 3-5 gates | -0.06 to -0.02 |

-0.338 in every band from 250 to 1100 km. Negative, and flat with range: the
signature of a deconvolution, not of SNR. A 480 us pulse cannot separate 4.5 km
gates, so the inversion trades power between neighbours.

**What that does to `fit_lpi`'s averaging.** `fit_lpi.py:805` combines gates as
`1/sum(1/var)`, assuming independence. Because the correlation is negative, the
reported variance of the average is too *large*, not too small. Applying the
code's own window, weights and variance rule to the archive and repeating the
item 2b scatter test on the averaged quantities:

| band | range_avg | reported/observed | reported sigma |
|---|---|---|---|
| misa-l 0-300 km | 0 -> 1 | 1.53 -> 1.78 | 1.75x -> 1.89x |
| misa-l 300-700 km | 0 -> 3 | 1.61 -> 3.22 | 1.80x -> 2.54x |
| misa-l 700-1500 km | 0 -> 5 | 1.86 -> 4.76 | 1.93x -> 3.09x |
| zenith-l 0-300 km | 0 -> 1 | 1.36 -> 1.65 | 1.65x -> 1.81x |
| zenith-l 300-700 km | 0 -> 3 | 1.34 -> 1.81 | 1.63x -> 1.90x |
| zenith-l 700-1500 km | 0 -> 5 | 1.70 -> 4.71 | 1.84x -> 3.07x |

Six of six move the same way: averaging compounds the overstatement.

**Consequences.** The hypothesis that the ~7% agreement reported in memo 6 was
two opposing errors cancelling is refuted; both effects push the same way. The
7% claim itself does not reconstruct from its own table and has been corrected
in the memo to 1.83x. Item 12 itself is not self-contained: changing gate widths
invalidates whatever calibration `var_scale` carries.

**Parameter level, measured directly.** Refitted an 8 h window on zenith-l
(fixed pointing) with and without range averaging, comparing reported errors
against the scatter between fits 400 s apart. First attempted on misa-l and
discarded: MISA alternates between two azimuths 89 deg apart between periods,
so consecutive fits see different volumes.

47 fits, 31 pairs 400 s apart:

| band | Te none / [1,3,5] | Ti | vi | ne |
|---|---|---|---|---|
| 250-400 km | 1.19 / 1.15 | 1.14 / 1.22 | 0.82 / 0.95 | 1.22 / 1.09 |
| 400-600 km | 1.08 / 1.13 | 1.06 / 1.32 | 0.75 / 1.20 | 1.25 / 1.26 |
| 600-800 km | 1.70 / 1.50 | 1.51 / 1.58 | 0.99 / 1.35 | 2.59 / 2.49 |
| 800-1100 km | 3.02 / 3.32 | 2.79 / 3.53 | 2.06 / 2.36 | 4.50 / 9.91 |

Below 600 km every entry is between 0.75 and 1.32, median 1.15: the reported
uncertainties are right to within about 30% there. Above it they are overstated,
badly so on the topside. Range averaging moves Te by at most 12% and ne below
800 km by at most 11%, Ti by up to 27% and vi by up to 60%, despite inflating
the ACF variance by 2-2.5x -- the nonlinear fit absorbs most of it.

Do: nothing to item 12 itself. The error argument for variable gates is gone,
the compute argument is worth ~2% of a run, and the cost is re-calibrating the
chain. Defer. What the measurement did turn up is worth its own item, below.

## 13. Hardcoded range gates in the plasma line clicker

`plasma_line_clicker.py:96`

"tbd fix hard coded range gates!" Affects the magic constant calibration
workflow.

## 14. Consider digital_metadata for the LPI output

`outlier_lpi.py:719`

"tbd: determine if this could be done better with digital_metadata". Worth a
decision, then either do it or delete the comment.

## 15. Stale comment

`outlier_lpi.py:210`

"this is experiment specific. need to read from configuration eventually" —
already done: `run_analysis.py` passes `max_range_delay_us` from the config.
Delete the comment.

## 16. Timing table is still a constant

`radar_timing.py`

Item 9 left one table instead of three, but it still states one experiment
configuration in source. A second experiment with different gates would need it
edited rather than configured.

Do: allow the config to supply or override the table, once there is a second
experiment to test it against. Reading it from the recording is not possible:
the metadata carries sweepids, not gate boundaries.

## 16b. Range averaging slices with a negative lower bound

`fit_lpi.py:805`

`avg_var[ri]=1/nansum(1/var_orig[(ri-ra):min(N,ri+ra)])` has no `max(0,...)` on
the lower bound, so for `ri < ra` the slice wraps and comes back empty, the sum
is zero and the averaged variance is `inf`. The line above it, which averages
the ACF, does guard with `max(0,...)`.

Latent only: the fit skips gates below 250 km and `ra` is at most 5, so the
affected gates are never fitted. Two characters to fix.

Note also that the window is `[ri-ra, ri+ra)`, which is `2*ra` gates, while
`range_avg_window_km` at line 957 reports `(2*range_avg+1)` gates. One of the
two is wrong.

## 12b. Topside uncertainties are overstated, and worse for ne

`fit_lpi.py`, all bands above 600 km

From the parameter-level test of item 12, on zenith-l with fixed pointing:
reported errors exceed the observed scatter by 1.5-1.7x at 600-800 km and
2.8-3.5x above 800 km for Te and Ti, and by 4.5-9.9x for ne above 800 km. Below
600 km everything agrees to within 30%, so this is a range-dependent effect,
not the overall scale.

Nothing measured so far explains a range dependence that strong. The candidates
worth checking: the r^2 range weighting, the ne_const scaling path (ne is much
the worst, which points here), and whether the fit is hitting parameter bounds
where SNR is low.

Do: investigate. This is where the remaining uncertainty error actually lives.

## 17. Singular normal matrix at the start of a recording

`outlier_lpi.py:691`

A 20 s LPI run at the very start of eclipse2024 raises
`LinAlgError: Singular matrix` sixteen times, preceded by
`invalid value encountered in divide` on `mm_g/sigma_lp_est`. The exception is
caught and the period is skipped, so nothing is corrupted, but the log is
alarming and integration periods are silently lost. Present before the item 9
merge as well, so it is pre-existing and not a regression.

Do: find which range gates give `sigma_lp_est` of zero or NaN, and either
exclude them from the design matrix or say plainly in the log that the period
was dropped and why. A rank check before the inverse would be cheaper than
catching the exception.

---

## 18. `outlier_lpi.py` biases T_sys low, by a signal dependent factor

`outlier_lpi.py:414-417`, `outlier_lpi.py:526-528`

The lag profile inversion estimates the DC offset as the mean of the 500
microsecond background window and subtracts it from both the background and the
noise injection window:

```python
z_dc_noise = n.mean(z_noise[(last_echo-500):last_echo])
bg_samples.append(          n.mean(n.abs(z_noise[(last_echo-500):last_echo] - z_dc_noise)**2.0) )
bg_plus_inj_samples.append( n.mean(n.abs(z_noise[(noise0):noise1]          - z_dc_noise)**2.0) )
```

That mean is not a DC offset. It is an estimate of one from `N` independent
samples, where `N` is the *filtered* bandwidth times the window length, not the
500 samples the window holds. `z_noise` has already been through
`lpf.lpf()` at `1.2 * pass_band`, so for the 18 kHz passband this experiment
used,

    N = 2 * 1.2 * 18000 Hz * 500 us = 21.6

Subtracting it does two things, both in the same direction:

1. it removes `1/N` of the background window's own variance, and
2. it adds an uncorrelated `|mean|^2 = P_bg/N` to the *injection* window, which
   the background window's mean has nothing to do with.

So the estimator returns

    T_est = T_inj * P_bg (1 - 1/N) / [ (P_inj - P_bg) + 2 P_bg / N ]

instead of `T_inj * P_bg / (P_inj - P_bg)`. Writing `y = T_inj / T_true` for the
true ratio `(P_inj - P_bg)/P_bg`, the inverse is

    T_true = T_inj / [ T_inj (1 - 1/N) / T_est  -  2/N ]

The bias is **not** a constant offset. It grows as the band narrows and as
T_sys itself grows, because `P_inj - P_bg` shrinks relative to `P_bg`.

### Measured

2024-04-08, `zenith-l`, 749 coded pulses at the Cygnus A transit and the same
processing at baseline. The same pulses, reprocessed at three bandwidths with
and without the DC subtraction:

| passband | DC removed | DC kept |
|---|---|---|
| full band (1 MHz) | 1573 K | 1584 K |
| 50 kHz | 1508 K | 1612 K |
| 18 kHz | 1365 K | 1607 K |

With the DC kept, all three bandwidths agree to 2 %, which is what a flat
spectrum source requires. The spread appears only when the DC is removed.

The formula predicts every value with no free parameter:

| case | predicted | measured |
|---|---|---|
| 18 kHz, peak | 1355 K | 1365 K |
| 50 kHz, peak | 1504 K | 1508 K |
| DC kept, peak | 1600 K | 1607 K |
| 18 kHz, baseline | 151 K | 159 K |

so the bias is 6 % at a 160 K baseline and 15 % at the 1600 K peak.

`avg_range_doppler_spec.py` is not affected: it subtracts a single global
constant estimated over many pulses (`estimate_dc`), which is a DC offset, and
does not re-estimate one per window.

### Why it matters beyond T_sys

`alpha` is stored alongside `T_sys` and is what scales power to temperature for
the electron density calibration (`fit_lpi.py`, `fit_lp.py`, and the magic
constant workflow). `alpha` carries the same defect, so `ne` inherits a signal
dependent scale error wherever the LPI product is used.

### Confirmation

Correcting the stored `lpi_30` values with the formula above and fitting them
against the 440 MHz sky model gives, for the zenith antenna,
`T_0 = 155.3 +- 1.7 K`, `eta = 0.654 +- 0.060`, against
`T_0 = 156.2 +- 1.6 K`, `eta = 0.661 +- 0.055` from
`range_doppler_300_outlier`. Two estimators with different pulse sets, gates
and passbands, agreeing to one per cent once the bias is undone.

### Not done here

The fix is to estimate the DC the way `avg_range_doppler_spec.py` does, over
many pulses, and never per window. That changes the inversion's output and the
`ne` calibration needs revalidating against it, so it does not belong on this
branch. `compare_tsys_sky.py:lpi_bias_factor` carries the derivation and the
correction in code in the meantime.

---

## 19. MISA pointing cannot be resolved per integration period

`tsys_harvest.py:attach_pointing`

MISA alternates its azimuth between integration periods. Interpolating
`misa_azimuth` from `antenna_control_metadata` at the centre of a ten second
period gives a value that swings by tens of kelvin of modelled sky temperature
from one period to the next, while the measured `T_sys` stays flat across the
same periods. The pointing assignment, not the measurement, is what is wrong.

Consequences, all visible in `compare_tsys_sky.py` output:

- error in the modelled sky temperature is error in the regressor, which drags
  a least squares slope towards zero;
- MISA's two estimators give `eta` 0.48 and 0.63, where the zenith antenna's
  agree to 1 %;
- `range_doppler_300_outlier` additionally spans elevation 37 to 45 degrees,
  mid slew, where `lpi_30` is parked at 45, and the two differ by 28 K in `T_0`.

Resolving it needs the pointing worked out per pulse from the raw metadata
timestamps and averaged over the pulses actually used, rather than interpolated
once per integration period. Not attempted. MISA is reported by
`compare_tsys_sky.py` but flagged as not a measurement.

---

## 20. The zenith antenna is 1.84 degrees off the geodetic zenith

`tsys_harvest.py:ZENITH_EL_DEG`

Not a defect in this repository, but it invalidates the obvious assumption and
nothing here recorded it. The hard target work measured the zenith beam at
elevation 88.16, azimuth 172.9 from 251 satellite matches, Rayleigh Z = 367,
with `misa-l` and its recorded encoder pointing as a control showing no offset.

The 440 MHz data confirm it independently and by a completely different route.
That tilt puts the boresight at declination

    42.619 + 1.84 * cos(172.9) = 40.79 degrees

and Cygnus A sits at 40.734, so the two agree to 0.06 degrees, a twelfth of a
beamwidth. Cyg A accordingly transits the beam centre once per sidereal day and
drives `T_sys` from 170 K to 1600 K. Assuming the geodetic zenith instead puts
the source 1.89 degrees off boresight, 2.7 beamwidths out, where a generous
sidelobe envelope leaves a few kelvin rather than 1400. The two transits in the
recording are 3 m 55 s apart against the 3 m 56 s sidereal day deficit, so it is
certainly celestial.

Whether the displacement is mechanical, a survey error or a convention in the
metadata is still unresolved.

---

## Documented approximations, not defects

These are deliberate simplifications that are worth knowing about when reading
results, listed so they are not mistaken for oversights:

- `fit_lpi.py:90`, `fit_lpi.py:103` — Debye length effects ignored.
- `il_interp.py:49` — monostatic geometry only.

## Hardware, out of software scope

`avg_range_doppler_spec.py:379`

The transmit sample should be interleaved into the echo channel with an analog
switch, so that `tx-h` and the echo share one receiver chain and the channel
delay correction is not needed at all. Until then the delay is measured from
the leakthrough (`tx_delay.py`), which fixes the relative alignment but leaves
a small systematic set by the coupling path length.

Related and still open: confirm the measured channel delay is stable across a
long recording. It was checked at the start of eclipse2024 only; the pipeline
measures once per run and applies that value to every integration period.
