# Spin-orbit symmetry convention in dft_tools `dmftproj` (issue #148)

A from-scratch derivation of the spinor symmetry convention used by the
`dmftproj` projector tool in [TRIQS/dft_tools](https://github.com/TRIQS/dft_tools),
with a no-Wien2k reproducer, written to investigate
[issue #148](https://github.com/TRIQS/dft_tools/issues/148) (wrong symmetries
for spin-orbit + spin-polarized calculations).

Companion regression test: [TRIQS/dft_tools#292](https://github.com/TRIQS/dft_tools/pull/292).

The Mathematica scripts under `derivation/` reproduce every claim below; run them
with `wolfram -script derivation/<name>.wl`.

## What the issue proposes

johanjoensson (2020) reports that for spin-polarized + spin-orbit Wien2k runs on
the full BZ, symmetrizing the already-symmetric local Hamiltonian/density matrix
shifts its eigenvalues. Proposed fix: flip the sign of `factor` at four sites in
`outputqmc.f`:

```
-               ephase=EXP(CMPLX(0.d0,factor))
+               ephase=EXP(CMPLX(0.d0,-factor))
```

Reporter's own caveats, verbatim: "we might just have been lucky"; "This fix does
nothing for the case when the user supplies their own (non-spin-diagonal) basis
transformation (using the 'fromfile' option)." The change was confirmed
empirically by two contributors on the thread but never derived or merged.

## Convention (pinned)

z-y-z Euler angles. The code's own reference spinor matrix (`setsym.f:863-869`,
attributed to Wien2k `SRC_lapwdm/sym.f`):

```
spmt = [[ e^{+i(a+g)/2} cos(b/2),  e^{-i(a-g)/2} sin(b/2) ],
        [-e^{+i(a-g)/2} sin(b/2),  e^{-i(a+g)/2} cos(b/2) ]]   (det 1, unitary)
```

`setsym.f:147-165`: for SP+SO, `time_inv=.TRUE.` iff the xy-block determinant
(`= cos b`) is negative (the spin-flip / b=pi ops), with `phase = g-a`; else
`time_inv=.FALSE.`, `phase = a+g`. The Python consumer (`symmetry.py:150-162`)
applies `M_sym = (1/N) sum_g mat_g . (M or conj M) . mat_g^dag`, conjugating M for
`time_inv` ops. So the antiunitary spin-flip K is carried by `time_inv`; the
written `mat` is the unitary remainder.

## Symbolic check (`part1b_clean.wl`)

`spinrotmat` (setsym) reproduces `spmt` exactly for b=0 and b=pi, including the
`-CONJG` at `setsym.f:833`. `outputqmc`'s non-mixing block matches `spmt` for b=0
but not for b=pi (diagonal vs antidiagonal); this is expected, since b=pi carries
time reversal factored into `time_inv`. The proposed flip turns the (correct) b=0
up/up block into its complex conjugate.

## Idempotency reproducer

A valid (anti)unitary representation makes the symmetrizer an idempotent
projector, so an already-symmetric M keeps its eigenvalues. A sign or structure
error breaks idempotency, which is exactly the reported failure. This needs no
Wien2k data.

- `mre2.wl` (naive orbital model, magnetic group D4 = 4 z-rotations + 4 in-plane
  C2): l=0 is idempotent for both signs (the flip is a no-op); l=1, l=2 are not
  idempotent for either sign.
- `diag.wl`: the orbital-only symmetrizer is idempotent and the Wigner rep closes
  exactly for l=1,2, so the orbital model is faithful. The l>=1 failure in
  `mre2.wl` localizes to the spin-orbital coupling of the time-reversal ops.
- `faithful_mre.wl`: replicates the real construction exactly: `dmat` + `d_matrix`
  (setsym.f:694-776), the orbital `timeinv_op` (`T_{m,-m}=(-1)^m`,
  `mat->T conj(mat)`, applied to rotl for time_inv ops per setsym.f:325), and
  outputqmc's non-mixing exported `mat`, then the symmetry.py symmetrizer. For D4,
  single correlated atom, l=1 and l=2:
  - CURRENT `exp(+i phase/2)`: idempotency 1e-16, eigenvalue shift 1e-16.
  - FLIP `exp(-i phase/2)`: idempotency 1e-16, eigenvalue shift 1e-15.

Both signs give an idempotent (correct) symmetrizer. Idempotency of
`(1/N) sum_g g.(M or conjM).g^dag` is invariant under a fixed unitary basis change
`reptrans`, so no choice of basis (cubic harmonics, `fromfile`) changes this for a
single atom.

## Conclusion

The single-atom spinor matrices are a valid representation for both sign
conventions, so the proposed four-site flip is a no-op on every case reproducible
in isolation.

This was then checked against full Wien2k calculations (Wien2k 14.2, full BZ,
SOC + spin-polarized, the reporter's no-op symmetrization test):

| case | symmetry | corr. atoms | time_inv ops | eigenvalue drift |
|------|----------|-------------|--------------|------------------|
| Sr2MgOsO6 | tetragonal | 1 | 0 | 4.5e-7 |
| CaOs2 (fluorite) | cubic | 2 (equivalent) | 8 | 1.3e-7 |

CaOs2 has every condition the bug should need: cubic symmetry, two
symmetry-equivalent correlated atoms (the inter-atomic mapping
`jorb = R[isym](iorb)`), eight time-reversal operations, and the non-mixing
export path. It still symmetrizes correctly. The decisive test: rebuilding
dmftproj with the proposed `factor -> -factor` flip and re-running on the same
data **changes** `case.symqmc` (the matrices differ) but leaves the symmetrization
result **unchanged** (drift still 1.3e-7). The two signs are equally valid
representations differing by a phase that cancels in `D M D^dag`.

So on real cubic multi-equivalent-atom SOC + spin-polarized data, the current
dmftproj symmetrization is correct and the proposed flip is physically inert. The
reported failure does not reproduce; it may be specific to Wien2k 18.2, or to the
`fromfile` (spin-mixing) basis path that the reporter said the flip does not fix
and which is not tested here.

One structural observation for maintainers: the internal symmetrizer
`symmetrize_mat.f:99-110,231-242` puts the spin phase on the off-diagonal
(up/dn, dn/up) blocks with `ephase=exp(+-i*phase)` and `ephase=1` on the diagonal
blocks, whereas `outputqmc.f:758-782` puts `ephase=exp(+i*phase/2)` on the diagonal
blocks. These are different representations of the same operation; the relation
between the two is worth a look.

## License

GPL-3.0-or-later, matching upstream TRIQS.
