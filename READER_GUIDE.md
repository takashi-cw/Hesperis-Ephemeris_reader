# Hesperis Ephemeris Reader — Reader Guide

This document is a navigation index for humans and AI readers alike.
It describes what each document covers, which sources are authoritative,
and what should not be treated as specification.

---

## What this project is

Hesperis Ephemeris Reader is a collection of BSP (SPK Type 2 & Type 3) readers
for NASA JPL planetary ephemeris files (.bsp), independently implemented
in multiple languages (JavaScript, Python, Swift, Rust, Julia).

It reads JPL DE-series kernels and returns Chebyshev-interpolated
barycentric position/velocity vectors. No external astronomy library is used.

---

## Primary documents

| Document | Role |
|----------|------|
| `README.md` | Project overview, purpose, implementation list, and basic usage |
| `READER_GUIDE.md` | This file — navigation index and authoritative source guide |

### Per-language entry points

| Language | README | Source |
|----------|--------|--------|
| JavaScript | `js/README.md` | `js/src/` |
| Python | `py/` | `py/` |
| Swift | `swift/` | `swift/` |
| Rust | `rust/` | `rust/` |
| Julia | `julia/` | `julia/` |

---

## Authoritative sources

- **For implementation behavior**: the current source code in each language's `src/` directory.
- **For project scope and purpose**: `README.md`.
- **For license questions**: `LICENSE` (in each language subdirectory where present).

---

## Do not treat as authoritative

- Any intermediate output files
- Scratch notes or temporary experiments
- Comments referencing external libraries — these appear only as comparison references, not as implementation sources

---

## Important caution

Hesperis returns raw barycentric position vectors (km, ICRF frame).
It does not compute:
- Zodiac sign positions
- Astrological house placements
- Apparent (topocentric or geocentric ecliptic) positions
- Any calendar or timezone conversions

Do not describe BSP output and astrological chart output as equivalent.
Coordinate transforms and reduction stages are the responsibility of the calling application.

---

## License

MIT. No Swiss Ephemeris lineage. Pure independent implementations.