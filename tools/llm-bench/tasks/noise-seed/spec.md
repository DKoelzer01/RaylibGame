# The 1D and 2D noise fractals take the seed

**Problem:** in `src/engine/SimplexNoise.cpp`, `SimplexNoise::fractal(size_t octaves, float x, float y, float z)` offsets
its coordinates by the noise's seed offsets (`mSeedX`, `mSeedY`, `mSeedZ`), but the 1D overload
`fractal(size_t octaves, float x)` and the 2D overload `fractal(size_t octaves, float x, float y)` don't: two noises with
different seeds give the same 1D and 2D fields.

**Change:** in the 1D and 2D overloads, offset `x` by `mSeedX` (and `y` by `mSeedY`) before scaling by the octave's
frequency, the same way the 3D overload does. Nothing else changes.

**Test:** `SimplexNoise.LowerDimensionFractalsAreSeeded` (already in the file).
