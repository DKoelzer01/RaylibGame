# Weather fields follow the world seed

**Problem:** in `src/space/weather.cpp`, the weather noise fields (`Fields fields`, filled by `Init(unsigned int seed)` in
the anonymous namespace) are built once per process: `Init` returns early when `fields.ready` is set. Loading a second
world with another seed in the same session keeps the first world's weather.

**Change:** `Init` rebuilds the fields when it is called with a different seed than the one they were built with, and
returns early only for the same seed. Remember the seed in `Fields` (an `unsigned int seed` member).

**Test:** `Weather.FieldsFollowTheWorldSeed` (already in the file): two seeds give different cloud values, the same seed
gives the same value again.
