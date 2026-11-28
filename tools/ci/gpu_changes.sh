#!/bin/bash
# (#1392) Whether the self-hosted GPU steps are needed for this run: --test-terrain-gpu (the GPU terrain bit-exact against the
# CPU field) and the fluid CPU = GPU parity. They are slow, and the runners are the bottleneck; they only change verdict when
# shaders, terrain generation, fluid or render code change. Writes gpu=true|false to $GITHUB_OUTPUT (and says why).
# Everything runs for a schedule, a manual run, a branch's first push (no base to diff), and when the diff can't be read.
# Env: EVENT (github.event_name), BEFORE (github.event.before), AFTER (github.sha).
set -u
out="${GITHUB_OUTPUT:-/dev/stdout}"
say() { echo "gpu=$1" >> "$out"; echo "GPU steps: $1 ($2)"; [ -n "${GITHUB_STEP_SUMMARY:-}" ] && echo "GPU steps: $1 ($2)" >> "$GITHUB_STEP_SUMMARY"; exit 0; }
[ "${EVENT:-}" = push ] || say true "a ${EVENT:-manual} run runs everything"
case "${BEFORE:-}" in ""|0000000000000000000000000000000000000000) say true "a new branch: no base to compare";; esac
# (#2009 incident) --depth only when the clone is already shallow: a depth fetch makes a full clone shallow
depth=""; [ "$(git rev-parse --is-shallow-repository 2>/dev/null)" = true ] && depth="--depth=1"
git cat-file -e "$BEFORE^{commit}" 2>/dev/null || git fetch --quiet $depth origin "$BEFORE" 2>/dev/null || say true "the base $BEFORE isn't in the checkout (the checkout fetches 50 commits)"
files=$(git diff --name-only "$BEFORE" "${AFTER:-HEAD}" 2>/dev/null) || say true "the diff couldn't be read"
# What the GPU steps exercise: every shader; the terrain field and its GPU twin (density, noise, erosion, rivers, strata,
# caves live in voxel), meshing, the world-generation data and bodies; the fluid code; the render targets; and the build.
gpuRe='^resources/.*\.(glsl|fs|vs|comp)$'
gpuRe+='|^src/engine/(voxel|SimplexNoise|erosion|canyons|rivers|strata|terraingpu|terrainarena|terraindeltas|terrainlayers|cubemarch|cubemarchconsts|gpucompute|gpumem|glcaps|worldgenparams|worldprofile|bodies|world|fluid|fluidkernel|waterfield|watersurface|scenetarget|skyphys|farsphere)\.(h|cpp)$'
gpuRe+='|^resources/data/(bodies|biomes|worldgen|strata|fluids)[^/]*\.txt$'
gpuRe+='|^(CMakeLists\.txt|tools/build_linux\.sh|\.github/workflows/selfhosted\.yml|tools/ci/gpu_changes\.sh)$'
hit=$(printf '%s\n' "$files" | grep -E "$gpuRe" | head -3 | tr '\n' ' ')
[ -n "$hit" ] && say true "changed: $hit"
say false "no shader, terrain, fluid, render or build file changed in ${BEFORE:0:7}..${AFTER:0:7}"
