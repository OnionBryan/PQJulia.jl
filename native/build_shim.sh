#!/usr/bin/env bash
# Build libforgeddec.dylib — the C-ABI bridge to Forged-lab's audited native DEC
# (geom build_hodge2d/3d circumcentric Hodge stars + geom_su2_wilson holonomy).
# Links the prebuilt libgeom.a + libm5_math.a from the Forged-lab native tree.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
FORGED="${FORGED_LAB:-/Users/bryan/PycharmProjects/Forged-lab}"

GEOM_A="$FORGED/native/build/src/geom/libgeom.a"
MATH_A="$FORGED/native/build/libm5_math.a"
for a in "$GEOM_A" "$MATH_A"; do
  [ -f "$a" ] || { echo "missing archive: $a (build the Forged-lab native tree first)"; exit 1; }
done

# Recompile the FEEC Galerkin source (now with the k=2 Whitney 2-form mass) and
# link it AHEAD of libm5_math.a so this object wins over the archive's k≤1 copy.
# (Avoids a full native rebuild; the native source itself is already updated.)
GALERKIN_SRC="$FORGED/native/src/math/src/core/hodge_galerkin.cpp"
clang++ -std=c++17 -O2 -fPIC -c "$GALERKIN_SRC" \
  -I"$FORGED/native/src/math/include" -I"$FORGED/native/src" \
  -o "$HERE/hodge_galerkin.o"

clang++ -std=c++17 -O2 -fPIC -dynamiclib \
  -I"$FORGED/native/src/geom/include" \
  -I"$FORGED/native/src/math/include" \
  -I"$FORGED/native/src" \
  "$HERE/forged_dec_shim.cpp" \
  "$HERE/hodge_galerkin.o" \
  "$GEOM_A" "$MATH_A" \
  -framework Accelerate -framework Metal -framework Foundation \
  -o "$HERE/libforgeddec.dylib"

echo "built $HERE/libforgeddec.dylib"
nm -gU "$HERE/libforgeddec.dylib" | grep forged

# Guard the k=2 Whitney 2-form mass added to native galerkin_hodge_star.
clang++ -std=c++17 -O2 \
  -I"$FORGED/native/src/math/include" -I"$FORGED/native/src" \
  "$HERE/test_galerkin_k2.cpp" "$GALERKIN_SRC" "$MATH_A" \
  -framework Accelerate -o "$HERE/test_galerkin_k2" && "$HERE/test_galerkin_k2"
