#!/bin/bash
# Build one Hanabi OS sim to wasmdist/<name>/  (usage: ./build.sh <name> [legacy])
# 既定は hanabi_os。ABI / EXPORT_NAME(createSim) / 出力名(sim.js) は universe_cpp と共通。
#
# 2番目に legacy を渡すと、**昔のトーンマップ**(成分ごとに r/(1+r)。色が抜ける)で
# もう1本ビルドする。出力名と EXPORT_NAME を分けるので、比較ページ(compare.html)から
# 新旧2つの wasm を同時に読み込んで並べられる。
#   ./build.sh hanabi_os          -> sim.js        (createSim)       新しい光
#   ./build.sh hanabi_os legacy   -> sim_legacy.js (createSimLegacy) 昔の光
set -e; cd "$(dirname "$0")"
EMSDK="${EMSDK:-$HOME/emsdk}"
export EM_CONFIG="$EMSDK/.emscripten"
export PATH="$EMSDK/upstream/emscripten:$EMSDK/upstream/bin:$PATH"
NAME="${1:-hanabi_os}"
if [ "$2" = "legacy" ]; then
  DEFS="-DLIGHT_LEGACY"; OUT="sim_legacy.js"; EXPNAME="createSimLegacy"
else
  DEFS="";               OUT="sim.js";        EXPNAME="createSim"
fi
mkdir -p "wasmdist/$NAME"
em++ -O3 -std=c++17 -I. $DEFS "src/$NAME.cpp" \
  -sMODULARIZE=1 -sEXPORT_NAME=$EXPNAME -sENVIRONMENT=web -sALLOW_MEMORY_GROWTH=1 \
  -sEXPORTED_FUNCTIONS=_sim_init,_sim_w,_sim_h,_sim_reset,_sim_step,_sim_render,_sim_click,_sim_set,_sim_action,_sim_get,_malloc,_free \
  -sEXPORTED_RUNTIME_METHODS=cwrap,HEAPU8 \
  -o "wasmdist/$NAME/$OUT"
echo "built wasmdist/$NAME/$OUT (+.wasm)  [$EXPNAME ${DEFS:-新しい光}]"
