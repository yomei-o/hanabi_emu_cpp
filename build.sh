#!/bin/bash
# Build one Hanabi OS sim to wasmdist/<name>/  (usage: ./build.sh <name> [vivid])
# 既定は hanabi_os。ABI / EXPORT_NAME(createSim) / 出力名(sim.js) は universe_cpp と共通。
#
# 2番目に vivid を渡すと、**色を保つトーンマップ**(明るさだけ圧縮)でもう1本ビルドする。
# 出力名と EXPORT_NAME を分けるので、比較ページ(compare.html)から2つの wasm を
# 同時に読み込んで並べられる。既定は成分ごとの r/(1+r)(落ち着いた見え方)。
#   ./build.sh hanabi_os         -> sim.js       (createSim)      既定
#   ./build.sh hanabi_os vivid   -> sim_vivid.js (createSimVivid) 色を保つほう
set -e; cd "$(dirname "$0")"
EMSDK="${EMSDK:-$HOME/emsdk}"
export EM_CONFIG="$EMSDK/.emscripten"
export PATH="$EMSDK/upstream/emscripten:$EMSDK/upstream/bin:$PATH"
NAME="${1:-hanabi_os}"
if [ "$2" = "vivid" ]; then
  DEFS="-DLIGHT_VIVID"; OUT="sim_vivid.js"; EXPNAME="createSimVivid"
else
  DEFS="";              OUT="sim.js";       EXPNAME="createSim"
fi
mkdir -p "wasmdist/$NAME"
em++ -O3 -std=c++17 -I. $DEFS "src/$NAME.cpp" \
  -sMODULARIZE=1 -sEXPORT_NAME=$EXPNAME -sENVIRONMENT=web -sALLOW_MEMORY_GROWTH=1 \
  -sEXPORTED_FUNCTIONS=_sim_init,_sim_w,_sim_h,_sim_reset,_sim_step,_sim_render,_sim_click,_sim_set,_sim_action,_sim_get,_malloc,_free \
  -sEXPORTED_RUNTIME_METHODS=cwrap,HEAPU8 \
  -o "wasmdist/$NAME/$OUT"
echo "built wasmdist/$NAME/$OUT (+.wasm)  [$EXPNAME ${DEFS:-既定}]"
