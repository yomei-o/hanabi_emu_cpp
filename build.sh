#!/bin/bash
# Build one Hanabi OS sim to wasmdist/<name>/  (usage: ./build.sh <name>)
# 既定は hanabi_os。ABI / EXPORT_NAME(createSim) / 出力名(sim.js) は universe_cpp と共通。
#
# 光り方(トーンマップ)は**実行時の切替**になったのでビルドを分ける必要はない。
# ページの「✨ 光り方」ボタン、または sim_set(15, 1/0) で切り替える。
set -e; cd "$(dirname "$0")"
EMSDK="${EMSDK:-$HOME/emsdk}"
export EM_CONFIG="$EMSDK/.emscripten"
export PATH="$EMSDK/upstream/emscripten:$EMSDK/upstream/bin:$PATH"
NAME="${1:-hanabi_os}"
OUT="sim.js"; EXPNAME="createSim"
mkdir -p "wasmdist/$NAME"
em++ -O3 -std=c++17 -I. "src/$NAME.cpp" \
  -sMODULARIZE=1 -sEXPORT_NAME=$EXPNAME -sENVIRONMENT=web -sALLOW_MEMORY_GROWTH=1 \
  -sEXPORTED_FUNCTIONS=_sim_init,_sim_w,_sim_h,_sim_reset,_sim_step,_sim_render,_sim_click,_sim_set,_sim_action,_sim_get,_malloc,_free \
  -sEXPORTED_RUNTIME_METHODS=cwrap,HEAPU8 \
  -o "wasmdist/$NAME/$OUT"
echo "built wasmdist/$NAME/$OUT (+.wasm)  [$EXPNAME]"
