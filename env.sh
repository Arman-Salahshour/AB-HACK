# Source me:  source env.sh
# Image tags and the docker helper used by every script in scripts/.

AB_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
export AB_ROOT

export Q=quay.io/apheris/oslo-hackathon

# --- design models (RFD4-Proteina, one image per weight) ---------------------
export RFD4_ANTIBODY="$Q:apheris-rfd4-proteina-v0.3.0-antibody-2026-10-30"
export RFD4_HER2="$Q:apheris-rfd4-proteina-v0.3.0-her2-2026-11-07"
export RFD4_HER2_NOVEL="$Q:apheris-rfd4-proteina-v0.3.0-her2-novel-2026-11-09"

# --- folding models ----------------------------------------------------------
# name                           image                                    --model        template
export FOLD_BOLTZ2="$Q:apheris-boltz-v1.3.0-boltz2"                   # boltz2            yes
export FOLD_BOLTZ1="$Q:apheris-boltz-v1.3.0-boltz1"                   # boltz1            no
export FOLD_OF3OB="$Q:apheris-openfold3-v1.3.0-of3ob"                 # openfold3-of3ob   yes
export FOLD_OF3P2="$Q:apheris-openfold3-v1.3.0-of3p2"                 # openfold3-of3p2   yes
export FOLD_PROTENIX_V2="$Q:apheris-protenix-v1.3.0-v2"               # protenix-v2       yes
export FOLD_PROTENIX_V1="$Q:apheris-protenix-v1.3.0-v1-default"       # protenix-v1       yes
export FOLD_MOL1="$Q:apheris-mol1-v1.3.0-base"                        # mol1-base         yes
export FOLD_ESMFOLD2="$Q:apheris-esmfold2-v1.4.0-release"             # esmfold2-release  no

# Map a short weight name to its image: $(rfd4_image her2-novel)
rfd4_image() {
  case "${1:-her2-novel}" in
    antibody|antibody-2026-10-30)        echo "$RFD4_ANTIBODY" ;;
    her2|her2-2026-11-07)                echo "$RFD4_HER2" ;;
    her2-novel|her2-novel-2026-11-09)    echo "$RFD4_HER2_NOVEL" ;;
    *) echo "unknown weight: $1" >&2; return 1 ;;
  esac
}

# Map a short folding name to "image|--model value"
fold_image() {
  case "${1:-boltz2}" in
    boltz2)       echo "$FOLD_BOLTZ2|boltz2" ;;
    boltz1)       echo "$FOLD_BOLTZ1|boltz1" ;;
    of3ob)        echo "$FOLD_OF3OB|openfold3-of3ob" ;;
    of3p2)        echo "$FOLD_OF3P2|openfold3-of3p2" ;;
    protenix-v2)  echo "$FOLD_PROTENIX_V2|protenix-v2" ;;
    protenix-v1)  echo "$FOLD_PROTENIX_V1|protenix-v1" ;;
    mol1)         echo "$FOLD_MOL1|mol1-base" ;;
    esmfold2)     echo "$FOLD_ESMFOLD2|esmfold2-release" ;;
    *) echo "unknown folding model: $1" >&2; return 1 ;;
  esac
}

# Models without template support (leave out --template-paths).
fold_has_template() {
  case "${1:-boltz2}" in boltz1|esmfold2) return 1 ;; *) return 0 ;; esac
}

# Run a Foundry image with the project root mounted at /work.
# Paths you pass are relative to the project root.
abrun() {
  docker run --rm --gpus all --shm-size 16g \
    --user "$(id -u):$(id -g)" -e HOME=/tmp \
    -v /etc/passwd:/etc/passwd:ro \
    -v "$AB_ROOT:/work" -w /work "$@"
}

export -f rfd4_image fold_image fold_has_template abrun 2>/dev/null || true

# --- External APIs (keys live in .env.private, gitignored) -------------------
# AMASS  : AMASS_API_KEY / AMASS_BASE_URL  (https://api.amass.tech)
# JapanFold: JAPANFOLD_API_KEY / JAPANFOLD_BASE_URL (https://api.japanfold.aiand.com)
if [ -f "$AB_ROOT/.env.private" ]; then
  # shellcheck source=/dev/null
  set -a; . "$AB_ROOT/.env.private"; set +a
fi
if [ -f "$HOME/.claude/skills/amass/SKILL.md" ] || [ -f "$HOME/.claude/skills/japanfold/SKILL.md" ]; then
  : # skills present; keys should be exported already via .env.private
fi
