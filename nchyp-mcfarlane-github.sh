#!/usr/bin/env bash

set -o errexit
set -o pipefail
set -o nounset

SCRIPTPATH="$(
	cd -- "$(dirname "$0")" >/dev/null 2>&1
	pwd -P
)"

PROG_NAME=$0

INPUT_FILENAME="Genotypes.txt"
CHOICES_FILENAME="Choices.txt"
OUTPUT_FILENAME="Output.txt"

usage() {
	cat <<EOF
usage: ${PROG_NAME} [-h] [--input INPUT] [--choices CHOICES] [--output OUTPUT]

positional arguments:
  -i, --input           Input filename [default: ${INPUT_FILENAME}]
  -c, --choices         Choices filename [default: ${CHOICES_FILENAME}]
  -o, --output          Output filename [default: ${OUTPUT_FILENAME}]

optional arguments:
  -h, --help            show this help message and exit
EOF
	exit 1
}

VALID_ARGS=$(getopt -o hi:c:o: --long help,input:,choices:,output: -- "$@")
if [[ $? -ne 0 ]]; then
	exit 1
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
	case "$1" in
	-h | --help)
		usage
		shift
		;;
	-i | --input)
		INPUT_FILENAME=$2
		shift 2
		;;
	-c | --choices)
		CHOICES_FILENAME=$2
		shift 2
		;;
	-o | --output)
		OUTPUT_FILENAME=$2
		shift 2
		;;
	--)
		shift
		break
		;;
	esac
done

if [[ ! -f "${INPUT_FILENAME}" ]]; then
	echo "[!] Input file doesn't exist:"
	realpath "${INPUT_FILENAME}"
	exit 1
fi
if [[ ! -f "${CHOICES_FILENAME}" ]]; then
	echo "[!] Choices file doesn't exist:"
	realpath "${CHOICES_FILENAME}"
	exit 1
fi
if [[ -e "${OUTPUT_FILENAME}" ]]; then
	echo "[!] Output file already exists:"
	realpath "${OUTPUT_FILENAME}"
	exit 1
fi

WORK_DIR=$(mktemp -d)

cp "${SCRIPTPATH}"/*.m "${WORK_DIR}/"
cp "${INPUT_FILENAME}" "${WORK_DIR}/Genotypes.txt"
cp "${CHOICES_FILENAME}" "${WORK_DIR}/Choices.txt"

pushd "${WORK_DIR}" >/dev/null
octave -q "NcHyper240229Octave.m"
popd >/dev/null

cp "${WORK_DIR}/Output.txt" "${OUTPUT_FILENAME}"

# Cleanup
rm -rf "${WORK_DIR}"
