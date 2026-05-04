#!/bin/bash
BASE="$(cd "$(dirname "$0")" && pwd)"
NP="${NP:-20}"
TESTDIR="${1:-$BASE/tests/al_qha_el_litfix_new}"

for v in "$TESTDIR"/v*/; do
    vname=$(basename "$v")
    echo "===== $vname ====="
    cd "$v"
    scf_input=$(find . -maxdepth 1 -type f -name '*.in' | head -1)
    if [ -z "$scf_input" ]; then
        echo "  SCF INPUT NOT FOUND"
        cd "$BASE"
        continue
    fi

    if grep -q 'JOB DONE' qe.out 2>/dev/null; then
        echo "  SCF: done"
    else
        echo "  SCF: running..."
        mpirun -np "$NP" pw.x -input "$scf_input" > qe.out 2>&1
        grep -q 'JOB DONE' qe.out || { echo "  SCF FAILED"; cd "$BASE"; continue; }
        echo "  SCF: done"
    fi

    for s in elastic/*/*/; do
        [ -d "$s" ] || continue
        infile=$(ls "$s"*.in 2>/dev/null | head -1); [ -z "$infile" ] && continue
        outfile="${infile%.in}.out"
        grep -q 'JOB DONE' "$outfile" 2>/dev/null && continue
        mpirun -np "$NP" pw.x < "$infile" > "$outfile" 2>&1
    done

    mkdir -p tmp/_ph0
    if grep -q 'JOB DONE' dfpt/ph.out 2>/dev/null; then
        echo "  DFPT: done"
    else
        echo "  DFPT: running..."
        mpirun -np "$NP" ph.x -input dfpt/ph.in > dfpt/ph.out 2>&1
        grep -q 'JOB DONE' dfpt/ph.out || { echo "  DFPT FAILED"; tail -10 dfpt/ph.out; cd "$BASE"; continue; }
        echo "  DFPT: done"
    fi

    if grep -q 'JOB DONE' dfpt/q2r.out 2>/dev/null; then
        echo "  q2r: done"
    else
        echo "  q2r: running..."
        q2r.x < dfpt/q2r.in > dfpt/q2r.out 2>&1
        grep -q 'JOB DONE' dfpt/q2r.out || { echo "  q2r FAILED"; cd "$BASE"; continue; }
        echo "  q2r: done"
    fi

    if grep -q 'JOB DONE' dfpt/matdyn_dos.out 2>/dev/null; then
        echo "  matdyn: done"
    else
        echo "  matdyn: running..."
        matdyn.x < dfpt/matdyn_dos.in > dfpt/matdyn_dos.out 2>&1
        grep -q 'JOB DONE' dfpt/matdyn_dos.out || { echo "  matdyn FAILED"; cd "$BASE"; continue; }
        echo "  matdyn: done"
    fi

    cd "$BASE"
    echo "  $vname: COMPLETE"
done

echo "===== ALL 9 VOLUMES DONE ====="
