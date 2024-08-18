#!/bin/sh

if [ -f "$BELLE_TOP_DIR/src/anal/tatami/example/example_toy.dat" ]; then
  test_tatami_jpsiks $BELLE_TOP_DIR/src/anal/tatami/example/example_toy.dat
elif [ -n "$MY_TOP_DIR"  ]; then
  if [ -f "$MY_TOP_DIR/src/tatami/example/example_toy.dat" ]; then
    test_tatami_jpsiks $MY_TOP_DIR/src/tatami/example/example_toy.dat
  fi
fi
