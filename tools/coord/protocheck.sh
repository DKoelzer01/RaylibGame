#!/bin/bash
# Print kNetProtocol on each given remote branch (fetches each first).
for b in "$@"; do
  git fetch -q origin "$b" 2>/dev/null || { echo "$b: (no branch)"; continue; }
  echo "$b: $(git show "origin/$b:src/engine/netproto.h" 2>/dev/null | grep -oE 'kNetProtocol = [0-9]+')"
done
