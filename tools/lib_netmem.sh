#!/bin/bash
# Shared by the multi-instance net tests (#837): a host and its clients on one machine each commit a few GB (a plain client
# ~5.9 GB, a light one ~2.9 GB, the host ~5.6 GB), and running out of commit kills one with std::bad_alloc mid-test.
#   NET_LIGHT_CLIENT  flags for test clients: no spawn preload ring, a small chunk cache (they still stream what they see)
#   net_wait_for_memory GB NAME   waits up to 60 s for GB of free commit; prints a SKIP line and returns 1 if it never comes
NET_LIGHT_CLIENT="--preload 0 --cache-verts 1"
net_free_commit_gb() {
    powershell -NoProfile -Command "[math]::Floor((Get-CimInstance Win32_OperatingSystem).FreeVirtualMemory / 1MB)" 2>/dev/null | tr -d '\r'
}
net_wait_for_memory() {
    local need=$1 name=$2 free t=0
    while :; do
        free=$(net_free_commit_gb); free=${free:-0}
        [ "$free" -ge "$need" ] && return 0
        [ $t -ge 60 ] && { echo "=== $name: SKIP: only ${free} GB of commit free (needs ${need} GB for its game instances; other runs on this machine?)"; return 1; }
        sleep 5; t=$((t + 5))
    done
}
