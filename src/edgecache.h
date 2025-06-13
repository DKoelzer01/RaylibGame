#ifndef EDGECACHE_H
#define EDGECACHE_H

#include <raylib.h>

struct EdgeCacheEntry {
    int vertexIdx;
    EdgeCacheEntry() : vertexIdx(-1) {}
};

#endif // EDGECACHE_H
