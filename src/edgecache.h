#ifndef EDGECACHE_H
#define EDGECACHE_H

#include <raylib.h>

struct EdgeCacheEntry {
    int vertexIdx;
    Vector3 normal;
    EdgeCacheEntry() : vertexIdx(-1), normal({0,0,0}) {}
};

#endif // EDGECACHE_H
